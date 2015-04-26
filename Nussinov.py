#!/usr/bin/python
# -*- coding: utf-8 -*-


#    <Nussinov Algorithm (RNA Secondary Structure Prediction) Commandline Programm>
#    Copyright (C) <2015>  <Mathias Kreuter>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Implementation of the Nussinov RNA Secondary structure predicting algorithm
# pairing is only allowed if i and j are at least 4 bases apart

# Usage: Nussinov.py <FASTAfile>

import sys
import resource
import argparse
from itertools import combinations, ifilter
import RNA, igraph

# this function parses the commandline input 
def CmdParser(argv):
    p = argparse.ArgumentParser(description='This programm computes the optimal (maximum basepair containing) secondary structure of RNA using the classical Nussinov Algorithm.')
    p.add_argument('FastaFile', type=str, help="The RNA sequence as a FASTA format file. The nucleotides can be in lowercase, and the use of T for U is permitted. The programm does not process special characters such as N or D. (Don't use them.) '#' \
    might be used to comment out a part of the sequence.")
    g = p.add_mutually_exclusive_group()
    g.add_argument('-a', '--all', action='store_true', help="When this argument is passed all structures which have the maximal base pair count will be calculated. This flag should only be set if the sequence is shorter than \
    25 nucleotides. OTHERWISE YOUR COMPUTER IS GOING TO CRASH. If the structure of the sequence is very ambiguous 15 nt might already become problematic. In all honesty DON'T USE THIS, if you don't understand why writing a recursive algorithm \
    in python is inherently problematic.")
    g.add_argument('-o', '--one', action='store_true', help="When this argument is passed only one structure is calculated, and it probably isn't the most likely one.")
    g.add_argument('-g', '--graphic', action='store_true', help="When this argument is passed a graph representation of the secondary structure is drawn and displayed. Due to the algorithm used to calculate the graph representation, it will \
    look different every time the programm is called. So retry until the graph looks fine. This does not work for sequences larger than 100 nt due to scaling isues.")
    args = p.parse_args()
    return vars(args) 

# this funciton creates the internal sequence representation
def InputParser(arg):
    file = open(arg, "r")
    
    flag = False
    seq = ""
    for line in file:
        if line[0] == ">":
            flag = True
        elif flag and line[0] in "AGCUT":
            seq += (line[:-1]).upper()
    Matrix = []
    for _ in range(len(seq)):
        Matrix.append(["X" for _ in (range(len(seq)))])
    
    return seq, Matrix                        
                      
def IsPair(i, j):
    concat = i + j
    if concat == "GC" or concat == "CG":
        return True
    elif concat == "UA" or concat == "AU" or concat == "TA"or concat == "AT":
        return True
    else:
        return False    

# constructs the Nussinov Matrix        
def Nussinov(Sequence, SolutionsMatrix):
    # Initialization of the Matrix
    if SolutionsMatrix[0][0] == "X":    
        for i in range(len(Sequence)):
            j = i + 1
            SolutionsMatrix[i][i] = 0
            if j < (len(Sequence)):
                if IsPair(Sequence[i], Sequence[j]) and (i + 4 <= j):
                    SolutionsMatrix[i][j] = 1
                else:
                    SolutionsMatrix[i][j] = 0
         
        return Nussinov(Sequence, SolutionsMatrix)
        
    # Return the Matrix if the Matrix is full
    elif SolutionsMatrix[0][len(Sequence) - 1] != "X":
        S = [SolutionsMatrix[0]]
        for i in range(1, len(Sequence)):
            currow = SolutionsMatrix[i]
            fliprow = []
            for char in currow:
                if char == "X":
                    fliprow.append(0)
                else:
                    fliprow.append(char)
            S.append(fliprow)
            
        return S
    
    # recursive main part
    else:
        sp = 0
        for field in SolutionsMatrix[0]:
            if field != "X":
                sp += 1
        #print sp
        for j in range(sp, len(Sequence)):
            i = j - sp
        
            #for j in range(sp, len(Sequence) - 1 - sp):
            
            # construct the options stack from which the maximum is filled
            # in the matrix field
            OptStack = []
            # i and j form a base pair / base pairs 1 + Sequence[i+1][j-1]
            if IsPair(Sequence[i], Sequence[j]) and (i + 4 <= j):
                OptStack.append(1 + SolutionsMatrix[i + 1][j - 1])
            else:
                OptStack.append(0 + SolutionsMatrix[i + 1][j - 1])
            # i is unpaired / base pairs same as for Sequence[i + 1][j]
            OptStack.append(SolutionsMatrix[i + 1][j])
            # j is unpaired / base pairs same as for Sequence[i][j - 1]
            OptStack.append(SolutionsMatrix[i][j - 1])
            # both i and j are paired but with k1 and k2 respectively
            # / base pairs of S[i][k] + S[k+1][j]
            jkOptStack = [] 
            ikOptStack = []
            for jk in range(i+1, j): 
                jkOptStack.append(SolutionsMatrix[i][jk])
            for ik in range(i+2, j+1):
                ikOptStack.append(SolutionsMatrix[ik][j])
            kOptStack = [jkOptStack[x]+ikOptStack[x] for x in xrange(len(jkOptStack))]
            OptStack.append(max(kOptStack))
            # chosing the base pair maximizing solution
            SolutionsMatrix[i][j] = max(OptStack)
        
        return Nussinov(Sequence, SolutionsMatrix)

# this function traces back and returns all basepairings that are possible in the matrix
# this function might reach the recursion limit of python
# assumes that there is pairing!!!!!!!!!!
# before calling the function it has to be checked whether 
# Mat[0][len(Seq) -1] != 0
def cTraceBack(Seq, Mat, Pairs = [], TBStack = []):
    # after initizialistation the first item in the TBStack is poped in every 
    # recursion
    toolazytofixthisflag = True
    if TBStack != []:
        Ind = TBStack[0]
        TBStack = TBStack[1:]
        toolazytofixthisflag = False
    # initializing     
    if TBStack == [] and Pairs == [] and toolazytofixthisflag:
        i = 0
        j = len(Seq) - 1
        Ind = (i, j)
        TBStack.append(Ind)
        return cTraceBack(Seq, Mat, Pairs, TBStack)

    
    # return the pairs list if the traceback finished
    if TBStack == [] and Ind[0] >= Ind[1]: 
        return Pairs, Mat[0][len(Seq) - 1]
        
        
    # dismiss a finished Traceback Index Pair:
    if TBStack != [] and Ind[0] >= Ind[1]: 
        return cTraceBack(Seq, Mat, Pairs, TBStack)
    
    ## recursive main part:
    # i unpaired
    if Mat[Ind[0]][Ind[1]] == Mat[Ind[0] + 1][Ind[1]]:
        nInd = ((Ind[0] + 1), Ind[1])
        if nInd not in TBStack:
            TBStack = [nInd] + TBStack[:]
    # j unpaired
    if Mat[Ind[0]][Ind[1]] == Mat[Ind[0]][Ind[1] - 1]:
        nInd = (Ind[0], (Ind[1] - 1))
        if nInd not in TBStack:
            TBStack = [nInd] + TBStack[:]
    # base pair (i,j)
    if Mat[Ind[0]][Ind[1]] == (Mat[Ind[0] + 1][Ind[1] - 1] + 1) and IsPair(Seq[Ind[0]], Seq[Ind[1]]) and (Ind[0] + 4 <= Ind[1]):
        if Ind not in Pairs:
            Pairs.append(Ind)
        nInd = ((Ind[0] + 1), (Ind[1] - 1))
        if nInd not in TBStack:
            TBStack = [nInd] + TBStack[:]
    # split search
    else:
        for k in range((Ind[0] + 1), Ind[1]):
            if Mat[Ind[0]][Ind[1]] == (Mat[Ind[0]][k] + Mat[k + 1][Ind[1]]):
                if Ind[0] < k and (k+1) < Ind[1]:
                    inInd = ((Ind[0]), k)
                    jnInd = ((k + 1), (Ind[1]))
                    if jnInd not in TBStack:
                        TBStack = [jnInd] + TBStack[:]
                    if inInd not in TBStack:
                        TBStack = [inInd] + TBStack[:]
    
    # this clause is needed!!!!!!!!!!!!!!!
    if TBStack == []:
            return sorted(Pairs), Mat[0][len(Seq) - 1]
                    
    return cTraceBack(Seq, Mat, Pairs, TBStack)

# this function returns only one solution    
def TraceBack(Seq, Mat, Pairs = [], TBStack = []):
    # after initizialistation the first item in the TBStack is poped in every 
    # recursion
    toolazytofixthisflag = True
    if TBStack != []:
        Ind = TBStack[0]
        TBStack = TBStack[1:]
        toolazytofixthisflag = False
    # initializing     
    if TBStack == [] and Pairs == [] and toolazytofixthisflag:
        i = 0
        j = len(Seq) - 1
        Ind = (i, j)
        TBStack.append(Ind)
        return TraceBack(Seq, Mat, Pairs, TBStack)

    
    # return the pairs list if the traceback finished
    if TBStack == [] and Ind[0] >= Ind[1]: 
        return Pairs, Mat[0][len(Seq) - 1]
        
        
    # dismiss a finished Traceback Index Pair:
    if TBStack != [] and Ind[0] >= Ind[1]: 
        return TraceBack(Seq, Mat, Pairs, TBStack)
    
    ## recursive main part:
    # i unpaired
    elif Mat[Ind[0]][Ind[1]] == Mat[Ind[0] + 1][Ind[1]]:
        nInd = ((Ind[0] + 1), Ind[1])
        if nInd not in TBStack:
            TBStack = [nInd] + TBStack[:]
    # j unpaired
    elif Mat[Ind[0]][Ind[1]] == Mat[Ind[0]][Ind[1] - 1]:
        nInd = (Ind[0], (Ind[1] - 1))
        if nInd not in TBStack:
            TBStack = [nInd] + TBStack[:]
    # base pair (i,j)
    elif Mat[Ind[0]][Ind[1]] == (Mat[Ind[0] + 1][Ind[1] - 1] + 1) and IsPair(Seq[Ind[0]], Seq[Ind[1]]) and (Ind[0] + 4 <= Ind[1]):
        if Ind not in Pairs:
            Pairs.append(Ind)
        nInd = ((Ind[0] + 1), (Ind[1] - 1))
        if nInd not in TBStack:
            TBStack = [nInd] + TBStack[:]
    # split search
    else:
        for k in range((Ind[0] + 1), Ind[1]):
            if Mat[Ind[0]][Ind[1]] == (Mat[Ind[0]][k] + Mat[k + 1][Ind[1]]):
                if Ind[0] < k and (k+1) < Ind[1]:
                    inInd = ((Ind[0]), k)
                    jnInd = ((k + 1), (Ind[1]))
                    if jnInd not in TBStack:
                        TBStack = [jnInd] + TBStack[:]
                    if inInd not in TBStack:
                        TBStack = [inInd] + TBStack[:]
                    break
    
    # this clause is needed!!!!!!!!!!!!!!!
    if TBStack == []:
            return sorted(Pairs), Mat[0][len(Seq) - 1]
                    
    return TraceBack(Seq, Mat, Pairs, TBStack)

# very much improved function
# figures out the act. structures with maximal pairing            
def PosSecStruc(Pairs, PairNum):
    
    # discards structures in which one base pairs  more than one time
    def DuplicatesCheck(p):
        subcharL = []
        for char in p:
            for subchar in char:
                subcharL.append(subchar)
        for subchar in subcharL:        
            if subcharL.count(subchar) > 1:
                    return False

        return True
        
    #discardes structures in which the basepairs would overlap
    def OverlapCheck(ps):
        struc = sorted(ps)
        for i in xrange(len(struc)):
            for j in xrange(i +1,len(struc)):
                if struc[i][0] < struc[j][0] < struc[i][1] < struc[j][1]:
                    return False
                if struc[j][0] < struc[i][0] < struc[j][1] < struc[i][1]:
                    return False
        return True
            
    # this is very efficient code!!!!!! <--- REUSE        
    strucs = list(ifilter(lambda x: DuplicatesCheck(x), list(combinations(Pairs, PairNum))))
    out = list(ifilter(lambda x: OverlapCheck(x), strucs))
    return out
               
        
def main(argv):
    ArgumentDic = CmdParser(argv)
    if ArgumentDic["all"]:
        # necessary to allow the Traceback to finish without reaching the recursionlimit
        # only needed for calculating all possible structures
        resource.setrlimit(resource.RLIMIT_STACK, (2**29,-1))
        sys.setrecursionlimit(10**6)    
    s, m = InputParser(ArgumentDic["FastaFile"])
    m = Nussinov(s, m)
    # this catches sequences with no secondary structure
    if m[0][len(s) - 1] == 0:
        print 
        print "  Your Sequence has no computable secondary structure."
        print 
    else:
        if ArgumentDic["all"]:
            p, n = cTraceBack(s, m)
        else:
            p, n = TraceBack(s, m)
        print 
        print "  This sequence has a maximum of " + str(n) + " base pairs."
        print     
        print s
        R = PosSecStruc(p, n)
        for res in R:
            basepairseq = ""
            for _ in range(len(s)):
                basepairseq +="."
            for pair in res:
                bps = basepairseq[:pair[0]] + "(" + basepairseq[pair[0] + 1: pair[1]] + ")" + basepairseq[pair[1] + 1:]
                basepairseq = bps
            print basepairseq
        print
        if ArgumentDic["graphic"] == True:
            RNA.gmlRNA(s, basepairseq, "Nussinov.gml", "A")
            g = igraph.read("Nussinov.gml")            
            layout = g.layout("kk")
            igraph.plot(g, layout = layout)            
    
if __name__ == "__main__":
    main(sys.argv)

