#!/usr/bin/python
# -*- coding: utf-8 -*-

#    <Programm shuffling Biological Sequences in FASTA format using various algorithms>
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

### Usage examples: python SequenceShufflerMathias.py -k HH9.fasta 3 10
###                 python SequenceShufflerMathias.py -d HH9.fasta 10
###                 python SequenceShufflerMathias.py -s HH9.fasta 10
###                 python SequenceShufflerMathias.py -t 100

from __future__ import division
import sys, random, itertools, argparse

def CmdParser(argv):
    """
    this function defines the commandline interface
    """
    p = argparse.ArgumentParser(description='This programm shuffles (Protein and Nucleotide) Sequences using various Algorithms.')

    p.add_argument('FastaFile', type=str, help="The sequence in FASTA format. Only one sequence FASTA files can be used. The actual sequence can contain any ASCII characterset, upper and lower case however will be ignored.", default=False, nargs='?')
    p.add_argument('k', type=int, nargs='?', help="The size of the k-lets shuffled when the klet option is passed. The default value is 2.", default=2) 
    p.add_argument('num', type=int, nargs=1, help="The number of shuffled sequences that should be returned. Because of the nature of the algorithms employed, this value is only an upper limit. The program will return as many unique and differing from the input, sequences as it can generate. Especially in the case of the \
    k-let option this might mean that no output will be generated. Due to the nature of using randomness in computing it might be the case that in a subsequent run a different number or any output sequences can be generated, however if upon multiple trials no output is generated this probably means that the shuffling \
    can not be done without changing the frequencies the algorithm preserves. In the test case this variable spezifies the repeats of the subtests. (It's recommended to use 100)")       

    g = p.add_mutually_exclusive_group()
    g.add_argument('-t', '--test', action='store_true', help="When this argument is passed the programs algorithms will be tested.")
    g.add_argument('-d', '--dice', action='store_true', help="When this argument is passed the returned sequences will be created by randomly selecting the characters from the sequence alphabet with the probability they occur in the original sequence. This algorithm does only preserve the original character frequencies in the statistical mean.")
    g.add_argument('-s', '--swap', action='store_true', help="When this argument is passed the sequences returned will be generated by swapping all letters of the original sequence with another letter in a random fashion. By these means the character frequencies will be preserved in the output sequences. This is an implementation of the Fisher-Yates / Knuth shuffling algorithm.")
    g.add_argument('-k', '--klet', action='store_true', help="When this argument is passed the output will be generated by walking a uniformly random Euler path in the sequence (k - 1)-let multigraph. Uniform randomness is guaranteed, because the path is generated with a spanning tree (which is calculated using the Wilson algorithm). \
    The algorithm preserves the frequencies of subsequences from size one up to k.")
    
    args = p.parse_args()
    return vars(args) 

def InputParser(filename):
    """
    takes a (one sequence) FASTA file as input and returns the sequence as a string 
    and the sequencename as given in the FASTA header as a string
    """
    file = open(filename, "r")
    Truth = False
    seq = ""
    for line in file:
        if line[0] == ">":
            Truth = True
            seqname = line[1:-1]
        elif Truth:
            if line[-1:] == "\n":
                seq += (line[:-1]).upper()
            else:
                seq += line.upper()
    return seq, seqname

def charFreq(str):
    """
    helper function
    returns a dictionary with the frequencies of the charcters of the string
    """
    FreqDic = {}
    for char in str:
        if char not in FreqDic:
            FreqDic[char] = 1
        else:
            FreqDic[char] += 1       
    for key in FreqDic.keys():
        FreqDic[key] /= len(str)
    return FreqDic

def StatShuffler(seq, times): 
    """
    returns a list of sequences randomly generated with the frequencies of the character 
    Dictionary
    """
    seqlen = len(seq)
    FreqDic = charFreq(seq)
    res = []
    chars = FreqDic.keys()
    triggers = [FreqDic[x] for x in chars]
    for ita in range(1, len(triggers)):
        triggers[ita] += triggers[ita - 1]
    for _ in range(times):
        genseq = ""
        for _ in range(seqlen):            
            rnum = random.random()
            for ita in range(len(triggers)-1):
                if rnum >= triggers[-2]:
                    genseq += chars[-1]
                    break
                elif rnum < triggers[ita]:
                    genseq += chars[ita]
                    break
        res.append(genseq)        
        
    return res
    
def OneByOneExchangeShuffler(seq, times):
    """
    returns times "times" shuffeled sequences 
    shuffels by exchanging dinucleotides
    """
    res = []
    seqind = range(len(seq))
    for _ in range(times):
        shufL = list(seq)
        seqind_minI = seqind[:]
        for i in seqind[:-1]:
            j = random.choice(seqind_minI)
            seqind_minI.remove(i)
            chari = shufL[i]
            shufL[i] = shufL[j]
            shufL[j] = chari
        res.append(''.join(shufL))
    return res

   
def seq2Graph(seq, k):
    """
    returns the directed multigraph used in the klet algorithm
    only the part of the sequence which can be divided into k-1 sized
    vertices is used for the graph the rest is stored in rseq and also returned
    the verticeslist is unordered
    """
    if (len(seq) % (k-1)) != 0:
        cutoff = ((len(seq) // (k-1) ) * (k-1))
        rseq = seq[cutoff:]
        seq = seq[:cutoff]
    else:
        rseq = ""
    
    vertices = list(set([seq[i:i+k-1] for i in range(0, (len(seq)-k+2), k-1)]))
    adjacencylist = []
    for vertex in vertices:
        adjacencylist.append([])
    for ita in range(0, len(seq) - k+2, k-1):
        K_let = seq[ita:ita+k-1]
        nextK_let = seq[ita+k-1:ita+k+k-2]
        for i in range(len(vertices)):
            if vertices[i] == K_let:
                for j in range(len(vertices)):
                    if vertices[j] == nextK_let:
                        adjacencylist[i].append(j)
    return vertices, adjacencylist, seq, rseq

def rootIndexer(seq, vertices, k):
    """
    returns the root vertex index
    """
    root = seq[len(seq) - (k-1):]
    ct = -1
    for vertex in vertices:
        ct += 1
        if vertex == root:
            return ct
            
def treecrownIndexer(seq, vertices, k):
    """
    returns the starting node in the directed graph
    """
    treecrown = seq[:k-1]
    ct = -1
    for vertex in vertices:
        ct += 1
        if vertex == treecrown:
            return ct
          

def RandomSpanningTree(root, vertices, adjacencylist):
    """
    returns a dictionary containing the endnodes of the exitingedges of nodes
    in the graph
    """
    def RandomSuccesor(u):
        return random.choice(adjacencylist[u])
    Next = {}
    InTree = {}
    for i in range(len(vertices)):
        InTree[i] = False
    Next[root] = None
    InTree[root] = True
    for i in range(len(vertices)):
        u = i
        while not InTree[u]:
            Next[u] = RandomSuccesor(u)
            u = Next[u]
        u = i
        while not InTree[u]:
            InTree[u] = True
            u = Next[u]
    return Next

def WalkThroughTheWoodsOnASpringAfternoon(PathsNumber, PathLength, IdeaOfATree, root, treecrown, leaves, branches):
    """
    creates a random Eulerpath in the multigraph
    and walks it Pathsnumber of times
    --------------------------------------------
    branches = adjacencylist
    leaves = vertices
    IdeaOfATree = Next
    """
    def PavingTheWay(IdeaOfATree, branches):
        """
        constructs the Eulerpath
        """
        trail = []
        for b in range(len(branches)):
            if len(branches[b]) > 1:
                LeavingTwig = IdeaOfATree[b]
                Twigs = branches[b][:]
                if LeavingTwig != None:
                 Twigs.remove(LeavingTwig)
                random.shuffle(Twigs)
                trail.append(Twigs + [LeavingTwig])
            else:
                LeavingTwig = IdeaOfATree[b]
                if LeavingTwig == None:
                    trail.append(branches[b] + [LeavingTwig])
                else:
                    trail.append(branches[b])
        return trail

    def TrailWalker(trail, root, treecrown, PathLength):
        """
        walks the eulerpath
        """
        path = []
        CurrentPlace = treecrown
        CurrentWay = trail[CurrentPlace][0]
        trail[CurrentPlace] = trail[CurrentPlace][1:]
        while len(path) != PathLength:
            path.append(CurrentPlace)
            CurrentPlace = CurrentWay
            CurrentWay = trail[CurrentPlace][0]
            trail[CurrentPlace] = trail[CurrentPlace][1:]
            if CurrentWay == None:
                path.append(root)
                return path
        return path
    
    Wood = []
    for _ in range(PathsNumber):
        trail = PavingTheWay(IdeaOfATree, branches[:])
        path = TrailWalker(trail[:], root, treecrown, PathLength)
        flower = ""
        for sight in path:
            flower += leaves[sight]
        Wood.append(flower)
    return Wood
    
def K_letShuffler(seq, k, times):
    """
    shuffles the sequence "times" times, preserving the klet frequencies
    """
    vertices, adjacencylist, seq, rseq = seq2Graph(seq, k)
    root = rootIndexer(seq, vertices, k)
    treecrown = treecrownIndexer(seq, vertices, k)
    IdeaOfATree = RandomSpanningTree(root, vertices, adjacencylist)
    PathLength = len(seq) / (k-1)
    Wood = WalkThroughTheWoodsOnASpringAfternoon(times, PathLength, IdeaOfATree, root, treecrown, vertices, adjacencylist)
    for ita in range(times):
        Wood[ita] += rseq
    return Wood
    
def TestOneByOneExchangeShuffler(seq, testsize):
    """
    tests the OneByOneExchangeShuffler function
    by testing for the character frequencies of shuffled sequences
    (these should renain the same as in the original sequence), 
    their length and the recurance of sequences
    ---------------------------------------------------------------------------
    takes the sequence used for testing as input
    returns False if the test fails else True
    """
    FreqDic = charFreq(seq)
    testsize = 100
    
    # actual test
    testingresults = [OneByOneExchangeShuffler(seq, testsize) for _ in range(testsize)]
    
    testres1 = True
    # testing length of the shuffeled strings
    for partresult in testingresults:
        for result in partresult:
            if len(result) != len(seq):
                testres1 = False
    if testres1:
        print "    testing output length...                       successful"
    else:
        print "    testing output length...                     unsuccessful"
    
    testres2 = True
    # testing for repetition of a sequence, one is allowed
    flattenedList = list(itertools.chain.from_iterable(testingresults))
    if abs(len(flattenedList) - len(set(flattenedList))) > 1:
        testres2 = False
        
    if testres2:
        print "    testing uniqueness of output...                successful"
    else:
        print "    testing uniqueness of output...              unsuccessful"
    
    testres3 = True    
    # testing for deviation of the character frequencies from the ones 
    # for the testing sequence
    freqviolations = 0    
    for partresult in testingresults:
        partstring = ""
        for string in partresult:
            partstring += string
        partstringdic = charFreq(partstring)
        for key in partstringdic.keys():
            if FreqDic[key] != partstringdic[key]:
                freqviolations += 1
    if freqviolations != 0:
        testres3 = False
    
    if testres3:
        print "    testing character frequencies...               successful"
    else:
        print "    testing character frequencies...             unsuccessful"
    
    if testres1 and testres2 and testres3:
        return True
        
def TestKletShuffler(seq, Type, testsize):
    """
    tests the KletShuffler function
    by testing for the frequencies of nlets up to k of shuffled sequences
    (these should remain the same as in the original sequence), 
    their length and the recurance of sequences
    ---------------------------------------------------------------------------
    takes the sequence used for testing as input
    returns False if the test fails else True
    """
    
    testingresults = []    
    
    # actual test
    if Type == "nuc":
        testingresults.append([K_letShuffler(seq, 2, testsize) for _ in range(testsize)])
        testingresults.append([K_letShuffler(seq, 3, testsize) for _ in range(testsize)])
        testingresults.append([K_letShuffler(seq, 4, testsize) for _ in range(testsize)])
        
    if Type == "prot":
        testingresults.append([K_letShuffler(seq, 2, testsize) for _ in range(testsize)])
    
    testres1 = True
    # testing length of the shuffeled strings
    for partresult in testingresults:
        for subresults in partresult:
            for result in subresults:
                if len(result) != len(seq):
                    testres1 = False
    if testres1:
        print "    testing output length...                       successful"
    else:
        print "    testing output length...                     unsuccessful"
    
    testres2 = True
    # testing for repetition of a sequence, one is allowed
    for subtestingresults in testingresults: 
        
        flattenedList = list(itertools.chain.from_iterable(subtestingresults))
        if abs(len(flattenedList) - len(set(flattenedList))) > 1:
            testres2 = False
        
    if testres2:
        print "    testing uniqueness of output...                successful"
    else:
        print "    testing uniqueness of output...              unsuccessful"
    
    def kletfreqcalc(seq, k):
        """
        this function calculates the freq of klets
        returns a dictionary
        """
        FreqDic = {}
        alphabet = list(set([seq[i:i+k-1] for i in range(len(seq)-k+2)]))
        
        for char in alphabet:
            if char not in FreqDic:
                FreqDic[char] = 1
            else:
                FreqDic[char] += 1       
        for key in FreqDic.keys():
            FreqDic[key] /= len(seq)
        return FreqDic
    
    testres3 = True    
    # testing for deviation of the character frequencies from the ones 
    # for the testing sequence
    kct = 3
    for r in testingresults:
        kct += 1
        freqviolations = 0 
        for k in range(1, kct):
            seqdic = kletfreqcalc(seq, k)
            for partresult in r:
                for string in partresult:
                    strdic = kletfreqcalc(string, k)
                    for key in seqdic.keys():
                        if seqdic[key] != strdic[key]:
                            freqviolations += 1
        if freqviolations != 0:
            testres3 = False
                
    
    if testres3:
        print "    testing character frequencies...               successful"
    else:
        print "    testing character frequencies...             unsuccessful"

    if testres1 and testres2 and testres3:
        return True

def TestStatShuffler(seq, testsize):
    """
    tests the StatShuffler function
    by testing for the character frequencies of shuffled sequences, their length 
    and the recurance of sequences
    ---------------------------------------------------------------------------
    takes the sequence used for testing as input
    returns False if the test fails else True
    """
    epsilon = 0.05 #this is high so small sequences can be used for testing
    FreqDic = charFreq(seq)
    
    # actual test
    testingresults = [StatShuffler(seq, testsize) for _ in range(testsize)]
    
    testres1 = True
    
    # testing length of the shuffeled strings
    for partresult in testingresults:
        for result in partresult:
            if len(result) != len(seq):
                testres1 = False
    if testres1:
        print "    testing output length...                       successful"
    else:
        print "    testing output length...                     unsuccessful"
        
    testres2 = True
    # testing for repetition of a sequence, one is allowed
    flattenedList = list(itertools.chain.from_iterable(testingresults))
    if abs(len(flattenedList) - len(set(flattenedList))) > 1:
        testres2 = False
        
    if testres2:
        print "    testing uniqueness of output...                successful"
    else:
        print "    testing uniqueness of output...              unsuccessful"
    
    testres3 = True
    # testing for deviation by more than 5% of the character frequencies of the
    # testing sequence
    freqviolations = 0    
    for partresult in testingresults:
        partstring = ""
        for string in partresult:
            partstring += string
        partstringdic = charFreq(partstring)
        for key in partstringdic.keys():
            if abs(FreqDic[key] - partstringdic[key]) > epsilon:
                freqviolations += 1
    if freqviolations > ((testsize**2)*0.05):
        testres3 = False
    
    if testres3:
        print "    testing character frequency distributions...   successful"
    else:
        print "    testing character frequency distributions... unsuccessful"
        
    if testres1 and testres2 and testres3:
        return True

def GlobalTest(t):
    seqNuc = "ATTATTTTACTCTATACTCTATTTATTTCATATATACCATCTCATGTACTGTACCATCATAGTATGTCCTTAAATATATGTATATCGTGCATTAATGGCGTGCCCCATGCATATAGGCATGTACATATTATGCTTGATCTTACATGAGGACTTACATCTCAAAAGTTTATTTCAAGTGTATAGTCTGTAAGCATGTATTTCACTTAGTCCGGGAGCTTAATCACCAGGCCTCGAGAAACCAGCAACCCTTG"
    seqProt = "MQAKVENPLKSLRTAINRIVLVKLKDGSEYIGKLEQTDGTMNLVLRDCTEIREGTSEPVAKYGRVLIRGSNILFISVDYETVMNSEK"
    TestRes = []
    print 
    print " Testing the Dice algorithm."
    print 
    print "                    ______"
    print "        .-------.  /\\     \\"
    print "       /   o   /| /o \\  o  \\"
    print "      /_______/o|/   o\\_____\\"
    print "      | o     | |\\o   /o    /"
    print "      |   o   |o/ \\ o/  o  /"
    print "      |     o |/   \\/____o/"
    print "      '-------'"
    print
    print " Testing for Nucleotide Sequences."
    TestRes.append(TestStatShuffler(seqNuc,t))
    print " Testing for Protein Sequences."
    TestRes.append(TestStatShuffler(seqProt,t))
    print
    print " Testing the Fisher-Yates / Knuth shuffling algorithm."  
    print
    print "                   __"
    print "             _..-''--'----_."
    print "           ,''.-''| .---/ _`-._"
    print "         ,' \\ \\  ;| | ,/ / `-._`-."
    print "       ,' ,',\\ \\( | |// /,-._  / /"
    print "       ;.`. `,\\ \\`| |/ / |   )/ /"
    print "      / /`_`.\\_\\ \\| /_.-.'-''/ /"
    print "     / /_|_:.`. \\ |;'`..')  / /"
    print "     `-._`-._`.`.;`.\\  ,'  / /"
    print "         `-._`.`/    ,'-._/ /"
    print "           : `-/     \\`-.._/"
    print "           |  :      ;._ ("
    print "           :  |      \\  ` \\"
    print "            \\         \\   |"
    print "             :        :   ;"
    print "             |           /"
    print "             ;         ,'"
    print "            /         /"
    print "           /         /"
    print "                    /"
    print
    print " Testing for Nucleotide Sequences."
    TestRes.append(TestOneByOneExchangeShuffler(seqNuc,t))
    print " Testing for Protein Sequences."
    TestRes.append(TestOneByOneExchangeShuffler(seqProt,t))
    print
    print " Testing the k-let algorithm"
    print
    print "   +---------------------+"
    print "   v                     |"
    print " +----+     +----+     +----+       +----+"
    print " |    | --> | AA | --> |    | ----> | UA |"
    print " |    |     +----+     |    |       +----+"
    print " |    |                |    |         |"
    print " |    |                |    | ---+    |"
    print " | GG |                | AU |    |    |"
    print " |    | -------------> |    | <--+    |"
    print " |    |                |    |         |"
    print " |    |                |    |         |"
    print " |    | -------------> |    |         |"
    print " +----+                +----+         |"
    print "   ^                                  |"
    print "   +----------------------------------+"
    print
    print " Testing for Nucleotide Sequences."
    TestRes.append(TestKletShuffler(seqNuc, "nuc",t))
    print " Testing for Protein Sequences."
    TestRes.append(TestKletShuffler(seqProt, "prot",t))
    print
    Truth = True
    for Test in TestRes:
        if not Test:
            Truth = False
    if Truth:
        print " All tests were successful"
    else:
        print " At least one of the tests failed."
        print " You might want to run the tests again."
        print " Due to the nature of randomness, tests of algorithms employing "
        print " random number generators will sometimes fail."

def main(argv):
    ArgumentDic = CmdParser(argv)
    if ArgumentDic['FastaFile'] == False and ArgumentDic['test'] != True:
        sys.exit( "Please provide a FASTA file to be shuffled.")
    elif not ArgumentDic['test'] and ArgumentDic['FastaFile']:
        seq, seqname = InputParser(ArgumentDic['FastaFile'])

        
    if ArgumentDic['test']:
        GlobalTest(ArgumentDic['num'][0])
    elif ArgumentDic['dice']:
        resultingseq = set(StatShuffler(seq, ArgumentDic['num'][0]))
        res = []
        for s in resultingseq:
            if s != seq:
                res.append(s)
        output = ""
        ct = 0
        for resseq in res:
            ct += 1
            output += ">" + "Shuffled_" + str(ct) + " " + seqname + "\n" + resseq + "\n"
        print output
    elif ArgumentDic['swap']:
        resultingseq = set(OneByOneExchangeShuffler(seq, ArgumentDic['num'][0]))
        res = []
        for s in resultingseq:
            if s != seq:
                res.append(s)
        output = ""
        ct = 0
        for resseq in res:
            ct += 1
            output += ">" + "Shuffled_" + str(ct) + " " + seqname + "\n" + resseq + "\n"
        print output
        pass
    elif ArgumentDic['klet']:
        resultingseq = set(K_letShuffler(seq, ArgumentDic['k'], ArgumentDic['num'][0]))
        res = []
        for s in resultingseq:
            if s != seq:
                res.append(s)
        output = ""
        ct = 0
        for resseq in res:
            ct += 1
            output += ">" + "Shuffled_" + str(ct) + " " + seqname + "\n" + resseq + "\n"
        print output

               
if __name__ == "__main__":
    main(sys.argv)
