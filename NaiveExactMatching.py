# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 20:52:05 2020

@author: Kenneth
"""

"""
This function turns the content of a FASTA file into a string.
"""
def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

"""
This function turns the content of a FASTQ file into two parallel lists.
One list, sequences, contains the reads in the file.
The other list, qualities, contains the base quality scores of the reads.
"""
def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

"""
This function turns a string representing a DNA sequence into its reverse complement.
"""
def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

"""
This function returns the offset of each alignment where p matches t.
p is a pattern of DNA sequence.
t is a longer DNA sequence than p, such as a genome.
The method encoded in the function is the naive exact matching algorithm.
"""
def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

"""
This is a modified version of the naive exact matching algorithm.
In addition to aligning p with t, it also aligns the reverse complement of p with t.
When p is its own reverse complement, it reverts back to the unmodified algorithm.
"""
def naive_with_rc(p, t):
    occurrences = []
    p_rc = reverseComplement(p)
    
    if p == p_rc:
        for i in range(len(t) - len(p) + 1):  # loop over alignments
            match = True
            for j in range(len(p)):  # loop over characters
                if t[i+j] != p[j]:  # compare characters
                    match = False
                    break
            if match:
                occurrences.append(i)  # all chars matched; record        
    else:   
        for i in range(len(t) - len(p) + 1):  # loop over alignments
            match = True
            for j in range(len(p)):  # loop over characters
                if t[i+j] != p[j]:  # compare characters
                    match = False
                    break
            if match:
                occurrences.append(i)  # all chars matched; record
            match = True
            for j in range(len(p_rc)):  # loop over characters
                if t[i+j] != p_rc[j]:  # compare characters
                    match = False
                    break        
            if match:
                occurrences.append(i)  # all chars matched; record
    return occurrences

"""
This is a modified version of the naive exact matching algorithm.
It accepts alignments with up to two mismatches.
Therefore, it identifies approximate rather than exact matches.
"""
def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        mm = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                mm += 1
                if mm > 2:
                    match = False
                    break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

"""
This function converts an ASCII character to a base quality score.
"""
def phred33ToQ(qual):
    return ord(qual) - 33

"""
How many times does a DNA sequence or its reverse complement occur in a genome?
"""
phix_genome = readGenome('lambda_virus.fa')
occurrences = naive_with_rc('AGGT', phix_genome)
len(occurrences)
occurrences = naive_with_rc('TTAA', phix_genome)
len(occurrences)

"""
What is the leftmost offset where a DNA sequence or its reverse complement occurs in a genome?
"""
phix_genome = readGenome('lambda_virus.fa')
occurrences = naive_with_rc('ACTAAGT', phix_genome)
min(occurrences)
occurrences = naive_with_rc('AGTCGA', phix_genome)
min(occurrences)

"""
How many times does a DNA sequence occur in a genome when up to two mismatches are allowed?
"""
phix_genome = readGenome('lambda_virus.fa')
occurrences = naive_2mm('TTCAAGCC', phix_genome)
len(occurrences)

"""
What is the leftmost offset where a DNA sequence or its reverse complement occurs in a genome when up to two mismatches are allowed?
"""
phix_genome = readGenome('lambda_virus.fa')
occurrences = naive_2mm('AGGAGGTT', phix_genome)
min(occurrences)

"""
The following code executes a series of tasks:
1. Extract the reads and the base quality scores of each read from a FASTQ file.
2. Calculate the number of reads and the length of each read.
3. Calculate the average quality score at each offset, i.e. in each sequencing cycle.
4. Identify the offset with the lowest average quality score.
"""
seqs, quals = readFastq('ERR037900_1.first1000.fastq')
read_count = len(quals)
read_length = len(quals[0])
qual_average = [0]*read_length
for qual in quals:
    for i in range(read_length):
        qual_average[i] += phred33ToQ(qual[i])
qual_average = [i/read_count for i in qual_average]
for i in range(read_length):
    if qual_average[i] == min(qual_average):
        print(i)