#!/usr/bin/python

import sys
import os
import itertools
import argparse

def CreateOptions():
    parser = argparse.ArgumentParser()
    fmt = parser.add_mutually_exclusive_group()
    parser.add_argument("-i", dest="input", help="The input sequence", required=True, metavar="seq.fasta")
    parser.add_argument("-o", dest="output", help="The file to write tab-delimited k-mer counts", required=True, metavar="counts.tsv")
    parser.add_argument("-k", help="The length of the oligonucleotides to count", required=True, type=int)
    fmt.add_argument("-f", dest="fa_flag", help="Specifying the input sequence is in FASTA format [DEFAULT]", action='store_true', default=True)
    fmt.add_argument("-s", dest="bare", help="Denote that the file contains a single sequence, not preceded by a header", action="store_true")
    parser.add_argument("-m", help="If included, contig-wise k-mer counts will be returned opposed to k-mer frequencies for the whole file", action="store_true", default=False)
    args = parser.parse_args()
    return args

def PatternToNumber(pattern):
    """Returns a number representing the position of pattern
    of length k in an array of all permutations of k-mers"""
    count = n = 0
    k = len(pattern)
    c = (k - 1)
    while (c >= 0):
        factor = (4 ** n)
        if pattern[c] == "A":
            count += (0 * factor)
        if pattern[c] == "C":
            count += (1 * factor)
        if pattern[c] == "G":
            count += (2 * factor)
        if pattern[c] == "T":
            count += (3 * factor)
        c = (c - 1)
        n += 1
    return count

#Algorithm from Coursera:
#def PatternToNumber(pattern):
#    k = len(pattern)
#    if ( k == 0 ):
#        return 0    
#    symbol = pattern[k - 1]
#    prefix = pattern[0:(k - 1)]
#    if symbol == "A":
#        count = 0
#    if symbol == "C":
#        count = 1
#    if symbol == "G":
#        count = 2
#    if symbol == "T":
#        count = 3
#    return 4 * PatternToNumber(prefix) + count

def NumberToPattern(freq_array, num):
    """Finds the k-mer at position num in array
    containing all permutations of all k-mers"""
    for pattern in freq_array:
        pos = PatternToNumber(pattern)
        if (pos == num):
            return pattern
        else:
            pass

def ComputingFrequencies(fasta, k):
    count_array = dict()
    freq_array = dict()

    kmers_list = itertools.product('ACGT', repeat=k)
    for kmers in kmers_list:
        kmer = "".join(kmers)
        barcode = PatternToNumber(kmer)
        freq_array[barcode] = 0
        count_array[kmer] = 0
    n = 0
    
    for header in fasta:
        while (k <= len(fasta[header])):
            oligo = fasta[header][n:k]
            freq_array[PatternToNumber(oligo)] += 1
            n += 1
            k += 1

    for kmer in count_array:
        num = PatternToNumber(kmer)
        count_array[kmer] = freq_array[num]
    return count_array
       
def format_fasta(fileIn):
    fasta = dict()
    line = fileIn.readline()
    header = ""
    while line:
        if (line[0] == '>'):
            header = line.strip()
            fasta[header] = ""
        else:
            fasta[header] += line.strip()
        line = fileIn.readline()
    return fasta

#def FindingFrequentWordsBySorting(fasta, k):
#    FrequentPatterns = list()
#    index = list()
#    count = list()
#    i = 0
#    n = k
#    for header in fasta:
#        while (n <= len(fasta[header])):
#            oligo = fasta[header][i:n]
#            index[i] = PatternToNumber(oligo)
#            count[i] = 1
#            n += 1
#            i += 1
#        SortedIndex = index.sort()
#        i = 0
#        for kmer in SortedIndex:
#            if (SortedIndex[i] == SortedIndex[i-1]):
#                count[i] = count[i-1] + 1
#        i = 0
#        maxCount = max(count)
#        while (i < (4**k)):
#            if (count[i] == maxCount):
#                pattern = NumberToPattern(SortedIndex[i], k)
#                FrequentPatterns.append(pattern)
#        return FrequentPatterns

def writeCounts(output, count_array):
    """Output the k-mer counts in a matrix format"""
    with open(output, 'w') as outF:
        header = sorted(count_array.keys())
        outF.write("\t" + "\t".join(header) + "\n")
        kmers = sorted(count_array[header[0]].keys())
        for i in kmers:
            line = [i]
            for j in header:
                line.append(count_array[j][i])
            line = map(str, line)
            outF.write("\t".join(line) + "\n")

def main():
    args = CreateOptions()
    fasta = dict()
    with open(args.input) as fileIn:
        if (args.bare == True):
            tmp = ""
            for line in fileIn:
                tmp += line.strip()
            fasta[">1"] = tmp
        else:
            fasta = format_fasta(fileIn)
    if args.m:
        count_array = dict()
        for header in fasta:
            tmp_fa = dict()
            tmp_fa[header] = fasta[header]
            freq_array = ComputingFrequencies(tmp_fa, args.k)
            count_array[header] = freq_array
        writeCounts(args.output, count_array) 
    else:
        count_array = ComputingFrequencies(fasta, args.k)
        with open(args.output, 'w') as outF:
            sort_kmers = sorted(count_array.keys())
            for kmer in sort_kmers:
                outF.write(kmer+"\t"+str(count_array[kmer])+"\n")
        
main()

