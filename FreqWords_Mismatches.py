#!/usr/bin/python3

#Used to solve the Frequenct Words with Mismatches Problem
#kmer_counter.py and freq_pattern_hamming.py must be in the same directory!

import sys
import itertools
from argparse import ArgumentParser

from kmer_counter import PatternToNumber, NumberToPattern, ComputingFrequencies
from freq_pattern_hamming import HammingDistance

def ParseOptions():
    parser = ArgumentParser()
    parser.add_argument("-i", dest="input", help="The input sequence", required=True, metavar="seq.txt")
    parser.add_argument("-k", help="The length of the oligonucleotides to count", required=True, type=int)
    parser.add_argument("-d", dest="dist", help="The number of acceptable mismatches between k-mers", required=True, type=int)
    parser.add_argument("-C", help="Indicating the most frequent words should include reverse complements as well", default=False, action='store_true')
    args = parser.parse_args()
    return args

def revComp(string):
    rString = string[::-1]
    RC = ""
    for n in rString:
        if (n == 'A'):
            RC += 'T'
        if (n == 'T'):
            RC += 'A'
        if (n == 'G'):
            RC += 'C'
        if (n == 'C'):
            RC += 'G'
    return RC

def AddMismatchFrequencies(freq_array, hamming_dict, dist, complement):
    if (complement == True):
         for fa_kmer, fa_freq in freq_array.items():
            if (fa_freq == 0):
                pass
            else:
                for hd_kmer in hamming_dict:
                    RC = revComp(hd_kmer)
                    if (HammingDistance(fa_kmer, hd_kmer) <= dist) and (fa_kmer != hd_kmer):
                        hamming_dict[hd_kmer] += fa_freq
                        hamming_dict[RC] += fa_freq
                    else:
                        pass
    else:
        for fa_kmer, fa_freq in freq_array.items():
            if (fa_freq > 0):
                for hd_kmer in hamming_dict:
                    if (HammingDistance(fa_kmer, hd_kmer) <= dist) and (fa_kmer != hd_kmer):
                        hamming_dict[hd_kmer] += fa_freq
                    else:
                        pass
            else:
                pass

def findMax(hamming_dict):
    Maximum = max(hamming_dict.values())
    most_frequent = list()
    for kmer in hamming_dict:
        if (hamming_dict[kmer] == Maximum) and kmer not in most_frequent:
            most_frequent.append(kmer)
            RC = revComp(kmer)
            most_frequent.append(RC)
    return most_frequent

def main():
    Args = ParseOptions()
    fasta = dict()
    hamming_dict = dict()
    with open(Args.input) as fileIn:
        tmp = ""
        for line in fileIn:
            tmp += line.strip()
        fasta[">1"] = tmp
    print("Counting k-mers...")
    count_array = ComputingFrequencies(fasta, Args.k)
#Needs to return the array of PatternToNumber representations
#Only use the numbers that have a value > 0
    hamming_dict.update(count_array)
    print("Counting Mismatched k-mers...")
    AddMismatchFrequencies(count_array, hamming_dict, Args.dist, Args.C)
    print("Finding the most frequent k-mers in the string...")
    most_frequent = findMax(hamming_dict)
    print("\n", " ".join(most_frequent))

main()
