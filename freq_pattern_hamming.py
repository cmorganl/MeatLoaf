#!/usr/bin/python

import sys
import argparse

def getOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", dest="sequence", help="The file containing the hypothetical genome on a second line", required=True)
    parser.add_argument("-p", dest="pattern", help="The pattern to be found in the genome", required=True)
    parser.add_argument("-d", dest="distance", help="The hamming distance between the sequence and pattern to be considered a match", required=True, type=int)
    args = parser.parse_args()
    return args

def HammingDistance(p1, p2):
    i = d = 0
    while (i < len(p1)):
        if (p1[i] != p2[i]):
            d += 1
        else:
            pass
        i += 1
    return d

def ApproximatePatternCount(sequence, pattern, distance):
    i = 0
    count = 0
    k = len(pattern)
    while (i <= (len(sequence) - k)):
            oligo = sequence[i:(i+k)]
            if ( HammingDistance(oligo, pattern) <= distance ):
                count += 1
            else:
                pass
            i += 1
    return count

def main():
    pass

if __name__ == "__main__":
    args = getOptions()
    genome = ""
    with open(args.sequence) as data:
        genome = data.readline().strip()
    count = ApproximatePatternCount(genome, args.pattern, args.distance)
    print("The mer",args.pattern,"appears approximately", str(count),"times in the genome.")

main()
