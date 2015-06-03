#!/usr/bin/env python
__author__ = "Connor Morgan-Lang"

import sys
import argparse

def set_arguments():
    parser = argparse.ArgumentParser()
    extract = parser.add_argument_group("extract")
    parser.add_argument("-i", "--fasta", help = "The fasta file to be subsetted.", required=True)
    parser.add_argument("-l", "--list", help="The list of contig headers to be subsetted from the fasta file", required=True)
    parser.add_argument("-o", "--output", help="The output fasta file [DEFAULT = subset.fasta]", required=False, default="subset.fasta")
    extract.add_argument("-s", "--slice", help="Flag indicating a specific sequence should be extracted from a specific header. \
                         The headers to be used for extraction will be pulled from the provided list", required=False, action="store_true")
    extract.add_argument("-p", "--pos", help="Positions (start,end) of the provided contig to extract", required=False)
    args = parser.parse_args()
    return args

def subset_fasta(fasta, headers):
    """
    :return: A dictionary with the contig headers of interest as keys and their respective sequences as their values.
    """
    subset = dict()
    header = ""
    add = 0
    with open(fasta) as fas:
        line = fas.readline()
        while line:
            if ( line[0] == '>' ):
                header = line.strip()
                if header in headers:
                    subset[header] = ""
                    add = 1
                else:
                    add = 0
            elif (add == 1):
                subset[header] += line.strip()
            else:
                pass
            line = fas.readline()
    return subset

def load_list(header_list):
    headers = list()
    with open(header_list) as LoH:
        for line in LoH:
            headers.append(line.strip())
    return headers

def extract_seq(contig_seq, pos):
    """
    :return: A sub-string of the contig sequence as determined by the pos tuple.
    """
    start, end = pos
    return str(contig_seq[start:end])


def main():
    args = set_arguments()
    headers = load_list(args.list)
    subset = subset_fasta(args.fasta, headers)
    if args.slice:
        with open(args.output, 'w') as fa_out:
            for contig in subset:
                seq = extract_seq(subset[contig], args.pos)
                fa_out.write('>' + str(contig) + '\n')
                fa_out.write(str(seq) + "\n")
    else:
        with open(args.output, 'w') as fa_out:
            for contig in subset:
                fa_out.write('>' + str(contig) + '\n')
                fa_out.write(subset[contig] + "\n")
    return 0

main()
