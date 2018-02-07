#!/usr/bin/env python

import argparse

__author__ = 'Connor Morgan-Lang'


def set_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--fasta",
                        help="The fasta file to be subsetted.",
                        required=True)
    parser.add_argument("-o", "--output",
                        help="The output fasta file [DEFAULT = duprm.fasta]",
                        required=False,
                        default="duprm.fasta")
    args = parser.parse_args()
    return args


def read_fasta(fasta):
    seqs = dict()
    header = ""
    sequence = ""
    with open(fasta) as fas:
        line = fas.readline()
        while line:
            if line[0] == '>':
                if header:
                    seqs[header] = sequence
                header = line.strip()
                sequence = ""
            else:
                sequence += line.strip()
            line = fas.readline()
    seqs[header] = sequence
    return seqs


def write_unique(seqs, output):
    fasta_out = open(output, 'w')
    written = list()
    for header in seqs:
        if header not in written:
            written.append(header)
            fasta_out.write(header + "\n")
            fasta_out.write(seqs[header] + "\n")
    fasta_out.close()
    return


def main():
    args = set_arguments()
    seqs = read_fasta(args.fasta)
    write_unique(seqs, args.output)

main()
