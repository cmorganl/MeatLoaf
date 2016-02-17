#!/usr/bin/env python

import argparse as ap
import os

__author__ = 'Connor Morgan-Lang'


def get_sampleID(filename):
    basename = filename.split(os.sep)[-1]
    return '.'.join(basename.split('.')[:-1])


def getOptions():
    parser = ap.ArgumentParser()
    parser.add_argument("-t", "--otuTable",
                        help="seqs_otu.txt output from QIIME containing picked OTUs",
                        required=True)
    parser.add_argument("-o", "--output",
                        help="The output otu table in 'tidy' format (i.e., one row for each OTU-sequence pair)",
                        required=False,
                        default=".tsv")
    args = parser.parse_args()
    if args.output == ".tsv":
        args.output = get_sampleID(args.otuTable) + args.output
    return args


def load_otuTable(otuTable):
    otus = dict()
    with open(otuTable) as table:
        for line in table:
            fields = line.strip().split("\t")
            otus[fields[0]] = fields[1:]
    return otus


def tidy_table(OTUs, output):
    with open(output, 'w') as otus_out:
        for otu in OTUs:
            for seq in OTUs[otu]:
                otus_out.write(otu + "\t" + seq + "\n")
    return


def main():
    args = getOptions()
    OTUs = load_otuTable(args.otuTable)
    tidy_table(OTUs, args.output)

main()
