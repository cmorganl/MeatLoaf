#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'

import sys
import os
import re
import argparse
from itertools import chain


def get_arguments():
    parser = argparse.ArgumentParser(add_help=False)
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument("-C", "--collection",
                               help="Product of anvi-export-collection."
                                    "A tab-separated table with contig name in first column and bin name in second",
                               required=True)
    required_args.add_argument("-f", "--fasta",
                               help="FASTA file containing contigs to be extracted by collection",
                               required=True)

    optopt = parser.add_argument_group("Optional options")
    optopt.add_argument('-p', "--prefix", default="", required=False,
                        help="Prefix for the FASTA bins to be created")
    optopt.add_argument("-h", "--help",
                        action="help",
                        help="Show this help message and exit")

    miscellaneous_opts = parser.add_argument_group("Miscellaneous options")
    miscellaneous_opts.add_argument('-v', '--verbose', action='store_true', default=False,
                                    help='Prints a more verbose runtime log')

    args = parser.parse_args()

    if sys.version_info > (2, 9):
        args.py_version = 3
    else:
        args.py_version = 2

    return args


def load_anvi_collection(args):
    binning_map = dict()
    split_re = re.compile("(.*)_split_(.*)$")

    try:
        col_tab = open(args.collection, 'r')
    except IOError:
        raise IOError("ERROR: Unable to open collection " + args.collection + " for reading!")

    line = col_tab.readline().strip()
    while line:
        header = ""
        bin_id = ""
        if not line:
            # Continue if the line is blank
            pass
        else:
            try:
                header, bin_id = line.split('\t')
            except ValueError:
                if len(line.split('\t')) != 2:
                    sys.stderr.write("ERROR: Collection table is not in tab-delimited, two-column format!\n"
                                     "Assumption failed for line:\n" + line + "\n")
                else:
                    raise ValueError
            if header == "":
                sys.stderr.write("ERROR: header element is empty in line:\n" + line + "\n")
                sys.exit()
            if bin_id == "":
                sys.stderr.write("ERROR: bin element is empty in line:\n" + line + "\n")
                sys.exit()

            if split_re.match(header):
                split = int(split_re.match(header).group(2))
                if split > 1:
                    pass
                else:
                    header = split_re.match(header).group(1)
                    if bin_id not in binning_map.keys():
                        binning_map[bin_id] = list()
                    binning_map[bin_id].append(header)

        line = col_tab.readline().strip()

    if args.verbose:
        total = 0
        sys.stdout.write("# Bin_name\tNumber of Sequences\n")
        for collection in sorted(binning_map):
            total += len(binning_map[collection])
            sys.stdout.write("\t" + collection + "\t" + str(len(binning_map[collection])) + "\n")
        sys.stdout.write("# Total = " + str(total) + "\n")

    col_tab.close()
    return binning_map


def subset_fasta(args, headers):
    """
    Finds each header in `headers` within the `fasta` argument.
    NOTE: headers within the fasta file are split on spaces!
    :return: A dictionary with the contig headers of interest as keys and their respective sequences as their values.
    """
    subset = dict()
    header = ""
    add = 0
    try:
        fas = open(args.fasta, 'r')
    except IOError:
        raise IOError("ERROR: Unable to open fasta " + args.fasta + " for reading!\n")

    if args.verbose:
        sys.stdout.write("Extracting binned sequences from " + args.fasta + "... ")
        sys.stdout.flush()

    line = fas.readline()
    while line:
        if line[0] == '>':
            header = line[1:].strip()
            if header in headers:
                subset[header] = ""
                add = 1
            else:
                add = 0
        elif add == 1:
            subset[header] += line.strip()
        else:
            pass
        line = fas.readline()

    if args.verbose:
        sys.stdout.write("done.\n")

    if args.verbose:
        sys.stdout.write("\tNumber of sequences parsed from FASTA file: " + str(len(subset)) + "\n")
    if len(headers) != len(subset):
        sys.stderr.write("WARNING: Number of sequences binned (" + str(len(headers))
                         + ") is not equal to number of sequences parsed (" + str(len(subset)) + ")!\n")

    return subset


def write_bins(args, binning_map, binned_subset):
    if args.verbose:
        sys.stdout.write("\tWriting FASTA files for each bin.\n")
        sys.stdout.write("\tOutput files will be in: " + args.prefix + "_*.fasta" + "\n")

    for bin_id in sorted(binning_map):
        bin_fasta = args.prefix + "_" + bin_id + ".fasta"
        fa_string = ""
        try:
            output_handler = open(bin_fasta, 'w')
        except IOError:
            sys.stderr.write("ERROR: Unable to open " + bin_fasta + " for writing!\n")
            raise IOError
        binned_seqs = binning_map[bin_id]
        for header in binned_seqs:
            fa_string += '>' + header + "\n"
            fa_string += binned_subset[header] + "\n"

        output_handler.write(fa_string)
        output_handler.close()
    return


def main():
    args = get_arguments()
    binning_map = load_anvi_collection(args)
    headers = list(chain.from_iterable(binning_map.values()))
    binned_subset = subset_fasta(args, headers)
    write_bins(args, binning_map, binned_subset)


main()
