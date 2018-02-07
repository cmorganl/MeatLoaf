#!/usr/bin/env python
__author__ = 'connor'

import argparse
import sys
import os


def set_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--list", help="The list of fasta files to be concatenated.", required=True)
    parser.add_argument("-o", "--output", help="The output fasta file [DEFAULT = concat.fasta]", required=False,
                        default="concat.fasta")
    parser.add_argument("-m", "--minLength", help="The minimum contig length [DEFAULT = 200].",
                        required=False, default=200, type=int)
    parser.add_argument("-s", "--scaffold", action="store_true", default=False,
                        help="Glue the contigs together into unordered scaffolds by writing a string of 100 'N's "
                             "between each contig")
    args = parser.parse_args()
    return args


def load_list(file_list):
    files = list()
    with open(file_list) as LoF:
        for line in LoF:
            if ( line[0] == '>'):
                files.append(line[1:].strip())
            else:
                files.append(line.strip())
    return files


def drop_ext_path(file):
    basename = file.split(os.sep)[-1]
    name = '.'.join(basename.split('.')[:-1])
    return name


def load_scaffold(open_fasta, minLength):
    """
    Glue the contigs together with 100-character strings of 'N's
    :param open_fasta: An opened file in FASTA-format
    :param minLength: Minimum length for a contig to be included in the output FASTA
    :return: A dictionary representation of the fasta file
    """
    scaffolds = dict()
    return scaffolds


def concatenate_fasta_files(files, output, minLength, scaffold):
    """
    Write the concatenated fasta file with headers for each file being the file name.
    :param files: List of FASTA files to concatenate
    :param output: The output file
    :param minLength: Minimum length for a contig to be included in the output FASTA
    """
    ordered_files = sorted(files)
    fa_out = open(output, 'w')
    head = ""
    seq = ""
    for file in ordered_files:
        acc = 1
        with open(file) as fas:
            line = fas.readline()
            while line:
                if line[0] == '>':
                    if len(seq) >= minLength:
                        if scaffold == True:
                            seq += ("N"*100)
                        else:
                            fa_out.write('>' + head + "\n")
                            fa_out.write(seq + "\n")
                            acc += 1
                            seq = ""
                    name = drop_ext_path(file)
                    head = name + '_' + str(acc)
                else:
                    seq += line.strip()
                line = fas.readline()
        if (len(seq) >= minLength):
            fa_out.write('>' + head + "\n")
            fa_out.write(seq + "\n")
        seq = ""
        head = ""
    fa_out.close()
    return 0


def main():
    print "Beginning", sys.argv[0]
    args = set_arguments()
    files = load_list(args.list)
    print "Writing the concatenated FASTA file to", args.output
    concatenate_fasta_files(files, args.output, args.minLength, args.scaffold)
    print "Finished."
    return 0

main()

