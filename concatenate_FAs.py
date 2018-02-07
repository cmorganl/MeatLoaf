#!/usr/bin/python

import os
import argparse


def get_options():
    parser = argparse.ArgumentParser(description="Script to concatenate multiple FASTA files and ensure no redundant"
                                                 "headers exist by prepending file names on the header.")
    parser.add_argument("-o", "--output", required=True,
                        help="The output FASTA file")
    parser.add_argument("-l", "--fasta_list", required=False, default=None,
                        help="A list of FASTA files to concatenate with paths included")
    args = parser.parse_args()

    if not os.path.isfile(args.fasta_list):
        raise IOError("ERROR: " + args.fasta_list + " doesn't exist!")

    return args


def write_data(fasta, sample, output_fasta):
    with open(fasta.strip()) as FA_in:
        for line in FA_in:
            if line.startswith('>'):
                new = ">" + sample + "__" + line[1:]
                output_fasta.write(new)
            else:
                output_fasta.write(line)
                    

def get_sample_name_from_file_name(file_name):
    base = os.path.basename(file_name)
    fragments = base.split('.')[:-1]
    sample = '.'.join(fragments)
    return sample


def main():
    args = get_options()
    if os.path.exists(args.output):
        response = raw_input("Output file already exists. Should it be removed? [y|n] ")
        if response == 'y':
            os.remove(args.output)
    try:
        output_fasta = open(args.output, 'w')
    except:
        raise IOError("Unable to open " + args.output + " for writing!")
    with open(args.fasta_list) as fastas:
        for line in fastas:
            sample = get_sample_name_from_file_name(line)
            write_data(line, sample, output_fasta)

main()
