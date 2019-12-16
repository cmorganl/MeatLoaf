#!/usr/bin/env python
__author__ = 'connor'

import argparse
import logging
import sys
import os


def set_arguments():
    parser = argparse.ArgumentParser(description="Script to concatenate multiple FASTA files and ensure no redundant"
                                                 " headers exist by prepending file names on the header.")
    parser.add_argument("-l", "--list",
                        help="The list of fasta files to be concatenated.",
                        required=True)
    parser.add_argument("-o", "--output",
                        help="The output fasta file [DEFAULT = concat.fasta]",
                        required=False, default="concat.fasta")
    parser.add_argument("-m", "--min_length",
                        help="The minimum sequence length [DEFAULT = 200].",
                        required=False, default=200, type=int)
    parser.add_argument("-s", "--scaffold",
                        help="Glue the contigs together into unordered scaffolds by writing a string of 100 'N's "
                             "between each contig",
                        action="store_true", default=False)
    args = parser.parse_args()
    return args


def load_list(file_list):
    files = list()
    with open(file_list) as list_handler:
        for line in list_handler:
            files.append(line.strip())
    return files


def get_file_name_prefix(file):
    basename = file.split(os.sep)[-1]
    name = '.'.join(basename.split('.')[:-1])
    return name


# No bioinformatic software would be complete without a contribution from Heng Li.
# Adapted from his readfq generator
def generate_fasta(fasta_handler):  # this is a generator function
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:
            for line in fasta_handler:  # search for the start of the next record
                if line[0] == '>':  # fasta header line
                    last = line[:-1]  # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:], [], None
        for line in fasta_handler:  # read the sequence
            if line[0] == '>':
                last = line[:-1]
                break
            seqs.append(line[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs)  # yield a fasta record
            if not last:
                break
        else:
            seq, seqs = ''.join(seqs), []
            for line in fasta_handler:  # read the quality
                seqs.append(line[:-1])
            if last:  # reach EOF before reading enough quality
                yield name, seq  # yield a fasta record instead
                break


def read_fasta_to_dict(fasta_file):
    """
    Reads any fasta file using a generator function (generate_fasta) into a dictionary collection

    :param fasta_file: Path to a FASTA file to be read into a dict
    :return: Dict where headers/record names are keys and sequences are the values
    """
    fasta_dict = dict()
    try:
        fasta_handler = open(fasta_file, 'r')
    except IOError:
        logging.error("Unable to open " + fasta_file + " for reading!\n")
        sys.exit(5)
    for record in generate_fasta(fasta_handler):
        name, sequence = record
        fasta_dict[name] = sequence.upper()
    return fasta_dict


def concatenate_fasta_files(files, output, min_length, scaffold, sep='_') -> None:
    """
    Write the concatenated fasta file with headers for each file being the file name.

    :param files: List of FASTA files to concatenate
    :param output: The output file
    :param min_length: Minimum length for a contig to be included in the output FASTA
    :param scaffold: Boolean indicating whether contigs should be connected by string of 'N's
    :param sep: A character to be used for separating file prefixes from the original FASTA headers
    """
    ordered_files = sorted(files)
    glue = "N" * 100
    buffer = ""

    try:
        fa_out = open(output, 'w')
    except IOError:
        logging.error("Unable to open output FASTA file '%s' for reading")
        sys.exit(3)

    for fasta_file in ordered_files:
        file_prefix = get_file_name_prefix(fasta_file)
        fasta_dict = read_fasta_to_dict(fasta_file)
        if scaffold:
            fa_out.write(">" + file_prefix + "\n")
            fa_out.write(glue.join(list(fasta_dict.values())))
            continue
        for seq_name in fasta_dict:
            if len(fasta_dict[seq_name]) < min_length:
                continue

            buffer += '>' + file_prefix + sep + seq_name + "\n"
            buffer += fasta_dict[seq_name] + "\n"

            if len(buffer) > 1E6:
                fa_out.write(buffer)
                buffer = ""

        # Write and clear what is remaining in the buffer
        fa_out.write(buffer)
        buffer = ""
    fa_out.close()
    return


def main():
    logging.info("\t\t** Beginning %s **\n" % sys.argv[0])
    args = set_arguments()
    files = load_list(args.list)
    logging.info("Writing the concatenated FASTA file to '%s'... " % args.output)
    concatenate_fasta_files(files, args.output, args.min_length, args.scaffold)
    logging.info("done.\n")
    return


main()
