import argparse
import sys
import os

__author__ = 'Connor Morgan-Lang'


def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta",
                        help="FASTA file with split libraries (e.g., output from split_libraries.py).",
                        required=True)
    parser.add_argument("-m", "--map",
                        help="Mapping file - first three columns being library, barcode, adaptor - for QIIME pipeline",
                        required=True)
    parser.add_argument("-o", "--output_file",
                        help="Name of output FASTA file",
                        required=False,
                        default="seqs_trimmed.fna")
    args = parser.parse_args()
    return args


def digest_map(map_file):
    if not os.path.exists(map_file):
        sys.exit("ERROR: path to --map does not exist!")
    barcode_adapter_map = dict()
    with open(map_file) as fasting_map:
        line = fasting_map.readline()
        while line[0] == "#":
            line = fasting_map.readline()
        while line:
            fields = line.strip().split('\t')
            lib, barcode, adaptor = fields[:3]
            barcode_adapter_map[lib] = barcode + adaptor
            line = fasting_map.readline()
    return barcode_adapter_map


def trim_seqs(args, barcode_adapter_map):
    """
    Trims the barcode and adaptor sequences from the beginning of the input FASTA file.
    Assumes that the header is on a single line and the paired sequence is on the following line
    :param args:
    :param barcode_adapter_map:
    :return:
    """
    with open(args.fasta) as seqs_fna:
        with open(args.output_file, 'w') as seqs_trimmed:
            num_unmatched = 0
            line = seqs_fna.readline()
            lib_name = ""
            while line:
                # Random threshold of 1000 - could use proportion instead
                if num_unmatched >= 1000:
                    sys.exit("ERROR:\t" + lib_name + " does not match barcode and adaptor sequence! "
                             "Have these sequences already been trimmed?")
                if line[0] == '>':
                    if line[1:].split('_')[0] != lib_name:
                        lib_name = line[1:].split('_')[0]
                        num_unmatched = 0
                    seqs_trimmed.write(line)
                else:
                    # This must be the sequence line and the first n characters will be stripped
                    # where n is the length of barcode_adapter_map[lib_name]
                    # Here, we 'taste' the first few characters and make sure they match
                    if line[:7] == barcode_adapter_map[lib_name]:
                        seqs_trimmed.write(line[len(barcode_adapter_map[lib_name]):])
                    else:
                        num_unmatched += 1
                        print "Warning: library name doesn't match barcode!"
            	line = seqs_fna.readline()
    return

def main():
    args = get_options()
    barcode_adapter_map = digest_map(args.map)
    trim_seqs(args, barcode_adapter_map)

main()
