#!/usr/bin/env python3

import sys
import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
import argparse

__author__ = "Connor Morgan-Lang"
__maintainer__ = "Connor Morgan-Lang"
__date__ = "03-07-2020"
__version__ = "0.0.1"


def set_arguments():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-g", "--genbank", required=True,
                        help="A genbank file")
    parser.add_argument("-o", "--output_dir", dest="output", required=True,
                        help="The directory to write output files")

    return parser


def format_feature(feature: SeqFeature) -> str:
    feature_entry = "\t".join([str(feature.location.start), str(feature.location.end), feature.type]) + "\n"
    for qual_key, qual_val in feature.qualifiers.items():
        feature_entry += "\t\t\t" + qual_key + "\t" + qual_val.pop() + "\n"
    return feature_entry


def feature_table_header_gen(seq_record: SeqRecord) -> str:
    return "\t".join([">Feature", seq_record.name]) + "\n"


def genbank_operator(sys_args):
    arg_parser = set_arguments()
    args = arg_parser.parse_args(sys_args)

    if not os.path.isdir(args.output):
        os.mkdir(args.output)

    with open(args.genbank, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):  # type: SeqRecord
            feature_tab = os.path.join(args.output, record.name + ".tbl")
            with open(feature_tab, 'w') as feat_handle:
                feat_handle.write(feature_table_header_gen(record))
                for feature in record.features:  # type: SeqFeature
                    feat_handle.write(format_feature(feature))

    return


if __name__ == "__main__":
    genbank_operator(sys.argv[1:])
