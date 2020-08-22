#!/usr/bin/env python

import os
import unittest
import sys
from time import time
import logging
import argparse

from pyfastx import Fasta, Fastq
from pyfastxcli import fastx_format_check


__MAX_BUFFER_LEN = 1E6


class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.fq = Fastq(file_name="test_data/test_TarA.1.fq", build_index=False)
        self.ont = Fastq(file_name="test_data/ZCS_1K.fq", build_index=False)
        self.output = open(file="test_data/tmp.txt", mode='w')
        return

    def tearDown(self) -> None:
        self.output.close()
        return

    def test_seq_chopper(self):
        seqs = seq_chopper("ACGTAGATCATC\n", 5)
        self.assertEqual(list, type(seqs))
        self.assertEqual(2, len(seqs))
        for seq in seqs:
            self.assertEqual(5, len(seq))
        return

    def test_level_fastq(self):
        ont_subseqs = level_fastq(self.ont, self.output, 1000)
        ilmn_subseqs = level_fastq(self.fq, self.output, 50)
        self.assertEqual(36, ilmn_subseqs)
        return

    def test_max_seq_level_fastq(self):
        ilmn_subseqs = level_fastq(self.fq, self.output, seq_len=50, num=30, num_type='s')
        self.assertEqual(30, ilmn_subseqs)

    def test_max_char_level_fastq(self):
        ilmn_subseqs = level_fastq(self.fq, self.output, seq_len=50, num=100, num_type='c')
        self.assertEqual(2, ilmn_subseqs)

    def test_prep_output(self):
        test_fh = prep_output(fastx_file="test_data/test_TarA.1.fq", output_fastx='', extension="fq")
        self.assertEqual("test_data/test_TarA.1_levelled.fq", test_fh.name)
        return


def get_options():
    parser = argparse.ArgumentParser(add_help=False)
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument("-f", "--fastx", dest="fastx",
                               help="Path to a fasta or fastq file with sequences to be levelled.",
                               required=True)

    optopt = parser.add_argument_group("Optional arguments")
    optopt.add_argument('-o', '--output', default='', required=False,
                        help='The output file [ DEFAULT = ./${fastx}_levelled.f[a|q] ]')
    optopt.add_argument("-l", "--read_length", default=1E3, required=False, type=int,
                        help="The minimum and maximum length to level the reads. [ DEFAULT = 1E3 ]")
    optopt.add_argument("-n", "--number", default=0, required=False, dest="num", type=int,
                        help="The amount of sequence information (e.g. contigs, reads, base-pairs) to write"
                             " [ DEFAULT = no limit ].")
    optopt.add_argument("-t", "--limit_type", required=False, default='s', choices=['s', 'c'], dest="lim_t",
                        help="The type to constrain the '--number' parameter by, either sequences (s) or characters (c)"
                             " [ DEFAULT = 's' ].")

    miscellaneous_opts = parser.add_argument_group("Miscellaneous options")
    miscellaneous_opts.add_argument('--overwrite', action='store_true', default=False,
                                    help='overwrites previously processed output folders')
    miscellaneous_opts.add_argument('-v', '--verbose', action='store_true', default=False,
                                    help='prints a more verbose runtime log')
    miscellaneous_opts.add_argument("-h", "--help",
                                    action="help",
                                    help="show this help message and exit")

    args = parser.parse_args()

    return args


def level_fasta(fa, out_fh, seq_len: int, num=0, num_type='s') -> int:
    buff = ""
    n_subseqs = 0
    n_char = 0
    for name, seq in fa:  # type: (str, str)
        if len(seq) < seq_len:
            continue
        else:
            sub_seqs = seq_chopper(seq, seq_len)
            i = 0
            while sub_seqs:
                subseq = sub_seqs.pop(0)
                buff += ">%s_%d\n%s\n" % (name, i, subseq)
                i += 1
                n_subseqs += 1
                n_char += len(subseq)

                if num > 0:
                    if (num_type == 's' and n_subseqs == num) or (num_type == 'c' and n_char >= num):
                        out_fh.write(buff)
                        return n_subseqs

        if len(buff) > __MAX_BUFFER_LEN:
            out_fh.write(buff)
            buff = ""
    out_fh.write(buff)
    return n_subseqs


def level_fastq(fq, out_fh, seq_len: int, num=0, num_type='s'):
    buff = ""
    n_subseqs = 0  # Number of sequences written
    n_char = 0  # The number of characters written
    for name, seq, qual in fq:
        if len(seq) < seq_len:
            continue
        else:
            sub_seqs = seq_chopper(seq, seq_len)
            sub_quals = seq_chopper(qual, seq_len)
            i = 0
            while sub_seqs:
                subseq = sub_seqs.pop(0)
                subqual = sub_quals.pop(0)
                buff += "@%s_%d\n%s\n+\n%s\n" % (name, i, subseq, subqual)
                i += 1
                n_subseqs += 1
                n_char += len(subseq)

                if num > 0:
                    if (num_type == 's' and n_subseqs == num) or (num_type == 'c' and n_char >= num):
                        out_fh.write(buff)
                        return n_subseqs

        if len(buff) > __MAX_BUFFER_LEN:
            out_fh.write(buff)
            buff = ""

    out_fh.write(buff)
    return n_subseqs


def seq_chopper(seq: str, seq_len: int) -> list:
    """
    Takes a string and chops it into substrings of length seq_len

    :param seq: A string of characters, e.g. 'ACCGATC"
    :param seq_len: The desired subsequence lengths
    :return: A list of subsequences
    """
    seq = seq.strip()
    i = 0
    max_idx = seq_len*(i+1)
    sub_seqs = []
    while max_idx <= len(seq):
        sub_seqs.append(seq[seq_len*i:max_idx])
        i += 1
        max_idx = seq_len*(i+1)
    return sub_seqs


def prep_output(fastx_file: str, output_fastx: str, extension: str, gzipped=False):
    if extension[0] != '.':
        extension = '.' + extension

    if len(output_fastx) == 0:
        file_name, suffix1 = os.path.splitext(fastx_file)

        if gzipped:
            file_name, suffix2 = os.path.splitext(file_name)
        output_fastx = file_name + "_levelled" + extension

    try:
        output_handler = open(output_fastx, 'w')
    except IOError:
        logging.error("Unable to open {} for writing.\n".format(output_fastx))
        sys.exit(3)
    return output_handler


def get_fastx(fastx_file) -> (str, str):
    fastx_type = fastx_format_check(fastx_file)
    if fastx_type == 'fasta':
        fx = Fasta(file_name=fastx_file, build_index=False)
        ext = '.fa'
    elif fastx_type == 'fastq':
        fx = Fastq(file_name=fastx_file, build_index=False)
        ext = '.fq'
    else:
        logging.error("Unknown fastx type: '{}'\n".format(fastx_type))
        sys.exit(3)
    return fx, ext


def main():
    args = get_options()
    fx, ext = get_fastx(args.fastx)
    fh = prep_output(fastx_file=args.fastx, output_fastx=args.output, extension=ext, gzipped=fx.is_gzip)

    start = time()
    if type(fx) is Fasta:
        level_fasta(fx, fh, args.read_length, args.num, args.lim_t)
    elif type(fx) is Fastq:
        level_fastq(fx, fh, args.read_length, args.num, args.lim_t)
    fh.close()
    end = time()
    logging.debug("{} completed levelling in {}s.\n".format(args.fastx, end - start))
    return


if __name__ == "__main__":
    main()
