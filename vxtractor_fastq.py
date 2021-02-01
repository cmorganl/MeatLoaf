#!/usr/bin/env python3

import os
import sys
import glob
import argparse
import subprocess

from Bio.Data import IUPACData
import pyfastx

__description__ = "A tool for extracting variable regions from FASTQ files guided by V-Xtractor alignment positions."
__author__ = "Connor Morgan-Lang"


def _maketrans(complement_mapping):
    """Make a python string translation table (PRIVATE).
    Arguments:
     - complement_mapping - a dictionary such as ambiguous_dna_complement
       and ambiguous_rna_complement from Data.IUPACData.
    Returns a translation table (a string of length 256) for use with the
    python string's translate method to use in a (reverse) complement.
    Compatible with lower case and upper case sequences.
    For internal use only.
    """
    keys = "".join(complement_mapping.keys()).encode("ASCII")
    values = "".join(complement_mapping.values()).encode("ASCII")
    return bytes.maketrans(keys + keys.lower(), values + values.lower())


_dna_complement_table = _maketrans(IUPACData.ambiguous_dna_complement)


def reverse_complement(seq: str) -> str:
    return seq.translate(_dna_complement_table)[::-1]


class AmpliconSample:
    def __init__(self, sample_name: str, prefix: str, fwd: str, rev: str):
        self.name = sample_name
        self.fwd_primer = fwd
        self.rev_primer = rev
        self.fastq_name = prefix

        self.paired_end = True
        self.raw_fwd = ""
        self.raw_rev = ""
        self.trim_path_f = ""
        self.trim_path_r = ""
        self.merge_path = ""
        self.var_pos_tbl = ""
        self.final_fq = ""
        return

    def fetch_raw_fastq(self, input_dir):
        sample_fq = glob.glob(input_dir + self.name + "*")
        if len(sample_fq) > 2:
            raise AssertionError("More than two FASTQ files found for sample '{}'.\n".format(self.name))
        if len(sample_fq) == 1:
            self.paired_end = False
        return

    def rc_primers(self) -> (str, str):
        return reverse_complement(self.fwd_primer), reverse_complement(self.rev_primer)


def get_options(sys_args):
    parser = argparse.ArgumentParser(description="A tool for primer-trimming, merging and extracting 16S rRNA gene "
                                                 "variable regions from FASTQ files.",
                                     add_help=False)

    req_args = parser.add_argument_group("Required arguments")
    opt_args = parser.add_argument_group("Optional arguments")
    mis_args = parser.add_argument_group("Miscellaneous arguments")

    req_args.add_argument("-p", "--primer_map", dest="primers", required=True,
                          help="Path to a CSV file listing the forward and reverse primers used for each sample.")
    req_args.add_argument("-d", "--fastq_path", dest="fastq_dir", required=True,
                          help="Path to the directory containing FASTQ files.")

    opt_args.add_argument('-o', '--output_dir', default='./out', required=False,
                          help="Path to a directory to write the merged, trimmed FASTQ files [ DEFAULT = './out' ]")
    opt_args.add_argument('-o', '--tmp_dir', default='./tmp', required=False,
                          help="Path to a directory containing intermediate files [ DEFAULT = './tmp' ]")

    mis_args.add_argument("-t", "--threads", default=4,
                          help="The number of threads available for PEAR and cutadapt.")
    mis_args.add_argument('-v', '--verbose', action='store_true', default=False,
                          help='prints a more verbose runtime log')
    mis_args.add_argument("-h", "--help",
                          action="help", help="show this help message and exit")

    args = parser.parse_args(sys_args)
    return args


def prep_for_analysis(input_dir, output_dir, temporary_dir) -> None:
    if not os.path.isdir(input_dir):
        raise IOError("Unable to find directory containing FASTQ files '{}'".format(input_dir))

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    if not os.path.isdir(temporary_dir):
        os.mkdir(temporary_dir)
    return


def read_sample_primer_table(primer_table: str) -> list:
    amplicon_samples = []
    with open(primer_table) as sample_primers:
        for line in sample_primers:
            sample_name, prefix, fwd, rev = line.strip().split(',')
            sample = AmpliconSample(sample_name, prefix, fwd, rev)
            amplicon_samples.append(sample)
    return amplicon_samples


def cutadapt_wrapper(amplicon_samples: list, temporary_dir: str, threads: int) -> None:
    for sample in amplicon_samples:  # type: AmpliconSample
        sample.trim_path_f = os.path.join(temporary_dir, sample.fastq_name + "_cutadapt_R1.fastq")
        sample.trim_path_r = os.path.join(temporary_dir, sample.fastq_name + "_cutadapt_R2.fastq")
        rc_fwd, rc_rev = sample.rc_primers()
        cut_proc = subprocess.Popen(["cutadapt",
                                     "-n", str(4),
                                     "-g", sample.fwd_primer, "-a", rc_rev,
                                     "-G", sample.rev_primer, "-A", rc_fwd,
                                     "-o", sample.trim_path_f,
                                     "-p", sample.trim_path_r,
                                     sample.raw_fwd, sample.raw_rev,
                                     "--quiet", "-j", str(threads),
                                     "--trim-n", "--minimum-length", str(50)],
                                    shell=False)
        cut_proc.wait()
    return


def bbmerge_wrapper(amplicon_samples: list, temporary_dir: str, threads: int) -> None:
    for sample in amplicon_samples:  # type: AmpliconSample
        if not sample.paired_end:
            continue
        sample.merge_path = os.path.join(temporary_dir, sample.fastq_name + "_merged.fastq")

        # Find the adapter sequences
        adapt_proc = subprocess.Popen(["bbmerge.sh", "t=" + str(threads),
                                       "in1=" + sample.trim_path_f,
                                       "in2=" + sample.trim_path_r,
                                       "outa=adapters.fa"])
        adapt_proc.wait()

        # Merge the FASTQ files
        merge_proc = subprocess.Popen(["bbmerge.sh", "t=" + str(threads),
                                       "in1=" + sample.trim_path_f,
                                       "in2=" + sample.trim_path_r,
                                       "out=" + sample.merge_path,
                                       "adapters=adapters.fa"])
        merge_proc.wait()

        os.remove(sample.trim_path_f)
        os.remove(sample.trim_path_r)
    return


def fq2fa(fastq: str, fa: str) -> None:

    fa_handler = open(fa, 'w')

    fq = pyfastx.Fastq(fastq, build_index=False, full_name=True)
    for name, seq, qual in fq:
        fa_handler.write(">{}\n{}\n".format(name, seq))

    fa_handler.close()

    return


def vxtractor_wrapper(amplicon_samples: list, hmms_path: str, temporary_dir: str) -> None:
    for sample in amplicon_samples:  # type: AmpliconSample
        # Covert fastq file to FASTA
        tmp_fa = os.path.join(temporary_dir, sample.name + ".fa")
        fq2fa(fastq=sample.merge_path, fa=tmp_fa)

        # Run V-Xtractor
        sample.var_pos_tbl = os.path.join(temporary_dir, sample.name + "_vxtractor.csv")
        vx_proc = subprocess.Popen(["vxtractor.pl",
                                    "-r", "V3.-.V5",
                                    "-h", hmms_path,
                                    "-c", sample.var_pos_tbl,
                                    tmp_fa])
        vx_proc.wait()

        os.remove(tmp_fa)

    return


def read_vxtracted(vxtract_csv: str, start_name: str, end_name: str) -> dict:
    variable_positions = {}
    start_field_pos = 0
    end_field_pos = 0
    with open(vxtract_csv, 'r'):
        for line in vxtract_csv:
            if start_name in line and end_name in line:
                # Find the field positions for the start and end data in the CSV
                fields = line.strip().split()
                pos_dict = dict(zip(fields, range(len(fields))))
                start_field_pos = pos_dict[start_name]
                end_field_pos = pos_dict[end_name]
            elif line[0] == '#':
                continue
            fields = line.strip().split(',')
            variable_positions[fields[0]] = (fields[start_field_pos], fields[end_field_pos])

    return variable_positions


def trim_fq(fq_in: str, trim_positions: dict, fq_out: str) -> None:
    fq_in_handler = pyfastx.Fastq(fq_in, build_index=False, full_name=True)
    fq_out_handler = open(fq_out, 'w')

    for name, seq, qual in fq_in_handler:
        start, end = trim_positions[name]
        fq_out_handler.write("@{}\n{}\n+\n{}\n".format(name,
                                                       seq[start:end],
                                                       qual[start:end]))
    fq_out_handler.close()
    return


def main(cli_args):
    args = get_options(cli_args)
    prep_for_analysis(args.fastq_dir, args.output_dir, args.tmp_dir)
    # Read the file mapping sample to primers
    amplicon_samples = read_sample_primer_table(args.primers)
    # Find the raw FASTQ files
    for sample in amplicon_samples:
        sample.fetch_raw_fastq(args.fastq_dir)

    # Run cutadapt
    cutadapt_wrapper(amplicon_samples, args.tmp_dir, args.threads)
    # Merge reads with BBMerge
    bbmerge_wrapper(amplicon_samples, args.tmp_dir, args.threads)
    # Run V-Extractor on the merged reads in FASTA format
    vxtractor_wrapper(amplicon_samples, args.hmm_path, args.tmp_dir)

    for sample in amplicon_samples:  # type: AmpliconSample
        # Read the V-Xtractor alignment positions
        v_positions = read_vxtracted(sample.var_pos_tbl, start_name="V3rightlong", end_name="V5leftlong")
        # Stream the fastq file, and write the extracted positions
        sample.final_fq = os.path.join(args.output_dir, sample.name + "_extracted.fq")
        trim_fq(sample.merge_path, v_positions, sample.final_fq)
    return


if __name__ == '__main__':
    main(sys.argv)
