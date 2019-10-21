#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'

import argparse
import sys
import re


def get_options():
    parser = argparse.ArgumentParser(add_help=False)
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument("-a", "--alignment_table",
                               help="tabular (tsv-delimited) file output by fr-hit",
                               required=True)
    required_args.add_argument("-r", "--reference",
                               help="name of the reference",
                               required=True)
    required_args.add_argument("-q", "--query",
                               help="name of the query",
                               required=True)
    required_args.add_argument("-f", "--fasta", dest="genome",
                               help="reference genome in FASTA format",
                               required=True)
    required_args.add_argument("-n", "--num_reads",
                               help="number of reads in the query (meta)genome that were aligned to the reference",
                               required=True,
                               type=int)

    optopt = parser.add_argument_group("Optional options")
    optopt.add_argument('-g', "--group", default="contig", choices=["contig", "genome"], required=False,
                        help="Reports summary stats for individual contigs or entire reference.")
    optopt.add_argument('-o', '--output', default='./output/', required=False,
                        help='output directory [DEFAULT = ./output/]')
    optopt.add_argument('-i', '--identity',
                        help='the sequence identity threshold. Anything lower than this value is removed. '
                             '[ DEFAULT = 50 ]',
                        default=50,
                        type=int)
    optopt.add_argument('-p', "--alignment_proportion",
                        help="minimum proportion threshold of the query aligned. A proportion below this is removed. "
                             "[ DEFAULT = 0.5 ]",
                        default=0.5,
                        type=float)

    miscellaneous_opts = parser.add_argument_group("Miscellaneous options")
    miscellaneous_opts.add_argument('--overwrite', action='store_true', default=False,
                                    help='overwrites previously processed output folders')
    miscellaneous_opts.add_argument('-v', '--verbose', action='store_true', default=False,
                                    help='prints a more verbose runtime log')
    miscellaneous_opts.add_argument("-h", "--help",
                                    action="help",
                                    help="show this help message and exit")

    args = parser.parse_args()

    if args.alignment_proportion > 1:
        sys.stderr.write("ERROR: `--alignment_proportion` x must be 0 < x < 1!\n")
        sys.stderr.write("Exiting...\n")
        sys.exit()

    return args


class Mapping:
    def __init__(self):
        self.proportion = 0
        self.identity = 0
        self.strand = ""
        self.reference = ""

    def load_dat(self, fields):
        if len(fields) != 11:
            sys.stderr.write("ERROR: fr-hit output line is not of the expected format!\n")
            sys.exit(3)
        read_len = int(re.sub('nt|aa', '', fields[1]))
        self.proportion = float(int(fields[3])/read_len)
        self.identity = float(re.sub('%', '', fields[7]))
        self.reference = fields[8]
        self.strand = fields[6]
        return


class RefSeq:
    def __init__(self, header, sequence):
        self.name = header
        self.sequence = sequence
        self.length = len(sequence)
        self.num_f = 0
        self.num_r = 0
        self.fpkm = 0
        self.rpk = 0
        self.tpm = 0

    def calc_fpkm(self, num_reads):
        mmr = float(num_reads/1000000)
        fragments = self.num_f + self.num_r
        if fragments == 0:
            self.fpkm = 0
        else:
            self.fpkm = float((fragments/self.length)/mmr)
        return

    def calc_tpm(self, denominator):
        """
        Divide the read counts by the length of each gene in kilobases.
        Count up all the RPK values in a sample and divide this number by 1,000,000.
        Divide the RPK values by the “per million” scaling factor.

        :param denominator: The per-million scaling factor
        :return:
        """
        self.tpm = self.rpk/denominator
        return

    def get_info(self):
        return ['"' + self.name + '"',
                str(self.length),
                str(self.num_f), str(self.num_r),
                str(self.fpkm), str(self.tpm)]


def read_frhit_table(args):
    hit_table = args.alignment_table
    min_proportion = args.alignment_proportion
    min_identity = args.identity
    hits_list = list()

    try:
        input_handler = open(hit_table, 'r')
    except IOError:
        raise IOError("ERROR: Unable to open " + hit_table + " for reading!\n")

    line = input_handler.readline()
    # ReadName	ReadLength	E-value	AlignmentLength	Begin	End	Strand	Identity	ReferenceSequenceName	Begin	End
    if not re.match(r'^\S+\t\S+\t\S+\t[0-9]+\t[0-9]+\t[0-9]+\t[+-]\t\S+\t\S+\t[0-9]+\t[0-9]+$', line):
        sys.stderr.write("ERROR: line in " + hit_table + " is not of the expected format!\n")
        sys.exit(2)
    x = 0
    # while line and x < 1000:
    while line:
        hit = Mapping()
        fields = line.split('\t')
        # Parse the relevant information from the line and load it into the Mapping instance
        hit.load_dat(fields)
        # Does it meet the thresholds?
        if hit.proportion >= min_proportion and hit.identity >= min_identity:
            hits_list.append(hit)
        # Move on to the next line
        line = input_handler.readline()
        x += 1

    input_handler.close()

    if args.verbose:
        sys.stdout.write(str(x) + " lines in alignment table.\n")
        sys.stdout.write(str(len(hits_list)) + " alignments passing thresholds.\n")
        sys.stdout.flush()

    return hits_list


def read_fasta(args):
    try:
        fasta_handler = open(args.genome, 'r')
    except IOError:
        raise IOError("ERROR: Unable to open the FASTA file " + args.genome + " for reading!")

    line = fasta_handler.readline()
    # Check FASTA format
    if line[0] != '>':
        sys.stderr.write("ERROR: First line in FASTA file is not a header!\n")
        sys.exit(4)
    genome_dict = dict()
    sequence = ""
    header = ""
    while line:
        line = line.strip()
        if line and line[0] == '>':
            if sequence:
                genome_dict[header] = RefSeq(header, sequence)
            header = line
            sequence = ""
        else:
            sequence += line
        line = fasta_handler.readline()

    genome_dict[header] = RefSeq(header, sequence)

    return genome_dict


def add_recruitments(hits_list, genome_dict):
    header_replacements = dict()

    for recruitment in hits_list:
        header = '>' + recruitment.reference
        if header not in genome_dict.keys():
            if header in header_replacements:
                header = header_replacements[header]
            else:
                found = False
                # Attempt to rescue by splitting on whitespace
                for fasta_header in genome_dict.keys():
                    if re.search(header, fasta_header):
                        # fr-hit split on whitespace
                        if fasta_header.split(' ')[0] == header:
                            header_replacements[header] = fasta_header
                            header = fasta_header
                            found = True
                            break

                if not found:
                    sys.stderr.write("ERROR: Unable to find " + recruitment.reference + " in provided FASTA headers!\n")
                    sys.exit(5)

        ref_seq = genome_dict[header]
        if recruitment.strand == '+':
            ref_seq.num_f += 1
        elif recruitment.strand == '-':
            ref_seq.num_r += 1
    return


def calculate_normalization_metrics(genome_dict, sampled_reads):
    rpk_sum = 0
    for header in sorted(genome_dict.keys()):
        ref_seq = genome_dict[header]
        ref_seq.calc_fpkm(sampled_reads)
        ref_seq.rpk = sampled_reads/(ref_seq.length/1000)
        rpk_sum += ref_seq.rpk

    denominator = rpk_sum / 1E6
    for header in genome_dict.keys():
        genome_dict[header].calc_tpm(denominator)

    return


def group_by(ref_seqs: dict, grouping=None):
    if grouping == "genome":
        genome_ref = RefSeq("Reference", "")
        for seq_name in ref_seqs:
            ref_seq = ref_seqs[seq_name]  # type: RefSeq
            genome_ref.length += ref_seq.length
            genome_ref.num_f += ref_seq.num_f
            genome_ref.num_r += ref_seq.num_r
            genome_ref.fpkm += ref_seq.fpkm
            genome_ref.tpm += ref_seq.tpm
        ref_seqs = {genome_ref.name: genome_ref}
    else:
        pass

    return ref_seqs


def write_summary(args, genome_dict):
    try:
        output_handler = open(args.output, 'w')
    except IOError:
        raise IOError("Unable to open " + args.output + " for writing!\n")

    output_buffer = "\t".join(["Reference", "Query"]) + "\n"
    for contig_header in genome_dict.keys():
        ref_seq = genome_dict[contig_header]  # type: RefSeq
        # print('\t'.join(ref_seq.get_info()))
        output_buffer += args.reference + '\t' + args.query + '\t' + '\t'.join(ref_seq.get_info()) + "\n"
        if len(output_buffer) >= 1E6:
            output_handler.write(output_buffer)
            output_buffer = ""

    output_handler.write(output_buffer)
    output_handler.close()

    return


def main():
    args = get_options()
    hits_list = read_frhit_table(args)
    genome_dict = read_fasta(args)
    add_recruitments(hits_list, genome_dict)
    calculate_normalization_metrics(genome_dict, args.num_reads)
    write_summary(args, genome_dict)


main()
