#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'

import argparse
import sys
import re
from os import path


def get_options():
    parser = argparse.ArgumentParser(add_help=False)
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument("-a", "--alignment_table", dest="aln_tab",
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
    optopt.add_argument('-g', "--group", default="contig", choices=["contig", "genome", "variable"], required=False,
                        help="reports summary stats for individual contigs or entire reference. [ DEFAULT = contig ]")
    optopt.add_argument('-x', "--variables", required=False, dest="grouping_vars",
                        help="If the output contains several genomes, this can be either a file or comma-separated list"
                             "of prefixes that relates contigs from the same genome to each other.",
                        default=None,
                        type=str)
    optopt.add_argument('-o', '--output', default='./fr-hit_output_summary.tsv', required=False,
                        help='the output file [ DEFAULT = ./fr-hit_output_summary.tsv ]')
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
    optopt.add_argument("-l", "--library",
                        help="The sequencing library type, either paired-end (pe) or single-end (se).",
                        default="pe", choices=["pe", "se"])

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
        sys.exit(3)

    return args


def read_group_vars(grouping_vars: str) -> list:
    """
    Turns the args.grouping_vars string into a list if it either points to a file or is a comma-separated list of
    strings provided on the command-line.

    :param grouping_vars: String representing either a file path or a list of strings
    :return: List of strings that will be used to group contigs
    """
    vars_list = []
    if grouping_vars:
        if path.isfile(grouping_vars):
            with open(grouping_vars) as groups_handler:
                vars_list = [line.strip() for line in groups_handler.readlines()]
        else:
            vars_list = grouping_vars.split(',')
    return vars_list


class Mapping:
    def __init__(self):
        self.proportion = 0
        self.identity = 0
        self.weight = 0
        self.strand = ""
        self.reference = ""
        self.query_name = ""

    def load_frhit_dat(self, fields):
        if len(fields) != 11:
            sys.stderr.write("ERROR: fr-hit output line is not of the expected format!\n")
            sys.exit(3)
        read_len = int(re.sub('nt|aa', '', fields[1]))
        self.proportion = float(int(fields[3])/read_len)
        self.identity = float(re.sub('%', '', fields[7]))
        self.reference = fields[8]
        self.query_name = fields[0]
        self.strand = fields[6]
        return


class RefSeq:
    def __init__(self, header, sequence):
        self.name = header
        self.sequence = sequence
        self.length = len(sequence)
        self.num_pos = 0
        self.num_neg = 0
        self.weight_sum = 0
        self.fpkm = 0
        self.rpk = 0
        self.tpm = 0

    def aggregate(self, ref_seq) -> None:
        self.length += ref_seq.length
        self.num_pos += ref_seq.num_pos
        self.num_neg += ref_seq.num_neg
        self.weight_sum += ref_seq.weight_sum
        self.fpkm += ref_seq.fpkm
        self.rpk += ref_seq.rpk
        self.tpm += ref_seq.tpm

    def calc_fpkm(self, num_reads):
        mmr = float(num_reads/1E6)
        if self.weight_sum == 0:
            self.fpkm = 0
        else:
            self.fpkm = float((self.weight_sum/self.length)/mmr)
        return

    def calc_tpm(self, denominator) -> None:
        """
        Divide the read counts by the length of each gene in kilobases.
        Count up all the RPK values in a sample and divide this number by 1,000,000.
        Divide the RPK values by the “per million” scaling factor.

        :param denominator: The per-million scaling factor
        :return: None
        """
        if self.weight_sum == 0:
            return
        self.tpm = self.rpk/denominator
        return

    def get_info(self) -> list:
        """
        Used to return a list of variables summarizing the reference sequence.

        :return: A list of the reference sequence's:
         name, length, number of forward and reverse reads mapped to it, FPKM and TPM
        """
        return [self.name,
                str(self.length),
                str(self.num_pos), str(self.num_neg),
                str(self.fpkm), str(self.tpm)]


def read_frhit_table(hit_table: str, min_proportion: float, min_identity: int, verbose=False) -> list:
    """
    Parses the mapped query sequences (e.g. reads) from an fr-hit output table and filters the alignments.
    All alignments that pass the min_identity and min_proportion thresholds become Mapping instances

    :param hit_table: Path to an alignment table generated by fr-hit
    :param min_proportion: Minimum proportion of a read aligned to be included
    :param min_identity: Minimum alignment similarity to the reference sequence to be included
    :param verbose: Boolean indicating whether status and summary statements should be printed to screen
    :return: A list of Mapping instances
    """
    hits_list = list()

    try:
        input_handler = open(hit_table, 'r')
    except IOError:
        raise IOError("ERROR: Unable to open " + hit_table + " for reading!\n")

    if verbose:
        sys.stdout.write("Reading fr-hit alignment table... ")
        sys.stdout.flush()

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
        hit.load_frhit_dat(fields)
        # Does it meet the thresholds?
        if hit.proportion >= min_proportion and hit.identity >= min_identity:
            hits_list.append(hit)
        # Move on to the next line
        line = input_handler.readline()
        x += 1

    input_handler.close()

    if verbose:
        sys.stdout.write("done.\n")
        sys.stdout.write(str(x) + " lines in alignment table.\n")
        sys.stdout.write(str(len(hits_list)) + " alignments passing thresholds.\n")
        sys.stdout.flush()

    return hits_list


def read_fasta(fasta_file, verbose=False) -> dict:
    """
    Its all in the title.

    :param fasta_file: Path to a file in FASTA format
    :param verbose: Boolean controlling the verbosity of standard output statements
    :return: Dictionary of RefSeq instances, one for each sequence read from FASTA, indexed by sequence names
    """
    try:
        fasta_handler = open(fasta_file, 'r')
    except IOError:
        raise IOError("ERROR: Unable to open the FASTA file " + fasta_file + " for reading!")

    if verbose:
        sys.stdout.write("Reading FASTA file... ")
        sys.stdout.flush()

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

    if verbose:
        sys.stdout.write("done.\n")
        sys.stdout.write("Read %d sequences from FASTA file.\n" % len(genome_dict.keys()))

    return genome_dict


def add_recruitments(hits_list, genome_dict) -> None:
    """
    Iterates through the list of alignments/hits and increments the number of reads aligned to either the
     positive (num_pos) or negative (num_neg) strand of the reference sequence. The RefSeq's weight_sum variable is also
     increased according to a particular alignment's weight, which is determined by the number of different reference
     sequences it aligned to and whether the library sequenced was paired-end or single-end.
    
    :param hits_list: A list of Mapping instances
    :param genome_dict: A dictionary of RefSeq instances indexed by their respective headers
    :return: None
    """
    header_replacements = dict()

    for recruitment in hits_list:  # type: Mapping
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

        ref_seq = genome_dict[header]  # type: RefSeq
        if recruitment.strand == '+':
            ref_seq.num_pos += 1
        elif recruitment.strand == '-':
            ref_seq.num_neg += 1
        ref_seq.weight_sum += recruitment.weight

    return


def calculate_normalization_metrics(genome_dict: dict, sampled_reads: int) -> None:
    """
    Calculates the normalized abundance values for each header's RefSeq instance in genome_dict
        1. Reads per kilobase (RPK) is calculated using the reference sequence's length and number of reads (provided
        by the user via CLI)
        2. Fragments per kilobase per million mappable reads (FPKM) is calculated from the number of fragments
        (this is different from reads by, in a paired-end library, forward and reverse pair makes up one fragment)
        normalized by the reference sequence length and the number of reads mapped.
        2. Transcripts per million (TPM) is calculated similarly to FPKM but the order of operations is different.

    :param genome_dict: A dictionary of RefSeq instances indexed by headers (sequence names)
    :param sampled_reads: The number of reads sequenced (not aligned). Required for normalizing by sequencing depth
    :return: None
    """
    rpk_sum = 0  # The total reads per kilobase (RPK) of all reference sequences
    for header in sorted(genome_dict.keys()):
        ref_seq = genome_dict[header]  # type: RefSeq
        if ref_seq.weight_sum == 0:
            continue
        ref_seq.calc_fpkm(sampled_reads)
        ref_seq.rpk = sampled_reads/(ref_seq.length/1E3)
        rpk_sum += ref_seq.rpk

    denominator = rpk_sum / 1E6
    for header in genome_dict.keys():
        ref_seq = genome_dict[header]  # type: RefSeq
        ref_seq.calc_tpm(denominator)

    return


def group_by(ref_seqs: dict, grouping=None, grouping_vars=None) -> dict:
    """
    Groups together RefSeq instances in ref_seqs either based on a variable provided or all reference sequences if
    group is set to 'genome'. When grouping, the numerical variables in RefSeq are summed so each group has the
    aggregate value across all RefSeq instances. Otherwise nothing is changed.

    :param ref_seqs: Dictionary of RefSeq instances that have been populated with at least the weights
    :param grouping: Specification of the groups, whether its 'genome', 'contig', or 'variable'
    :param grouping_vars: If grouping is 'variable' these strings are used in regular expressions to find the matching
    RefSeq instances by their sequence names (headers)
    :return: Dictionary of RefSeq instances for each of the new groups.
    """
    grouped_ref_seqs = {}
    if grouping_vars is None:
        grouping_vars = []

    if grouping == "genome":
        genome_ref = RefSeq("Reference", "")
        for seq_name in ref_seqs:
            ref_seq = ref_seqs[seq_name]  # type: RefSeq
            genome_ref.aggregate(ref_seq)
        grouped_ref_seqs[genome_ref.name] = genome_ref
        return grouped_ref_seqs
    elif grouping == "variable" and grouping_vars:
        ref_names = list(ref_seqs.keys())
        for var in grouping_vars:
            genome_ref = RefSeq(var, "")
            found = False
            i = 0
            while i < len(ref_names):
                ref = ref_names[i]
                if re.search(var, ref):
                    ref_seq = ref_seqs[ref]  # type: RefSeq
                    genome_ref.aggregate(ref_seq)
                    ref_names.pop(i)
                    found = True
                else:
                    i += 1
            if not found:
                sys.stderr.write("ERROR: Unable to find grouping variable '%s' in reference sequence names.\n" % var)
                sys.exit(3)
            grouped_ref_seqs[genome_ref.name] = genome_ref
        return grouped_ref_seqs
    else:
        pass

    return ref_seqs


def distribute_read_weights(hits_list: list, chemistry="pe", verbose=False) -> None:
    """
    Calculates a read's weight variable based on whether it aligned to multiple reference sequences

    :param hits_list: A list of Mapping instances
    :param chemistry: Either "pe" or "se" indicating whether the reads have mate-pairs (paired-end chemistry) or not
    :param verbose: Boolean indicating whether status and summary statements should be printed to screen
    :return: None
    """
    query_hits_dict = dict()
    if chemistry == "pe":
        numerator = 0.5
    elif chemistry == "se":
        numerator = 1
    else:
        sys.stderr.write("ERROR: Unrecognized chemistry passed to distribute_read_weights: '%s'\n" % chemistry)
        sys.exit(3)

    for hit in hits_list:  # type: Mapping
        if hit.query_name not in query_hits_dict:
            query_hits_dict[hit.query_name] = []
        query_hits_dict[hit.query_name].append(hit)
    for query_name in query_hits_dict:
        num_alignments = len(query_hits_dict[query_name])
        for hit in query_hits_dict[query_name]:  # type: Mapping
            hit.weight = numerator/num_alignments

    if verbose:
        sys.stdout.write("Found %d unique query names in the alignments\n" % len(query_hits_dict))
    return


def write_summary(args, genome_dict):
    try:
        output_handler = open(args.output, 'w')
    except IOError:
        raise IOError("Unable to open " + args.output + " for writing!\n")

    output_buffer = "\t".join(["Reference", "Query", "SeqName", "SeqLen", "Positive", "Negative", "FPKM", "TPM"]) + "\n"
    for contig_header in sorted(genome_dict, key=lambda x: genome_dict[x].weight_sum, reverse=True):
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
    vars_list = read_group_vars(args.grouping_vars)
    hits_list = read_frhit_table(args.aln_tab, args.alignment_proportion, args.identity, args.verbose)
    genome_dict = read_fasta(args.genome, args.verbose)
    distribute_read_weights(hits_list, args.library, args.verbose)
    add_recruitments(hits_list, genome_dict)
    genome_dict = group_by(genome_dict, args.group, vars_list)
    calculate_normalization_metrics(genome_dict, args.num_reads)
    write_summary(args, genome_dict)


main()
