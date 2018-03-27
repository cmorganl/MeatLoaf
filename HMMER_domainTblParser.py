#!/usr/bin/env python3

import os
import sys
import re
import argparse
import math


def get_options():
    parser = argparse.ArgumentParser(description="Parses a domtbl file generated by HMMER, applying its own set of"
                                                 "filters and summarizing those high-quality matches."
                                                 "Optionally, sequences of the high-quality matches can also be written"
                                                 "to a FASTA file.", add_help=False)
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument("-i",
                               dest="input",
                               help="Input domain table from HMMER to be parsed",
                               required=True)

    opt_args = parser.add_argument_group("Optional arguments")
    opt_args.add_argument("-f",
                          dest="fasta_in", required=False,
                          help="Path to a FASTA file containing the sequence that were searched against the HMM(s). "
                               "If one isn't provided, the script just prints summary stats.")
    opt_args.add_argument("-o",
                          dest="output", required=False, default="hmm_purified.fasta",
                          help="The name of the FASTA file to write containing sequences of the high-quality hits. "
                               "[default=hmm_purified.fasta]")
    opt_args.add_argument("-p",
                          dest="perc_aligned",
                          type=int,
                          default=90,
                          help="The minimum percentage of the HMM that was covered by the target sequence (ORF) "
                               "for the COG hit to be included [default=90]")
    opt_args.add_argument("-e",
                          dest="min_e",
                          type=float,
                          default=0.0001,
                          help="The largest E-value for the search to be accepted as significant [default=1E-3]")
    opt_args.add_argument("-a",
                          dest="min_acc",
                          type=float,
                          default=0.90,
                          help="The minimum acc threshold of the HMM search for reporting [default=0.90]")
    miscellaneous_opts = parser.add_argument_group("Miscellaneous options")
    miscellaneous_opts.add_argument("-v", "--verbose", action='store_true', default=False,
                                    help='Prints a more verbose runtime log')
    miscellaneous_opts.add_argument("-h", "--help",
                                    action="help",
                                    help="Show this help message and exit")
    args = parser.parse_args() 
    return args


def format_hmmer_domtbl_line(line):
    stats = []
    stat = ""
    for c in line:
        if c == ' ':
            if len(stat) > 0:
                stats.append(str(stat))
                stat = ""
            else:
                pass
        else:
            stat += c
    stats.append(str(stat))
    return stats


class DomainTableParser(object):

    def __init__(self, dom_tbl):
        self.alignments = {}
        self.i = 0
        self.lines = []
        self.size = 0
        try:
            self.commentPattern = re.compile(r'^#')
            self.src = open(dom_tbl)
        except IOError:
            sys.stderr.write("Could not open " + dom_tbl + " or file is not available for reading.\n")
            sys.exit(0)

    def __iter__(self):
        return self

    def read_domtbl_lines(self):
        """
        Function to read the lines in the domain table file,
        skipping those matching the comment pattern
        :return: self.lines is a list populated with the lines
        """
        line = self.src.readline()
        while line:
            comment = self.commentPattern.match(line)
            if not comment:
                self.lines.append(line.strip())
            if not line:
                break
            line = self.src.readline()
        self.size = len(self.lines)

    def next(self):
        """
        Reformat the raw lines of the domain table into
        an easily accessible hmm_domainTable format and perform
        QC to validate the significance of the alignments
        """
        if self.i < self.size:
            hit = format_hmmer_domtbl_line(self.lines[self.i])
            self.prepare_data(hit)
            self.i += 1
            
            try:
                return self.alignments
            except ValueError:
                return None
        else:
            self.src.close()
            return None
        
    def prepare_data(self, hit):
        self.alignments['query'] = str(hit[0])
        self.alignments['query_len'] = int(hit[2])
        self.alignments['hmm_name'] = str(hit[3])
        self.alignments['hmm_len'] = str(hit[5])
        self.alignments['Eval'] = float(hit[6])  # Full-sequence E-value (in the case a sequence alignment is split)
        self.alignments['num'] = int(hit[9])  # HMMER is able to detect whether there are multi-hits
        self.alignments['of'] = int(hit[10])  # This is the number of multi-hits for a query
        self.alignments['qstart'] = int(hit[17])
        self.alignments['qend'] = int(hit[18])
        self.alignments['acc'] = float(hit[21])
        self.alignments['desc'] = ' '.join(hit[22:])


def filter_poor_hits(args, dom_table):
    """
    Filters the homology matches based on their E-values and mean posterior probability of aligned residues from
    the maximum expected accuracy (MEA) calculation.
    Takes into account multiple homology matches of an ORF to a single gene and determines the total length of the
    alignment instead of treating them as individual alignments. This information is used in the next filtering step.
    """
    min_acc = float(args.min_acc)
    min_e = float(args.min_e)

    num_dropped = 0  # Obvious
    lines_parsed = 0  # Obvious
    multimatches = 0  # matches of the same query to the same HMM with multiple possible alignment positions (>1 lines)
    num_fragmented = 0  # matches that are broken into multiple alignments (indicated by num and of)
    cumulative_len = 0

    purified_matches = dict()
    hmm_size_dict = dict()
    previous_target = ""
    previous_query_len = ("", 0)

    if args.verbose:
        sys.stdout.write("\tMinimum E-value = " + str(min_e) + "\n")
        sys.stdout.write("\tMinimum acc = " + str(min_acc) + "\n")
        sys.stdout.write("\tPercentage of the HMM covered = " + str(args.perc_aligned) + "\n")
        sys.stdout.write("\tParsing matches from HMMER domain table... ")
        sys.stdout.flush()

    while dom_table.next():
        lines_parsed += 1
        data = dom_table.alignments
        if data['hmm_name'] not in hmm_size_dict.keys():
            hmm_size_dict[data['hmm_name']] = int(data['hmm_len'])
        if data['acc'] >= min_acc and data['Eval'] <= min_e:
                # Construct tuples with valuable information for the second filtering stage
                query_header = (data['query'], data['desc'])
                ali_len = data['qend'] - data['qstart']
                target_length = (data['hmm_name'], ali_len)
                if query_header != previous_query_len:
                    # This is a match between a new combination of target and query OR
                    # Same HMM (target) and a new ORF (query)
                    purified_matches[query_header] = list()
                    cumulative_len = 0
                    if data['of'] > 1:
                        num_fragmented += 1
                        # The target is added to purified_matches[query_header] in case the following alignment
                        # between the same target and query does not pass E-value or acc thresholds
                        match_start = int(data['qstart'])
                        match_end = int(data['qend'])
                        cumulative_len = ali_len
                    else:
                        purified_matches[query_header].append(target_length)

                elif data['hmm_name'] == previous_target and query_header == previous_query_len:
                    # Same HMM (target) and the same ORF (query)
                    if match_end > data['qstart']:
                        cumulative_len = data['qend'] - match_start
                    elif match_end < data['qstart']:
                        cumulative_len += ali_len

                    if data['of'] == 1:
                        multimatches += 1
                        purified_matches[query_header].append(target_length)
                    elif data['num'] < data['of']:
                        # More alignments should follow
                        match_end = data['qend']
                    else:
                        purified_matches[query_header].append((data['hmm_name'], cumulative_len))
                        
                elif data['hmm_name'] != previous_target and query_header == previous_query_len:
                    # New HMM (target), same ORF (query)
                    cumulative_len = 0
                    if data['of'] > 1:
                        # More alignments should follow
                        match_start = int(data['qstart'])
                        match_end = int(data['qend'])
                    else:
                        purified_matches[query_header].append((data['hmm_name'], ali_len))

                previous_query_len = query_header
                previous_target = data['hmm_name']
        else:
            num_dropped += 1

    if args.verbose:
        sys.stdout.write("done.\n")
    sys.stdout.write("Number of lines parsed:\t\t" + str(lines_parsed) + "\n")
    sys.stdout.write("Number of multi-matches:\t" + str(multimatches) + "\n")
    sys.stdout.write("Number of fragmented matches:\t" + str(num_fragmented) + "\n")

    return purified_matches, num_dropped, hmm_size_dict


def filter_incomplete_hits(args, purified_matches, num_dropped, hmm_size_dict):
    if args.verbose:
        sys.stdout.write("\tFiltering out short matches... ")
        sys.stdout.flush()
    complete_gene_hits = dict()

    for query in purified_matches:
        for target in purified_matches[query]:
            hmm_name, ali_len = target
            perc_aligned = (float(hmm_size_dict[hmm_name])/int(ali_len)*100)
            if perc_aligned >= args.perc_aligned:
                if hmm_name not in complete_gene_hits.keys():
                    complete_gene_hits[hmm_name] = []
                complete_gene_hits[hmm_name].append(query)
            else:
                num_dropped += 1

    if args.verbose:
        sys.stdout.write("done.\n")
    sys.stdout.write("Number of matches discarded:\t" + str(num_dropped) + "\n")
    sys.stdout.write("Matches for each HMM:\n")
    for hmm in complete_gene_hits:
        sys.stdout.write("\t" + hmm + "\t" + str(len(complete_gene_hits[hmm])) + "\n")
    sys.stdout.flush()
    return complete_gene_hits


def read_fasta(args):
    """
    :return: A dictionary with the contig headers of interest as keys and their respective sequences as their values.
    """
    fasta_dict = dict()
    try:
        fas = open(args.fasta_in, 'r')
    except IOError:
        raise IOError("ERROR: Unable to open fasta " + args.fasta_in + " for reading!\n")

    if args.verbose:
        sys.stdout.write("Reading " + args.fasta_in + "... ")
        sys.stdout.flush()

    header = ""
    line = fas.readline()
    while line:
        if line[0] == '>':
            header = line[1:].strip()
            fasta_dict[header] = ""
        else:
            if header == "":
                sys.stderr.write("ERROR: Header line is blank. FASTA file is not formatted correctly!\n")
                raise AssertionError
            fasta_dict[header] += line.strip()
        line = fas.readline()

    if args.verbose:
        sys.stdout.write("done.\n")

    if args.verbose:
        sys.stdout.write("\tNumber of sequences parsed from FASTA file: " + str(len(fasta_dict)) + "\n")

    return fasta_dict


def write_quality_matches_fasta(args, complete_gene_hits, fasta_dict):
    for marker in complete_gene_hits:
        if not args.output:
            output_file = marker + "_hmm_purified.fasta"
        else:
            output_file = args.output

        if args.verbose:
            sys.stdout.write("Writing " + marker + " sequences to " + output_file + "\n")

        try:
            if len(complete_gene_hits) > 1:
                output_handler = open(output_file, 'a')
            else:
                output_handler = open(output_file, 'w')
        except IOError:
            sys.stderr.write("ERROR: Unable to open " + output_file + " for appending!\n")
            raise IOError

        for query in complete_gene_hits[marker]:
            header = ' '.join(query)
            output_handler.write('>' + header + "\n")
            output_handler.write(fasta_dict[header] + "\n")
        output_handler.close()

    return


def main():
    """
    Drive the instantiation of the class and 
    further filtering of each of the hits
    """
    args = get_options()
    dom_table = DomainTableParser(args.input)
    if args.verbose:
        sys.stdout.write("\tReading HMMER domain table... ")
        sys.stdout.flush()
    dom_table.read_domtbl_lines()
    if args.verbose:
        sys.stdout.write("done.\n")
    purified_matches, num_dropped, hmm_size_dict = filter_poor_hits(args, dom_table)
    complete_gene_hits = filter_incomplete_hits(args, purified_matches, num_dropped, hmm_size_dict)
    if args.fasta_in:
        fasta_dict = read_fasta(args)
        write_quality_matches_fasta(args, complete_gene_hits, fasta_dict)


main()
