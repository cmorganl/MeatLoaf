#!/usr/bin/env python
__author__ = "Connor Morgan-Lang"

import sys
import argparse

def set_arguments():
    parser = argparse.ArgumentParser()
    extract = parser.add_argument_group("extract")
    parser.add_argument("-i", "--fasta", help = "The fasta file to be subsetted.", required=True)
    parser.add_argument("-l", "--list", help="The list of contig headers to be subsetted from the fasta file", required=True)
    parser.add_argument("-o", "--output", help="The output fasta file [DEFAULT = subset.fasta]", required=False, default="subset.fasta")
    parser.add_argument("-N", "--Nsplit", help="Split the contigs on ambiguity character 'N's [DEFAULT = False]", \
                        default=False, action="store_true")
    parser.add_argument("-m", "--minLength", help="The minimum contig length [DEFAULT = 200].", \
                        required=False, default=200, type=int)
    parser.add_argument("-x", "--longest", help="Take the 'x' longest contigs. Writes all contigs by default.", default=0, type=int)
    extract.add_argument("-s", "--slice", help="Flag indicating a specific sequence should be extracted from a specific header." ,\
                         required=False, action="store_true")
    extract.add_argument("-p", "--pos", help="Positions (start,end) of the provided contig to extract", required=False)
    args = parser.parse_args()
    return args

def subset_fasta(fasta, headers):
    """
    :return: A dictionary with the contig headers of interest as keys and their respective sequences as their values.
    """
    subset = dict()
    header = ""
    add = 0
    with open(fasta) as fas:
        line = fas.readline()
        while line:
            if ( line[0] == '>' ):
                header = line[1:].strip()
                if headers[0] == "all":
                    headers.append(header)
                if header in headers:
                    subset[header] = ""
                    add = 1
                else:
                    add = 0
            elif (add == 1):
                subset[header] += line.strip()
            else:
                pass
            line = fas.readline()
    return subset

def load_list(header_list):
    headers = list()
    with open(header_list) as LoH:
        for line in LoH:
            if ( line[0] == '>'):
                headers.append(line[1:].strip())
            else:
                headers.append(line.strip())
    return headers

def extract_seq(contig_seq, pos):
    """
    :return: A sub-string of the contig sequence as determined by the pos tuple.
    """
    start, end = pos
    return str(contig_seq[start:end])

def Nsplit_scaffolds(subset, minLength):
    """
    Splits all scaffolds containing 'N's on 'N' and
    adds these contigs to the dictionary to be returned.
    Removes all sub-contigs that are smaller than minLength.
    """
    unambiguous_subset = dict()
    acc = 1
    for head in subset.keys():
        if 'N' in subset[head]:
            contigs = subset[head].split('N')
            for c in contigs:
                if len(c) > minLength:
                    contig_name = head + "_" + str(acc)
                    acc += 1
                    unambiguous_subset[contig_name] = c
        else:
            contig_name = head + "_" + str(acc)
            acc += 1
            unambiguous_subset[contig_name] = subset[head]
        acc = 1
    return unambiguous_subset

def slice_contigs(subset, pos):
    contigs = dict()
    for contig in subset:
        seq = extract_seq(subset[contig], pos)
        contigs[contig] = seq
    return sorted(contigs)

def write_subset_to_output(subset, output, minLength, longest):
    n = 1
    with open(output, 'w') as fa_out:
        for contig in sorted(subset, key=lambda contig: len(subset[contig]), reverse=True):
            if (len(subset[contig]) > minLength and n <= longest):
                fa_out.write('>' + str(contig) + "\n")
                fa_out.write(subset[contig] + "\n")
                n += 1
    return 0

def main():
    print "Beginning", sys.argv[0]
    args = set_arguments()
    headers = ["all"]
    if args.list != "all":
        headers = load_list(args.list)
    print "Finding the sequences in", args.list
    subset = subset_fasta(args.fasta, headers)
    if args.Nsplit:
        print "Splitting the sequences on 'N's..."
        subset = Nsplit_scaffolds(subset, args.minLength)
    if args.slice:
        subset = slice_contigs(subset, args.slice, args.pos)
    print "Writing the contig subset to", args.output
    write_subset_to_output(subset, args.output, args.minLength, args.longest)

    return 0

main()
