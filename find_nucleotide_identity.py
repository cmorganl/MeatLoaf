#!/usr/bin/env python

import argparse
import subprocess
import itertools
import glob
import os
import sys


def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--genome_directory", required=True,
                        help="Path to directory containing FASTA files")
    parser.add_argument("-l", "--list", required=True,
                        help="A file listing the FASTA files to be included in all-vs-all whole genome alignment")
    parser.add_argument("-o", "--output_matrix", required=False, default="genome_similarity_matrix.tsv",
                        help="Name of file to write the identity matrix")
    args = parser.parse_args()
    return args


def get_genome_size(fasta):
    genome_cat = ""
    with open(fasta) as genome:
        for line in genome:
            if line[0] != '>':
                genome_cat += line.strip()
    return len(genome_cat)


def align_genomes(args):
    """
    Call nucmer for aligning genomes and reads the genomes to determine their size
    :param args:
    :return:
    """
    sys.stdout.write("Running nucmer... ")
    sys.stdout.flush()
    genome_size_dict = dict()
    fasta_files = list()
    nucmer_command = ["nucmer", "-l", str(100), "--coords", "-p"]

    # Read in FASTA files from list
    with open(args.list) as fasta_list:
        for line in fasta_list:
            fasta = line.strip()
            fasta_path = args.genome_directory + fasta
            if not os.path.isfile(fasta_path):
                raise IOError(fasta_path + "does not exist!")
            else:
                fasta_files.append(fasta_path)

    # Perform all-vs-all alignments
    for combo in itertools.permutations(fasta_files, 2):
        ref, query = combo
        ref_prefix = '.'.join(os.path.basename(ref).split('.')[0:-1])
        query_prefix = '.'.join(os.path.basename(query).split('.')[0:-1])

        # Get the number of bases in the FASTA file
        if ref_prefix not in genome_size_dict.keys():
            genome_size_dict[ref_prefix] = get_genome_size(ref)

        nucmer_prefix = args.genome_directory + ref_prefix + ".TO." + query_prefix
        launch_nucmer = nucmer_command + [nucmer_prefix, ref, query]
        launch_nucmer += ["1>", "/dev/null", "2>", "/dev/null"]
        if not os.path.isfile(nucmer_prefix + ".coords"):
            p_nucmer = subprocess.Popen(' '.join(launch_nucmer), shell=True, preexec_fn=os.setsid)
            p_nucmer.wait()
            if p_nucmer.returncode != 0:
                raise RuntimeError("nucmer did not complete successfully!")

    sys.stdout.write("done\n")
    sys.stdout.flush()
    return genome_size_dict


def get_coords_files(args):
    coords_files = glob.glob(args.genome_directory + "*coords")
    if len(coords_files) == 0:
        raise RuntimeError("Unable to find nucmer outputs")
    delta_files = glob.glob(args.genome_directory + "*delta")
    for delta in delta_files:
        os.remove(delta)
    return coords_files


def whitespace_delimiter_parser(text):
    """
    Function to return a list of strings parsed from text containing random whitespace as a delimiter
    :param text: string
    :return: list of strings
    """
    text_strings = list()
    curr_string = ""
    for c in text:
        if c == " " and curr_string:
            text_strings.append(curr_string)
            curr_string = ""
        elif c == " ":
            pass
        else:
            curr_string += c
    return text_strings


def parse_coords(coords_files):
    """
    Function to parse the reference length aligned and identities from all files in coords_files
    :param coords_files:
    :return:
    """
    sys.stdout.write("Parsing nucmer coords files... ")
    sys.stdout.flush()
    identity_dict = dict()

    for coords_f in coords_files:
        ref, query = coords_f.split(".TO.")
        ref = os.path.basename(ref)
        query = str(query).replace(".coords", "")
        if ref not in identity_dict.keys():
            identity_dict[ref] = dict()
        coords = open(coords_f, 'r')
        line = coords.readline()
        while line[0] != "=":
            line = coords.readline()
        line = coords.readline().strip()
        total_identical_bp = 0
        while line:
            fields = whitespace_delimiter_parser(line)
            ref_align_len = fields[6]
            identity = fields[9]
            total_identical_bp += (float(ref_align_len) * float(identity)) / 100
            line = coords.readline().strip()

        identity_dict[ref][query] = total_identical_bp
        coords.close()

    sys.stdout.write("done\n")
    sys.stdout.flush()
    return identity_dict


def write_identity_matrix(identity_dict, genome_size_dict, output_matrix):
    matrix = open(output_matrix, 'w')
    for reference in identity_dict:
        number_bases = genome_size_dict[reference]
        for query in identity_dict[reference]:
            identity = identity_dict[reference][query] / number_bases
            matrix.write(reference + "\t" + query + "\t" + str(identity) + "\n")
    matrix.close()
    return


def main():
    args = get_options()
    genome_size_dict = align_genomes(args)
    coords_files = get_coords_files(args)
    identity_dict = parse_coords(coords_files)
    write_identity_matrix(identity_dict, genome_size_dict, args.output_matrix)

main()
