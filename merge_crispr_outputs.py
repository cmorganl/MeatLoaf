#!/usr/bin/python

import argparse
import os
import sys
import glob


def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("--taxa_map", required=True)
    parser.add_argument("--sag_16s_accessions", required=True)
    parser.add_argument("--crispr_outputs_list", required=True)
    parser.add_argument("-p", "--assembly_path", required=True,
                        help="The directory containing the assemblies in FASTA format")
    parser.add_argument("-o", "--output", required=True,
                        help="The table containing mapped information separated by commas")
    args = parser.parse_args()
    return args


def read_taxa(taxa_map):
    taxa_dict = dict()
    with open(taxa_map) as taxa_file_map:
        for line in taxa_file_map:
            count, id, taxonomy = line.split(',')
            taxa_dict[id] = taxonomy
    return taxa_dict


def read_16s_accessions(sag_accession_file):
    sag_16s_id = dict()
    with open(sag_accession_file) as headers:
        for line in headers:
            info = line.split(' ')[0][1:]
            full_sag_id, full_accession = info.split('|')
            sag = full_sag_id.split('.')[1][:-1]
            accession = full_accession.split('.')[0]
            sag_16s_id[sag] = accession
    return sag_16s_id


def read_crisprfinder_output(crisprfinder_output, headers):
    crispr_info = dict()
    contig = ""
    with open(crisprfinder_output) as txt:
        for line in txt:
            if line.startswith("# Id:"):
                contig_name = line.lstrip("# Id:")
                contig = contig_name.strip()
                if contig not in headers:
                    pass
                elif contig not in crispr_info:
                    crispr_info[contig] = list()
            if line.startswith("# DR:"):
                info = line.lstrip("# DR: ")
                dr_seq = info.split("\t")[0]
                num_spacers = 0
                dr_length = len(dr_seq)
            if not line[0] == '#':
                if contig in headers:
                    fields = line.split("\t")
                    if len(fields) > 1 and fields[0] != "Spacer_begin_position":
                        num_spacers += 1
                        spacer_begin = fields[0].lstrip()
                        spacer_seq = fields[2].lstrip()
                        crispr_info[contig].append([dr_seq, dr_length, num_spacers, spacer_begin, spacer_seq])
    return crispr_info


def merge_outputs(taxa_dict, sag_16s_id, crispr_info, sample, output_file):
    with open(output_file, 'a') as crispr_summary:
        accession = sag_16s_id[sample]
        taxonomy = taxa_dict[accession].strip()
        for contig in crispr_info:
            for spacer_info in crispr_info[contig]:
                dr_seq, dr_length, num_spacers, spacer_begin, spacer_seq = spacer_info
                line = [sample, accession, taxonomy, contig, dr_seq, str(dr_length),
                        str(num_spacers), str(spacer_begin), spacer_seq]
                crispr_summary.write("\t".join(line))
    return


def read_fasta(fasta):
    headers = list()
    line = fasta.readline()
    while line:
        if line[0] == ">":
            headers.append(line[1:].strip())
        else:
            pass
        line = fasta.readline()
    return headers


def main():
    args = get_options()
    taxa_dict = read_taxa(args.taxa_map)
    sag_16s_id = read_16s_accessions(args.sag_16s_accessions)
    with open(args.crispr_outputs_list) as crispr_outputs_list:
        for crisprfinder_output_path in crispr_outputs_list:
            crisprfinder_output = os.path.basename(crisprfinder_output_path)
            sample = '.'.join(crisprfinder_output.split('.')[:-1])
            sample_files = glob.glob(args.assembly_path + os.sep + sample + "*")
            if len(sample_files) > 1:
                sys.exit("ERROR: More than one valid file found for " + sample + " in " + args.assembly_path)
            fasta_file = sample_files[0]
            try:
                fasta = open(fasta_file, 'r')
            except:
                raise IOError("Unable to open " + fasta_file)
            headers = read_fasta(fasta)
            crispr_info = read_crisprfinder_output(crisprfinder_output_path.strip(), headers)
            merge_outputs(taxa_dict, sag_16s_id, crispr_info, sample, args.output)
main()
