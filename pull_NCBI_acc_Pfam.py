#!/usr/bin/env python3

import argparse
import sys
import re

__author__ = 'Connor Morgan-Lang'


def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        required=True,
                        help="File containing NCBI accession IDs to be downloaded")
    parser.add_argument("-p", "--PfamId",
                        required=True, dest="pfam_id",
                        help="Comma-separated list of the PFam IDs to pull hits from")
    parser.add_argument("-o", "--output_prefix",
                        required=False, default="./",
                        help="Prefix name for the output file(s) - one for each PfamId provided")
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False,
                        help="Print what is happening at every stage.")

    args = parser.parse_args()

    return args


def parse_pfam(pfam_ncbi_file: str, pfam_ids: list):
    pfam_hit_dict = dict()

    # Strip away the version numbers for the PFam identifiers as these are not matched to
    no_ver_pfams = [pfam_id.split('.')[0] for pfam_id in pfam_ids]
    for pfam_id in no_ver_pfams:
        pfam_hit_dict[pfam_id] = []

    # Regular expressions for differentiating between the lines with hits, and those with other data
    p_acc_re = re.compile("^#=GF AC\s+(PF[0-9]+)\.[0-9]{1,2}$")
    hit_line_re = re.compile("^#=GS\s(.*)")

    # Open up the PFam hit file
    try:
        pfam_file_handler = open(pfam_ncbi_file, encoding="ISO-8859-1")
    except OSError:
        sys.stderr.write("ERROR: Unable to open " + pfam_ncbi_file + " for reading!\n")
        sys.exit(1)

    sys.stdout.write("Searching for a PFam ID match... ")
    sys.stdout.flush()
    for line in pfam_file_handler:
        if p_acc_re.match(line):
            if p_acc_re.match(line).group(1) in no_ver_pfams:
                # Update the current PFam ID being parsed
                curr_id = p_acc_re.match(line).group(1)
                # Remove the current PFam ID from the list of PFam IDs
                no_ver_pfams.pop(no_ver_pfams.index(curr_id))

                sys.stdout.write("found " + curr_id + ".\n")
                sys.stdout.flush()
                sys.stdout.write("Reading accessions... ")
                sys.stdout.flush()
                # Read all of the PFam HMM hits and parse out the NCBI accessions
                for sub_line in pfam_file_handler:
                    hit_line_match = hit_line_re.match(sub_line)
                    if hit_line_match:
                        annots = hit_line_match.group(1)
                        # Strip away the region that matched the PFam domain
                        ncbi_acc = re.sub("/[0-9]+-[0-9]+", '', annots.split(" ")[0])
                        pfam_hit_dict[curr_id].append(ncbi_acc)
                    elif p_acc_re.match(sub_line):
                        sys.stdout.write("done.\n")
                        if len(no_ver_pfams) == 0:
                            pfam_file_handler.close()
                            return pfam_hit_dict
                        sys.stdout.write("Searching for a PFam ID match... ")
                        sys.stdout.flush()
                        break
                    else:
                        pass
    pfam_file_handler.close()

    return pfam_hit_dict


def write_pfam_hit_lists(output_prefix: str, pfam_hit_dict: dict):
    sys.stdout.write("PFam-ID\tHits\n")
    for pfam_id in pfam_hit_dict:
        p_list_file = output_prefix + pfam_id + "_ncbi.txt"
        sys.stdout.write(pfam_id + "\t" + str(len(pfam_hit_dict[pfam_id])) + "\n")
        try:
            p_list_handler = open(p_list_file, 'w')
        except IOError:
            sys.stderr.write("ERROR: Unable to open " + p_list_file + " for writing!\n")
            sys.exit(1)

        p_list_handler.write("\n".join(pfam_hit_dict[pfam_id]) + "\n")
    return


def main():
    args = get_options()
    pfam_id_list = args.pfam_id.split(',')
    pfam_hit_dict = parse_pfam(args.input, pfam_id_list)
    write_pfam_hit_lists(args.output_prefix, pfam_hit_dict)
    return


main()
