#!/usr/bin/env python3

__author__ = 'Connor Morgan-Lang'


import argparse
import json
import re


def get_options():
    parser = argparse.ArgumentParser(add_help=False)
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument("-m", "--marker_genes",
                               help="marker_gene_stats.tsv file written by checkM containing all genome stats",
                               required=True)

    optopt = parser.add_argument_group("Optional options")
    optopt.add_argument("-c", "--complete",
                        help="The target completeness of selected genomes [ DEFAULT = 0.90 (90%) ] ",
                        required=False, default=0.90, type=float)

    miscellaneous_opts = parser.add_argument_group("Miscellaneous options")
    miscellaneous_opts.add_argument('-v', '--verbose', action='store_true', default=False,
                                    help='prints a more verbose runtime log')
    miscellaneous_opts.add_argument("-h", "--help",
                                    action="help",
                                    help="show this help message and exit")

    args = parser.parse_args()

    return args


def read_marker_gene_stats(marker_genes: str) -> dict:
    genome_marker_dict = dict()
    # Parse the file, loading the markers associated with each genome into the dictionary
    with open(marker_genes, 'r') as mgs_handler:
        for line in mgs_handler:
            genome, json_stats = line.split("\t")
            genome_marker_dict[genome] = set()
            json_stats = re.sub("'", "\"", json_stats)
            marker_stats = json.loads(json_stats)
            for contig in marker_stats:
                for mg in marker_stats[contig]:
                    genome_marker_dict[genome].add(mg)
    return genome_marker_dict


def find_optimal_set(genome_marker_dict, complete_target=0.90) -> list:
    optimal_genomes = list()
    optimal_mg = set()
    unique_markers = set()
    max_mg = 0
    centroid = ""
    # Find the number of unique markers identified
    for genome in genome_marker_dict:
        unique_markers.update(genome_marker_dict[genome])
        # Identify the genome that is most complete
        if len(genome_marker_dict[genome]) > max_mg:
            centroid = genome
            max_mg = len(genome_marker_dict[genome])
    print("%d unique marker genes identified across all genomes." % len(unique_markers))
    target = round(complete_target * len(unique_markers))
    print("Identifying optimal set of genomes to reach %d marker genes." % target)

    print("Genome '%s' was selected as the centroid with %d marker genes." % (centroid, max_mg))
    # Iteratively add the next genome that best complements the genome centroid
    optimal_genomes.append(centroid)
    optimal_mg.update(genome_marker_dict[centroid])

    while len(optimal_mg) < target and len(optimal_genomes) < len(genome_marker_dict):
        max_new_mg = 0
        to_absorb = ""
        for genome in genome_marker_dict:
            if genome in optimal_genomes:
                continue
            query_set = genome_marker_dict[genome]  # type: set
            new_mg = len(query_set.union(optimal_mg)) - len(optimal_mg)
            if new_mg > max_new_mg:
                max_new_mg = new_mg
                to_absorb = genome
        optimal_mg.update(genome_marker_dict[to_absorb])
        optimal_genomes.append(to_absorb)

    print("Optimal set includes %d unique marker genes." % len(optimal_mg))

    return optimal_genomes


def print_summary(optimal_set) -> None:
    print("The following genomes represent the fewest genomes required to reach the target completeness:")
    print("\n".join(optimal_set))
    return


def main():
    args = get_options()
    genome_marker_dict = read_marker_gene_stats(args.marker_genes)
    optimal_set = find_optimal_set(genome_marker_dict, args.complete)
    print_summary(optimal_set)


main()
