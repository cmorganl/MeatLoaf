#!/usr/bin/env python

import os
import argparse
import sys


def get_options():
    parser = argparse.ArgumentParser(description="Script to pull pathway-names, pathway-common-names, rxn-common-name, "
                                                 "num-reactions, num-covered-reactions, orf-count, orf, taxonomy, and "
                                                 "RPKM values for all sequences from various MetaPathways outputs.")
    parser.add_argument("-i", "--hierarchy", required=True,
                        help="The hierarchy text file for mapping components to classes")
    parser.add_argument("-o", "--output", required=True,
                        help="The output file containing mapped values from the hierarchy")
    parser.add_argument("-l", "--comp_list", required=False, default=None,
                        help="An optional list of hierarchy components to return information on")
    parser.add_argument("--pathways", required=False, default=False, action="store_true",
                        help="Flag to indicate the output should be at pathway level [DEFAULT = reaction-level]")
    parser.add_argument("-p", "--pathway_annotations", required=True,
                        help="The output from extract_pathway_table_from_pgdb.pl")
    args = parser.parse_args()

    if not os.path.isfile(args.hierarchy):
        raise IOError("ERROR: " + args.hierarchy + " doesn't exist!")

    return args


def check_inputs(args):
    """
    Checks to make sure the paths and files exist and if the individual files are provided,
    then the metapathways directory is not provided
    :param args: command-line arguments parsed by argparse
    :return: nothing
    """
    # TODO: check input files
    return


def round_up_mp_inputs(args):
    args.sample_id = args.mp_inputs.rstrip(os.sep).split(os.sep)[-1]
    pwy_txt = args.mp_inputs + os.sep + "results" + os.sep + "pgdb" + os.sep + args.sample_id + ".pwy.txt"
    if not os.path.isfile(pwy_txt):
        raise IOError("Unable to find pathway annotations file (pwy.txt) in MetaPathway directory!\n"
                      "Looked for " + pwy_txt + " but file doesn't exist.")
    ptools_input = args.mp_inputs + os.sep + "ptools" + os.sep + "0.pf"
    if not os.path.isfile(ptools_input):
        raise IOError("Unable to find PGDB input file (0.pf) in MetaPathway directory!\n"
                      "Looked for " + ptools_input + " but file doesn't exist.")
    args.pathway_annotations = pwy_txt
    args.pgdb_input = ptools_input
    return args


def get_depth(fields):
    d = 0
    while fields[d] == "":
        d += 1
    return d


def read_hierarchy(hierarchy):
    prev_depth = -1
    path = list()
    paths_dict = dict()
    for line in hierarchy:
        fields = line.rstrip().split('\t')
        curr_depth = get_depth(fields)
        if curr_depth > prev_depth:
            path.append(fields[curr_depth])
        if prev_depth >= curr_depth:
            if path[0] not in paths_dict.keys():
                paths_dict[path[0]] = list()
            paths_dict[path[0]].append(path[1:])
            path = path[0:curr_depth]
            path.append(fields[curr_depth])
            # path.append("--".join(fields[curr_depth:]))
        prev_depth = curr_depth
    paths_dict[path[0]].append(path[1:])
    return paths_dict


def write_long_hierarchy(hierarchy_paths, output):
    for super_class in hierarchy_paths:
        for reaction_path in hierarchy_paths[super_class]:
            output.write(super_class + "\t" + "\t".join(reaction_path) + "\n")
    return


def prune_reactions(hierarchy_paths):
    pathways_dict = dict()
    for super_class in hierarchy_paths:
        if super_class not in pathways_dict.keys():
            pathways_dict[super_class] = list()
        for reaction_path in hierarchy_paths[super_class]:
            pathway_path = reaction_path[0:-1]
            if pathway_path not in pathways_dict[super_class]:
                pathways_dict[super_class].append(pathway_path)
    return pathways_dict


def is_annotated(fields, metacyc=False):
    # orf, contig = fields[0:2]
    annotations = fields[2:]
    for entry in annotations:
        if entry != "" and entry != "\n":
            if metacyc and fields[-1] != "\n":
                return True
            else:
                return False


def get_process_paths(comp_list, hierarchy_paths):
    # Read the list, orf annotations and pathway annotations
    components = list()
    with open(comp_list) as components_list:
        for comp in components_list:
            # print comp
            if "\t" in comp:
                code, common = comp.strip().split("\t")
            else:
                code = comp.strip()
                common = ""
            components.append(code)

    # Map the processes in components to the hierarchy
    process_hierarchy_map = dict()
    for process in components:
        process_hierarchy_map[process] = list()
        for sub_class in hierarchy_paths.values():
            for path in sub_class:
                if process in path:
                    start = path.index(process)
                    process_hierarchy_map[process].append(path[start+1:])

    # Remove processes with no entry in the hierarchy
    no_hierarchy = list()
    for process in process_hierarchy_map:
        if len(process_hierarchy_map[process]) == 0:
            print "No entries found for", process, "in hierarchy"
            no_hierarchy.append(process)
    for process in no_hierarchy:
        process_hierarchy_map.pop(process)
    return process_hierarchy_map


def write_process_annotations(orf_annotations, output, num_processes):
    output.write("Process\tSample\tPathway\tPathway Common Name\tNumber of reactions"
                 "\tNumber covered\tReaction\tORF name\n")
    annotations_found = 0
    for process in orf_annotations.keys():
        if len(orf_annotations[process]) == 0:
            print "No entries found for", process, "in annotations"
        elif len(orf_annotations[process]) > 0:
            annotations_found += 1
            previous_annotation = ""
            process_annotations = sorted(orf_annotations[process])
            for annotation in process_annotations:
                if annotation != previous_annotation:
                    output.write(process + "\t" + "\t".join(annotation) + "\n")
                previous_annotation = annotation
    print annotations_found, "/", num_processes, "processes mapped to annotations"
    return


def search_dict(query_one, query_two, dictionary):
    """
    Matches a word with a key in the dictionary
    :param word: A pathway name
    :param dictionary: A collection of process names (which are anything in the metacyc hierarchy) and all
    processes associated with it from the hierarchy
    :return: a list packing the process and path associated with word or nothing
    """
    process_annotations = list()

    for process, path in dictionary.items():
        for pairs in path:
            if query_two in pairs:
                process_annotations.append([process, pairs[0]])
            if query_one in pairs:
                process_annotations.append([process, path[0]])

    if len(process_annotations) == 0:
        return False
    else:
        return process_annotations


def read_pwy_txt(pathway_annotations, process_hierarchy_map):
    orf_annotations = dict()
    with open(pathway_annotations) as pwy_txt:
        line = pwy_txt.readline()
        while line:
            fields = line.split("\t")
            if not fields[0] == "SAMPLE":
                # Match the pathway name with the process of interest, if possible
                process_annotations = search_dict(fields[1], fields[3], process_hierarchy_map)
                if process_annotations:
                    for path in process_annotations:
                        process = path[0]
                        if process not in orf_annotations.keys():
                            orf_annotations[process] = list()
                        orf = fields[-1].strip()
                        sample = fields[0]
                        n_reactions = fields[5]
                        n_covered = fields[6]
                        pwy_name = fields[1]
                        pwy_common = fields[2]
                        rxn_name = fields[3]
                        orf_annotations[process].append([sample, pwy_name, pwy_common,
                                                         n_reactions, n_covered, rxn_name, orf])
            line = pwy_txt.readline()

    # Add the missing processes to the orf_annotations with empty lists
    for process in process_hierarchy_map:
        if process not in orf_annotations.keys():
            orf_annotations[process] = list()
    return orf_annotations


def main():
    args = get_options()
    check_inputs(args)

    hierarchy = open(args.hierarchy, 'r')
    hierarchy_paths = read_hierarchy(hierarchy)
    hierarchy.close()
    if args.pathways:
        hierarchy_paths = prune_reactions(hierarchy_paths)

    try:
        output = open(args.output, 'w')
    except:
        raise IOError("ERROR: cannot make " + args.output + "!")
    if args.comp_list is None:
        write_long_hierarchy(hierarchy_paths, output)
    else:
        process_hierarchy_map = get_process_paths(args.comp_list, hierarchy_paths)
        orf_annotations = read_pwy_txt(args.pathway_annotations, process_hierarchy_map)
        write_process_annotations(orf_annotations, output, len(process_hierarchy_map.keys()))
    output.close()
main()
