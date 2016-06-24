#!/usr/bin/env python

import os
import argparse


def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--hierarchy", required=True,
                        help="The hierarchy text file for mapping components to classes")
    parser.add_argument("-o", "--output", required=True,
                        help="The output file containing mapped values from the hierarchy")
    parser.add_argument("-l", "--comp_list", required=False, default=None,
                        help="An optional list of hierarchy components to return information on")
    parser.add_argument("--pathways", required=False, default=False, action="store_true",
                        help="Flag to indicate the output should be at pathway level [DEFAULT = reaction-level]")
    # parser.add_argument("-a", "--orf_annotations", required=False,
    #                     help="If using a list to determine the proteins captured within pathways,"
    #                          "an associated MetaPathways ORF_annotation_table is provided here")
    parser.add_argument("-p", "--pathway_annotations", required=True,
                        help="The pwy.txt file from metapathways")
    parser.add_argument("-f", "--pgdb_input", required=True,
                        help="The 0.pf file from metapathways")
    args = parser.parse_args()

    if not os.path.isfile(args.hierarchy):
        raise IOError("ERROR: " + args.hierarchy + " doesn't exist!")

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


def map_proteins(process_hierarchy_map, ptools_input, orf_annotations):
    # Map the reactions and/or pathways for each process in process_hierarchy_dict to the ORF annotations
    print "Annotations:"
    annotation_process_map = dict()
    for process, paths in process_hierarchy_map.items():
        annotation_process_map[process] = list()
        for path in paths:
            for part in path:
                for name, product in ptools_input.items():
                    if part in product:
                        print part, product
                        annotation_process_map[process].append(product)
                    elif name in orf_annotations.keys():
                        print part, product
    annotations_found = 0
    for process in annotation_process_map:
        if len(annotation_process_map[process]) == 0:
            print "No entries found for", process, "in annotations"
        else:
            annotations_found += 1
    print annotations_found, "/", len(annotation_process_map.keys()), "processes mapped to annotations"
    return annotation_process_map


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
                    process_hierarchy_map[process].append(path)

    # Remove processes with no entry in the hierarchy
    no_hierarchy = list()
    for process in process_hierarchy_map:
        if len(process_hierarchy_map[process]) == 0:
            print "No entries found for", process, "in hierarchy"
            no_hierarchy.append(process)
    for process in no_hierarchy:
        process_hierarchy_map.pop(process)
    return process_hierarchy_map


def write_process_annotations(annotation_process_map, output):
    output.write("Process\tProtein annotation\tCommon name\tNumber of reactions\tNumber covered\n")
    for process in annotation_process_map.keys():
        if len(annotation_process_map[process]) > 0:
            for annotation in annotation_process_map[process]:
                output.write(process + "\t" + annotation + "\n")
    return


def read_pgdb_input(pgdb_input):
    ptools_input = dict()
    name = ""
    with open(pgdb_input) as pf:
        for line in pf:
            if line.strip() == "//":
                pass
            else:
                key, value = line.strip().split("\t")
                if key == "NAME":
                    name = value
                if key == "PRODUCT":
                    ptools_input[name] = value
    return ptools_input


def search_dict(word, dictionary):
    for key in dictionary:
        if word == key:
            return dictionary[key]
        for value in dictionary[key]:
            if word in value:
                return value
    return False


def read_pwy_txt(pathway_annotations, process_hierarchy_map):
    orf_annotations = dict()
    with open(pathway_annotations) as pwy_txt:
        for line in pwy_txt:
            fields = line.split("\t")
            if not fields[0] == "SAMPLE":
                status = search_dict(fields[1], process_hierarchy_map)
                if status:
                    orfs = fields[-1][1:-2].split(',')
                    for orf in orfs:
                        orf_annotations[orf] = status
    return orf_annotations


def main():
    args = get_options()
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
        ptools_input = read_pgdb_input(args.pgdb_input)
        process_hierarchy_map = get_process_paths(args.comp_list, hierarchy_paths)
        orf_annotations = read_pwy_txt(args.pathway_annotations, process_hierarchy_map)
        annotation_process_map = map_proteins(process_hierarchy_map, ptools_input, orf_annotations)
        write_process_annotations(annotation_process_map, output)
    output.close()
main()
