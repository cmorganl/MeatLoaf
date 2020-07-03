#!/usr/bin/env python

import os
import sys
import json
import logging
import argparse

_HIERARCHY = {0: '', 1: '', 2: '', 3: ''}


class KoId:
    def __init__(self, ko_record: str):
        words = ko_record.split()
        self.ko = words[0]
        self.desc = ' '.join(words[1:]).strip()
        self.a = ""
        self.b = ""
        self.c = ""
        self.d = ""
        return

    def to_tsv(self):
        return "\t".join([self.ko, self.desc, self.a, self.b, self.c])

    def add_hierachy_data(self):
        self.a = _HIERARCHY[1]
        self.b = _HIERARCHY[2]
        self.c = _HIERARCHY[3]


def get_arguments():
    parser = argparse.ArgumentParser(description="Converts a KO hierarchy JSON file (ko0001) to a table mapping "
                                                 "unique KEGG Orthology IDs to each level in the hierarchy. "
                                                 "The input file can be downloaded from here:\n"
                                                 "https://www.genome.jp/kegg-bin/get_htext?ko00001.keg")
    parser.add_argument("-k", "--hierarchy", help="KEGG Hierarchy file in JSON format", required=True)
    args = parser.parse_args()
    return args


def read_hierarchy(ko_json):
    if not os.path.isfile(ko_json):
        logging.error("File '{}' doesn't exist.\n".format(ko_json))
        sys.exit(3)
    try:
        ko_h = open(ko_json, 'r')
    except IOError:
        logging.error("Unable to open '{}' for reading.\n".format(ko_json))
        sys.exit(5)

    ko_dict = json.load(ko_h)
    ko_h.close()
    return ko_dict


def get_leaves(item, level=0) -> list:
    _HIERARCHY[level] = ' '.join(item["name"].split()[1:])
    level += 1
    leaves = []
    try:
        for i in item["children"]:
            leaves.extend(get_leaves(i, level))
    except KeyError:
        ko = KoId(item["name"])
        ko.add_hierachy_data()
        return [ko]
    return leaves


def load_ko(ko_dict: dict) -> list:
    ko_list = get_leaves(ko_dict)
    return ko_list


def write_table(ko_data: list, output_file: str) -> None:
    header = "\t".join(["KO", "Description", "L1", "L2", "L3"]) + "\n"
    try:
        handler = open(output_file, 'w')
    except IOError:
        logging.error("Unable to open output file '{}' for writing.\n".format(output_file))
        sys.exit(7)

    handler.write(header)
    for ko in ko_data:  # type: KoId
        handler.write(ko.to_tsv() + "\n")
    handler.close()
    return


def main():
    args = get_arguments()
    path, ext = os.path.splitext(args.hierarchy)
    output_tbl = path + ".tsv"
    ko_dict = read_hierarchy(args.hierarchy)
    ko_data = load_ko(ko_dict)
    write_table(ko_data, output_tbl)


main()
