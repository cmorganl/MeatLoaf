import sys
import argparse
import re
import itertools
import numpy as np
import unittest

from Bio.SeqFeature import SeqFeature, FeatureLocation


class GBFTTestClass(unittest.TestCase):
    def setUp(self) -> None:
        self.base_dir = "/mnt/sdb/Hallam_projects/ANME/ProcessedData/NCBI_deposition"
        self.test_arrays = (['1', '2', '3', '400', '2'], ['2', '3', '400'])
        self.real_arrays = ([1277, 1277, -839, -839, -1388, -1388, -1616, -1616, -1283, -1283, -851, -851,
                             830, 830, -365, -365, -692, -692, -818, -818, -464, -464, 1541, 1541, -73, -73,
                             -974, -974, 533, 533, 650, 650, 758, 758, 920, 920],
                            [293, 293, 296, 296, 152, 152,
                             1583, 1583, -839, -839, -1388, -1388, -1616, -1616, -1283, -1283, -851, -851,
                             821, 821, -365, -365, -737, -737, -818, -818, -464, -464, 1541, 1541, -75, -75,
                             533, 533, 650, 650, 758, 758, 920, 920])
        return

    def test_gbft_merge(self):
        main("--ft_ref {0}/fosmid_resequence_update/AY714824.1.tbl "
             "--ft_two {0}/fosmid_resequence_update/Prokka/GZfos17G11_3436005/AY714824.tbl "
             "--prefix {1}".format(self.base_dir, "GZ17G11").split())
        return

    def test_smith_waterman(self):
        start, end = smith_waterman([str(i) for i in self.real_arrays[1]], [str(i) for i in self.real_arrays[0]])
        self.assertEqual(2, start)
        return


# The following code was copied directly from Daniel Tiefenauer's blog at:
# https://tiefenauer.github.io/blog/smith-waterman/#usage-and-tests

def matrix(a, b, match_score=3, gap_cost=2):
    H = np.zeros((len(a) + 1, len(b) + 1), np.int)

    for i, j in itertools.product(range(1, H.shape[0]), range(1, H.shape[1])):
        match = H[i - 1, j - 1] + (match_score if a[i - 1] == b[j - 1] else - match_score)
        delete = H[i - 1, j] - gap_cost
        insert = H[i, j - 1] - gap_cost
        H[i, j] = max(match, delete, insert, 0)
    return H


def traceback(H, b, b_='', old_i=0):
    # flip H to get index of **last** occurrence of H.max() with np.argmax()
    H_flip = np.flip(np.flip(H, 0), 1)
    i_, j_ = np.unravel_index(H_flip.argmax(), H_flip.shape)
    i, j = np.subtract(H.shape, (i_ + 1, j_ + 1))  # (i, j) are **last** indexes of H.max()
    if H[i, j] == 0:
        return b_, j
    b_ = b[j - 1] + '-' + b_ if old_i - i > 1 else b[j - 1] + b_
    return traceback(H[0:i, 0:j], b, b_, i)


def smith_waterman(a, b, match_score=1, gap_cost=2):
    a, b = [str(i) for i in a], [str(j) for j in b]
    H = matrix(a, b, match_score, gap_cost)
    b_, pos = traceback(H, b)
    return pos, pos + len(b_)
# End


def get_arg_parser(sys_args):
    parser = argparse.ArgumentParser(add_help=False)
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument("--ft_ref", required=True,
                               help="The first feature table, whose 'locus_tag's, 'protein_id's and other qualifiers "
                                    "will be transfered onto the other feature table")
    required_args.add_argument("--ft_two", required=True,
                               help="The new feature table that will receive the current identifiers from ft_ref")

    optopt = parser.add_argument_group("Optional arguments")
    optopt.add_argument("-o", "--output", dest="output_ft", default="genbank_features.tbl",
                        help="Path to a new feature table to write to",
                        required=False)
    optopt.add_argument("-f", "--fasta", dest="fasta",
                        help="Path to a fasta file representing the new feature table",
                        required=False)
    optopt.add_argument("-p", "--prefix", default=None,
                        help="The prefix to use for the 'locus_tag' and 'protein_id' feature qualifiers")

    miscellaneous_opts = parser.add_argument_group("Miscellaneous options")
    miscellaneous_opts.add_argument('--overwrite', action='store_true', default=False,
                                    help='overwrites previously processed output folders')
    miscellaneous_opts.add_argument("-h", "--help",
                                    action="help",
                                    help="show this help message and exit")
    args = parser.parse_args(sys_args)

    return args


def exit_gracefully(exit_statement: str, errorcode=1) -> SystemExit:
    print(exit_statement)
    return sys.exit(errorcode)


def parse_features_from_table(feature_table: str) -> (str, list):
    header = ""
    features = []
    with open(feature_table) as ft_handler:
        header_line = ft_handler.readline()
        if header_line[0] != '>':
            exit_gracefully("Error: first line in feature table file isn't formatted properly:\n{}".format(header_line))
        else:
            header = header_line.strip()
        for line in [l.rstrip() for l in ft_handler.readlines()]:
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) == 3:
                fp, sp = [int(i) for i in fields[:2]]
                if fp > sp:
                    start = sp
                    end = fp
                    strand = -1
                else:
                    start = fp
                    end = sp
                    strand = 1
                feature = SeqFeature(type=fields[2], location=FeatureLocation(start, end), strand=strand)
                feature.id = None
                features.append(feature)
            elif len(fields) == 5:
                qualifier, value = fields[3:]
                feature.qualifiers[qualifier] = value
                if qualifier == "locus_tag":
                    feature.id = value
            else:
                exit_gracefully("Unexpected line format in feature table '{}':\n{}".format(feature_table, line))

    return header, features


def get_feature_length(feature: SeqFeature) -> str:
    return (feature.location.end - feature.location.start) * feature.strand


def merge_qualifiers_from_feature(ref_feature: SeqFeature, putative_feature: SeqFeature) -> SeqFeature:
    # Check to be sure the strand and type are identical
    if ref_feature.type != putative_feature.type:
        exit_gracefully("SeqFeature type mismatch during merge:\n{}\n{}".format(ref_feature, putative_feature))
    elif ref_feature.strand != putative_feature.strand:
        exit_gracefully("SeqFeature strand mismatch during merge:\n{}\n{}".format(ref_feature, putative_feature))

    new_feature = SeqFeature(type=ref_feature.type, location=putative_feature.location, strand=putative_feature.strand)
    for qual, val in putative_feature.qualifiers.items():
        # Populate the feature qualifiers, giving priority to the ref_feature values
        if qual not in ref_feature.qualifiers:
            new_feature.qualifiers[qual] = val
        else:
            new_feature.qualifiers[qual] = ref_feature.qualifiers[qual]

    return new_feature


def reconcile_feature_lists(feat_list_one: list, feat_list_two: list) -> list:
    """
    With feat_list_one representing the current reference features (i.e. those with the identifiers to be maintained)
    find the overlapping features, by adding missing qualifiers from feat_list_two and augmenting the former reference
    set with new features from feat_list_two.

    :param feat_list_one: The reference list of features
    :param feat_list_two:
     Used to prevent redundant locus identifiers by adding to the existing 'locus_id' not present in reference features.
    :return: A list of SeqFeature instances
    """
    updated_features = []
    # Determine the orientation of the two feature lists by comparing the lengths of their features

    feat_one_lens = [get_feature_length(f) for f in feat_list_one]
    feat_two_lens = [get_feature_length(f) for f in feat_list_two]

    start_one, _ = smith_waterman(feat_two_lens, feat_one_lens)
    start_two, _ = smith_waterman(feat_one_lens, feat_two_lens)

    i = 0
    j = 0
    while start_one > 0:
        print("Reference sequence feature skipped:\n{}".format(seq_feature_to_string(feat_list_one.pop(0))))
        start_one -= 1
        i += 1
    while start_two > 0:
        updated_features.append(feat_list_two.pop(0))
        start_two -= 1
        j += 1

    while feat_list_one and feat_list_two:
        while not -50 < (feat_one_lens[i] - feat_two_lens[j]) < 50:
            print("Reference sequence feature skipped:\n{}".format(seq_feature_to_string(feat_list_one.pop(0))))
            i += 1
        updated_features.append(merge_qualifiers_from_feature(feat_list_one.pop(0),
                                                              feat_list_two.pop(0)))
        i += 1
        j += 1

    # Ensure all features remaining are removed from both feature lists
    while feat_list_one:
        print("Reference sequence feature skipped:\n{}".format(seq_feature_to_string(feat_list_one.pop(0))))
    while feat_list_two:
        updated_features.append(feat_list_two.pop(0))

    if len(updated_features) != len(feat_two_lens):
        exit_gracefully("Number of updated features ({}) is different from the number of new features ({})."
                        "".format(len(updated_features), len(feat_two_lens)))

    return updated_features


def seq_feature_to_string(seq_feature: SeqFeature) -> str:
    qual_order = ["EC_number", "gene", "locus_tag",
                  "inference", "product", "transl_tabl", "protein_id", "note"]
    feature_str = ""
    if seq_feature.strand == 1:
        start = seq_feature.location.start
        end = seq_feature.location.end
    else:
        start = seq_feature.location.end
        end = seq_feature.location.start
    feature_str += "{}\t{}\t{}\n".format(start, end, seq_feature.type)
    for qual in qual_order:
        try:
            feature_str += "\t\t\t{}\t{}\n".format(qual, seq_feature.qualifiers[qual])
        except KeyError:
            continue
    return feature_str


def fix_feature_tags(seq_features: list, locus_prefix: str, ref_locus_offset=0) -> None:
    """

    :param seq_features:
    :param locus_prefix:
    :param ref_locus_offset: The number of loci in the reference feature table. Default = 0 (no offset).
    :return:
    """
    for feat in seq_features:  # type: SeqFeature
        for qual_name in ["locus_tag", "protein_id"]:
            try:
                name = feat.qualifiers[qual_name]
                if re.search(locus_prefix, name):
                    continue
                fields = name.split('_')
                # This could be a GenBank accession
                if len(fields) == 1:
                    continue
                orf_num = int(fields[-1]) + ref_locus_offset
                feat.qualifiers[qual_name] = "{}_{}".format(locus_prefix, orf_num)
            except KeyError:
                continue
    return


def write_feature_table(output_ft: str, header: str, feature_list: list) -> None:
    with open(output_ft, 'w') as out_handler:
        out_handler.write(header + "\n")
        for feature in feature_list:  # type: SeqFeature
            out_handler.write(seq_feature_to_string(feature))

    return


def count_loci(seq_features: list, sep='_') -> int:
    acc = 0
    for seq_feat in seq_features:  # type: SeqFeature
        if not seq_feat.id:
            continue
        # Split the locus_id into it's component parts
        parts = seq_feat.id.split(sep)
        try:
            num = int(parts[-1])
        except TypeError:
            exit_gracefully("Error: Unexpected format of locus identifier ({}). "
                            "Last element assumed to be an integer.".format(seq_feat))
        if num > acc:
            acc = num
    return acc


def main(sys_args):
    args = get_arg_parser(sys_args)

    ref_head, ref_features = parse_features_from_table(args.ft_ref)
    locus_acc_offset = count_loci(ref_features)
    two_head, two_features = parse_features_from_table(args.ft_two)

    updated_features = reconcile_feature_lists(ref_features, two_features)

    if args.prefix:
        fix_feature_tags(updated_features, args.prefix, locus_acc_offset)

    write_feature_table(args.output_ft, ref_head, updated_features)

    return


if __name__ == "__main__":
    main(sys.argv[1:])
