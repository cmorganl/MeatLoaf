import sys
import argparse
import re
import itertools
import numpy as np
import unittest

from Bio.SeqFeature import SeqFeature, FeatureLocation

_DIST_LIMIT = 30


class GBFTTestClass(unittest.TestCase):
    def setUp(self) -> None:
        self.test_data_dir = "./test_data/"
        self.original_dove = "{}AY714826.1.tbl".format(self.test_data_dir)
        self.update_dove = "{}AY714826.tbl".format(self.test_data_dir)
        self.original_ragged = "{}AY714854.1.tbl".format(self.test_data_dir)
        self.update_ragged = "{}AY714854.tbl".format(self.test_data_dir)
        self.original_miss = "{}AY714835.1.tbl".format(self.test_data_dir)
        self.update_miss = "{}AY714835.tbl".format(self.test_data_dir)
        self.original_sub = "{}AY714824.1.tbl".format(self.test_data_dir)
        self.update_sub = "{}AY714824.tbl".format(self.test_data_dir)

        self.test_arrays = (['1', '2', '3', '400', '2'], ['2', '3', '400'])
        self.real_arrays = ([1277, 1277, -839, -839, -1388, -1388, -1616, -1616, -1283, -1283, -851, -851,
                             830, 830, -365, -365, -692, -692, -818, -818, -464, -464, 1541, 1541, -73, -73,
                             -974, -974, 533, 533, 650, 650, 758, 758, 920, 920],
                            [293, 293, 296, 296, 152, 152,
                             1583, 1583, -839, -839, -1388, -1388, -1616, -1616, -1283, -1283, -851, -851,
                             821, 821, -365, -365, -737, -737, -818, -818, -464, -464, 1541, 1541, -75, -75,
                             533, 533, 650, 650, 758, 758, 920, 920])
        return

    def test_parse_features_from_table(self):
        _, loci = parse_loci_from_table(feature_table=self.update_sub)
        self.assertEqual(41, len(loci))
        self.assertEqual(2, set([lc.n_features() for lc in loci]).pop())
        return

    def test_count_loci(self):
        _, og_loci = parse_loci_from_table(feature_table=self.original_sub)
        self.assertEqual(18, count_loci(og_loci))

    def test_reconcile_feature_lists(self):
        combined = reconcile_feature_lists(locus_list_one=parse_loci_from_table(feature_table=self.original_sub)[1],
                                           locus_list_two=parse_loci_from_table(feature_table=self.update_sub)[1])
        self.assertEqual(82, len(combined))

    def test_gbft_merge_subseq(self):
        main("--ft_ref {0} --ft_two {1} --prefix {2}".format(self.original_sub, self.update_sub, "GZ17G11").split())
        return

    def test_gbft_merge_dovetail(self):
        main("--ft_ref {0} --ft_two {1} --prefix {2}".format(self.original_dove, self.update_dove, "GZ18C8").split())
        return

    def test_gbft_merge_ragged(self):
        main("--ft_ref {0} --ft_two {1} --prefix {2}".format(self.original_ragged, self.update_ragged, "GZ").split())
        return
    
    def test_gbft_merge_miss(self):
        main("--ft_ref {0} --ft_two {1} --prefix {2}".format(self.original_miss, self.update_miss, "GZ").split())
        return


class Locus:
    def __init__(self, start, stop, strand):
        self.start = start
        self.stop = stop
        self.strand = strand
        self.features = []
        self.tag = ""

    def n_features(self) -> int:
        return len(self.features)

    def get_info(self) -> str:
        return "Locus tag:\t'{}'\nStart-Stop\t{}-{}\nStrand\t{}\n".format(self.tag, self.start, self.stop, self.strand)

    def locus_len(self) -> int:
        length = abs_dist(self.stop, self.start) * self.strand
        if length != 0:
            return length
        else:
            exit_gracefully("A locus of length 0 was encountered:\n{}\n".format(self.get_info()))

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
    h = matrix(a, b, match_score, gap_cost)
    b_, pos = traceback(h, b)
    return pos, pos + len(b_)
# End


def get_arg_parser(sys_args):
    parser = argparse.ArgumentParser(add_help=False,
                                     description="The GenBank Feature Table (GBFT) merging script. "
                                                 "This can be used for comparing the features of two feature tables, "
                                                 "a 'reference' table and a 'new' table and merge the feature IDs "
                                                 "found in the reference table while maintaining all the new feature "
                                                 "table's locations and features that are not found in the reference.")
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


def abs_dist(n1: int, n2: int) -> int:
    return max(n1, n2) - min(n1, n2)


def get_feature_coordinates(pos_fields) -> (int, int, int):
    fp, sp = [int(re.sub(r'[<>]', '', i)) for i in pos_fields]
    if fp > sp:
        start = sp
        end = fp
        strand = -1
    else:
        start = fp
        end = sp
        strand = 1
    return start, end, strand


def load_lines_to_feature(feat_tbl_lines: list, sep="\t") -> SeqFeature:
    # Load the feature's description and positions
    desc_line = feat_tbl_lines.pop(0)
    fields = desc_line.split(sep)
    start, stop, strand = get_feature_coordinates(fields[:2])
    seq_feat = SeqFeature(type=fields[2], location=FeatureLocation(start, stop), strand=strand)
    seq_feat.id = None

    # Load the feature's qualifier lines
    while feat_tbl_lines and len(feat_tbl_lines[0].split(sep)) == 5:
        fields = feat_tbl_lines.pop(0).split(sep)
        qualifier, value = fields[3:]
        seq_feat.qualifiers[qualifier] = value
        if qualifier == "locus_tag":
            seq_feat.id = value
    return seq_feat


def load_lines_to_locus(feat_tbl_lines: list, sep="\t") -> Locus:
    """
    Loads a list of lines into a Locus instance, the features attribute of which contains a list of SeqFeature objects.

    :param feat_tbl_lines: A list of lines following the format of a valid GenBank feature table
    :param sep: The field separator, default being tab. This should never change, really.
    :return: A locus instance
    """

    # Load the feature's description and positions
    fields = feat_tbl_lines[0].split(sep)
    start, stop, strand = get_feature_coordinates(fields[:2])
    locus = Locus(start, stop, strand)

    # While the start, stop and strand values remain unchanged, continue loading new features into locus
    while feat_tbl_lines:
        if len(feat_tbl_lines[0]) == 0:
            break
        fields = feat_tbl_lines[0].split(sep)

        coords = get_feature_coordinates(fields[:2])
        if start != coords[0] or stop != coords[1]:
            break
        locus.features.append(load_lines_to_feature(feat_tbl_lines))

    return locus


def parse_loci_from_table(feature_table: str) -> (str, list):
    header = ""
    loci = []
    with open(feature_table) as ft_handler:
        header_line = ft_handler.readline()
        if header_line[0] != '>':
            exit_gracefully("Error: first line in feature table file isn't formatted properly:\n{}".format(header_line))
        else:
            header = header_line.strip()
        feature_lines = [line.rstrip() for line in ft_handler.readlines()]

    line = feature_lines[0]
    while line:
        loci.append(load_lines_to_locus(feature_lines))
        try:
            line = feature_lines[0]
        except IndexError:
            break
        if line and len(line.split("\t")) != 3:
            exit_gracefully("Unexpected line format in feature table '{}':\n{}".format(feature_table, line))

    return header, loci


def get_feature_length(feature: SeqFeature) -> str:
    return (feature.location.end - feature.location.start) * feature.strand


def validate_feature_merge(ref_feature: SeqFeature, putative_feature: SeqFeature) -> None:
    # Check to be sure the strand and type are identical
    if ref_feature.type != putative_feature.type:
        exit_gracefully("SeqFeature type mismatch during merge:\n{}\n{}".format(ref_feature, putative_feature))
    elif ref_feature.strand != putative_feature.strand:
        exit_gracefully("SeqFeature strand mismatch during merge:\n{}\n{}".format(ref_feature, putative_feature))
    return


def merge_features_from_locus(ref_locus: Locus, putative_locus: Locus) -> list:
    merged_features = []
    while ref_locus.features and putative_locus.features:
        ref_feature = ref_locus.features.pop(0)
        putative_feature = putative_locus.features.pop(0)
        validate_feature_merge(ref_feature, putative_feature)
        new_feature = SeqFeature(type=ref_feature.type, location=putative_feature.location, strand=ref_locus.strand)
        for qual, val in putative_feature.qualifiers.items():
            # Populate the feature qualifiers, giving priority to the ref_feature values
            if qual not in ref_feature.qualifiers:
                new_feature.qualifiers[qual] = val
            else:
                new_feature.qualifiers[qual] = ref_feature.qualifiers[qual]
        merged_features.append(new_feature)

    return merged_features


def align_start_positions(feat_one_lens, feat_two_lens, feat_list_one, feat_list_two) -> (list, int, int):
    """
    Ensure that the beginning first SeqFeature instance in each of the two lists of SeqFeatures are the same length.

    :param feat_one_lens:
    :param feat_two_lens:
    :param feat_list_one:
    :param feat_list_two:
    :return: A list of SeqFeature instances from feat_list_two (the new features) that were popped to align the lists
    """
    updated_features = []
    skipped_features = ""

    start_one, _ = smith_waterman(feat_two_lens, feat_one_lens)
    start_two, _ = smith_waterman(feat_one_lens, feat_two_lens)

    while start_one or start_two:
        if min(len(feat_one_lens), len(feat_two_lens)) == 0 or abs_dist(feat_one_lens[0], feat_two_lens[0]) < _DIST_LIMIT:
            break
        if start_one > 0:
            skipped_locus = feat_list_one.pop(0)
            skipped_features += locus_features_to_string(skipped_locus) + "\n"
            feat_one_lens.pop(0)
            start_one -= 1
        if start_two > 0:
            skipped_locus = feat_list_two.pop(0)  # type: Locus
            for seq_feature in skipped_locus.features:
                updated_features.append(seq_feature)
            feat_two_lens.pop(0)
            start_two -= 1

    if skipped_features:
        print("Reference feature(s) skipped:\n{}".format(skipped_features))

    try:
        return updated_features, feat_one_lens.pop(0), feat_two_lens.pop(0)
    except IndexError:
        return updated_features, 0, 0


def reconcile_feature_lists(locus_list_one: list, locus_list_two: list) -> list:
    """
    With feat_list_one representing the current reference features (i.e. those with the identifiers to be maintained)
    find the overlapping features, by adding missing qualifiers from feat_list_two and augmenting the former reference
    set with new features from feat_list_two.

    :param locus_list_one: The reference list of Locus instances, containing their corresponding SeqFeature instances
    :param locus_list_two:
     Used to prevent redundant locus identifiers by adding to the existing 'locus_id' not present in reference features.
    :return: A list of SeqFeature instances
    """
    updated_features = []
    # Determine the orientation of the two feature lists by comparing the lengths of their features

    feat_one_lens = [lc.locus_len() for lc in locus_list_one]
    feat_two_lens = [lc.locus_len() for lc in locus_list_two]
    desired_features = sum([lc.n_features() for lc in locus_list_two])

    while feat_one_lens and feat_two_lens:
        ref_len = feat_one_lens.pop(0)
        new_len = feat_two_lens.pop(0)
        # Ensure the feature lengths are similar between the old and reference features, ORF discrepancies are abundant
        if abs_dist(ref_len, new_len) > _DIST_LIMIT:
            feat_one_lens = [ref_len] + feat_one_lens
            feat_two_lens = [new_len] + feat_two_lens
            skipped_features, ref_len, new_len = align_start_positions(feat_one_lens, feat_two_lens,
                                                                       locus_list_one, locus_list_two)
            if len(skipped_features) == desired_features:
                print("No matching loci were found between the original and new feature tables.")

            updated_features += skipped_features
        if ref_len == new_len == 0:
            break
        updated_features += merge_features_from_locus(locus_list_one.pop(0), locus_list_two.pop(0))

    # Ensure all features remaining are removed from both feature lists
    while locus_list_one:
        print("Reference sequence feature skipped:\n{}".format(locus_features_to_string(locus_list_one.pop(0))))
    while locus_list_two:
        updated_features.append(locus_list_two.pop(0))

    # Test the difference between the original number of features and the updated number
    if len(updated_features) != desired_features:
        exit_gracefully("Number of updated features ({}) is different from the number of new features ({})."
                        "".format(len(updated_features), desired_features))

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


def locus_features_to_string(locus: Locus) -> str:
    feature_str = ""
    for seq_feature in locus.features:  # type: SeqFeature
        feature_str += seq_feature_to_string(seq_feature)
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


def count_loci(loci: list, sep='_') -> int:
    """
    When merging two feature tables, the majority of the features will be common between them.
    However, there may be some features that are not shared between them and if they are only present in the newer table
    they will need to be assigned new identifiers.

    This function parses the original feature table (any, really but the original is the one that should be passed in)
    and determines the maximum feature table value (i.e. the last feature's number).
    When new features are encountered, their identifiers will begin incrementing from this number.

    :param loci: A list of Locus instances
    :param sep: A character separating the feature table's prefix string and the feature's unique numeric identifier
    in the SeqFeature.id attribute.
    :return: An integer value representing the last feature
    """
    acc = 0
    num = 0
    for locus in loci:  # type: Locus
        for seq_feat in locus.features:  # type: SeqFeature
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

    ref_head, ref_loci = parse_loci_from_table(args.ft_ref)
    locus_acc_offset = count_loci(ref_loci)
    _, two_loci = parse_loci_from_table(args.ft_two)

    updated_features = reconcile_feature_lists(ref_loci, two_loci)

    if args.prefix:
        fix_feature_tags(updated_features, args.prefix, locus_acc_offset)

    write_feature_table(args.output_ft, ref_head, updated_features)

    return


if __name__ == "__main__":
    main(sys.argv[1:])
