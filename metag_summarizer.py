#!/usr/bin/env python3

import os
import csv
import unittest
import sys
import logging
import argparse

import numpy as np
import matplotlib.pyplot as plt
from pyfastx import Fasta, Fastq
from pyfastxcli import fastx_format_check

# from sklearn import svm, datasets
# from sklearn.metrics import roc_curve, auc
# from sklearn.model_selection import train_test_split
# from sklearn.preprocessing import label_binarize
# from sklearn.multiclass import OneVsRestClassifier
# from sklearn.metrics import roc_auc_score


class MetagenomeROCTester(unittest.TestCase):
    def setUp(self) -> None:
        self.bracken_tbl = os.path.join("test_data", "ZCS100k_levelled.bracken")
        self.positive_taxa = os.path.join("test_data", "taxa_names.txt")
        self.test_fq = os.path.join("test_data", "test_TarA.1.fq")
        self.taxa_abund = {"Salmonella enterica": 1000, "Staphylococcus aureus": 500, "Pseudomonas": 500}
        return

    def test_read_bracken_classifications(self):
        taxa_map = read_bracken_classifications(self.bracken_tbl)
        self.assertEqual(1076, len(taxa_map))
        self.assertEqual(758, taxa_map["Salmonella bongori"])
        return

    def test_read_lines_to_set(self):
        species = read_lines_to_set(self.positive_taxa)
        self.assertEqual(10, len(species))
        return

    def test_num_fastx_records(self):
        num_reads = get_num_fastx_records(self.test_fq)
        self.assertEqual(12, num_reads)

    def test_bin_classes(self):
        tp, fp, tn, fn = bin_classes(taxa_abund=self.taxa_abund, total_queries=3000,
                                     positives={"Salmonella enterica", "Pseudomonas"})
        self.assertEqual(1500, tp)


class MyFormatter(logging.Formatter):

    error_fmt = "%(levelname)s - %(module)s, line %(lineno)d:\n%(message)s"
    warning_fmt = "%(levelname)s:\n%(message)s"
    debug_fmt = "%(asctime)s\n%(message)s"
    info_fmt = "%(message)s"

    def __init__(self):
        super().__init__(fmt="%(levelname)s: %(message)s",
                         datefmt="%d/%m %H:%M:%S")

    def format(self, record):

        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._style._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.DEBUG:
            self._style._fmt = MyFormatter.debug_fmt

        elif record.levelno == logging.INFO:
            self._style._fmt = MyFormatter.info_fmt

        elif record.levelno == logging.ERROR:
            self._style._fmt = MyFormatter.error_fmt

        elif record.levelno == logging.WARNING:
            self._style._fmt = MyFormatter.warning_fmt

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._style._fmt = format_orig

        return result


def prep_logging(log_file_name=None, verbosity=False) -> None:
    """
    Allows for multiple file handlers to be added to the root logger, but only a single stream handler.
    The new file handlers must be removed outside of this function explicitly

    :param log_file_name:
    :param verbosity:
    :return: None
    """
    if verbosity:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.INFO

    # Detect whether a handlers are already present and return if true
    logger = logging.getLogger()
    if len(logger.handlers):
        return

    formatter = MyFormatter()
    # Set the console handler normally writing to stdout/stderr
    ch = logging.StreamHandler()
    ch.setLevel(logging_level)
    ch.terminator = ''
    ch.setFormatter(formatter)

    if log_file_name:
        output_dir = os.path.dirname(log_file_name)
        try:
            if output_dir and not os.path.isdir(output_dir):
                os.makedirs(output_dir)
        except (IOError, OSError):
            sys.stderr.write("ERROR: Unable to make directory '" + output_dir + "'.\n")
            sys.exit(3)
        logging.basicConfig(level=logging.DEBUG,
                            filename=log_file_name,
                            filemode='w',
                            datefmt="%d/%m %H:%M:%S",
                            format="%(asctime)s %(levelname)s:\n%(message)s")
        logging.getLogger('').addHandler(ch)
        logging.getLogger('').propagate = False
    else:
        logging.basicConfig(level=logging_level,
                            datefmt="%d/%m %H:%M:%S",
                            format="%(asctime)s %(levelname)s:\n%(message)s")
    return


def get_options(sys_args):
    parser = argparse.ArgumentParser(description="A script for generating Receiver Operating Characteristic (ROC)"
                                                 " curves from metagenome taxonomic classifications.",
                                     add_help=False)

    req_args = parser.add_argument_group("Required arguments")
    opt_args = parser.add_argument_group("Optional arguments")
    mis_args = parser.add_argument_group("Miscellaneous arguments")

    req_args.add_argument("-c", "--classification_tbl", dest="c_tbl", required=True,
                          help="Path to a classification table")
    req_args.add_argument("-p", "--positive_names", dest="p_names", required=True,
                          help="Path to a file listing the names of true positive organism names -"
                               " one line per organism")
    req_args.add_argument("-i", "--fastx_file", dest="fastx", required=True,
                          help="Path to the FASTA or FASTQ file containing the classified sequences")

    opt_args.add_argument('-o', '--output', default='', required=False,
                          help='The output [ DEFAULT =  ]')
    opt_args.add_argument('-r', "--tax_rank", default="species", required=False,
                          help="The taxonomic rank threshold for correct classifications. [ DEFAULT = 'species' ]")

    mis_args.add_argument('--overwrite', action='store_true', default=False,
                          help='overwrites previously processed output folders')
    mis_args.add_argument('-v', '--verbose', action='store_true', default=False,
                          help='prints a more verbose runtime log')
    mis_args.add_argument("-h", "--help",
                          action="help", help="show this help message and exit")

    args = parser.parse_args(sys_args)
    return args


def read_lines_to_set(file_path: str) -> set:
    taxa = set()
    try:
        taxa_file = open(file_path, 'r')
    except IOError:
        logging.error("Unable to open file '{}' for reading.\n".format(file_path))
        sys.exit(3)

    for line in taxa_file:
        if not line:
            continue
        taxa.add(line.strip())
    taxa_file.close()

    return taxa


def guess_tbl_format(tbl_path: str) -> str:
    tbl_format = ""
    return tbl_format


def read_bracken_classifications(bracken_table: str) -> dict:
    # header = ["name", "taxonomy_id", "taxonomy_lvl",
    #           "kraken_assigned_reads", "added_reads", "new_est_reads", "fraction_total_reads"]
    taxa_read_map = {}
    with open(bracken_table) as tbl:
        reader = csv.DictReader(tbl, delimiter="\t")
        for row in reader:
            taxa_read_map[row["name"]] = int(row["new_est_reads"])
    return taxa_read_map


def get_classifications(tbl_path: str, tbl_format=None) -> dict:
    logging.info("Reading classifications from '{}'... ".format(tbl_path))
    # TODO: Determine the file format is none was provided
    if tbl_format is None:
        tbl_format = guess_tbl_format(tbl_path)
    if tbl_format == "bracken":
        taxa_map = read_bracken_classifications(tbl_path)
    else:
        logging.error("Unable to read file format '{}'. Exiting.\n")
        sys.exit(3)

    logging.info("done.\n")

    logging.debug("Read {} classified sequences across {} unique taxa.\n".format(sum(taxa_map.values()),
                                                                                 len(taxa_map.keys())))

    return taxa_map


def get_num_fastx_records(fastx_file: str) -> int:
    logging.info("Calculating the number of sequence records in '{}'... ".format(fastx_file))
    fastx_type = fastx_format_check(fastx_file)
    if fastx_type == 'fasta':
        fx = Fasta(file_name=fastx_file, build_index=False)
    elif fastx_type == 'fastq':
        fx = Fastq(file_name=fastx_file, build_index=False)
    else:
        logging.error("Unknown fastx type: '{}'\n".format(fastx_type))
        sys.exit(3)

    i = 0
    for _ in fx:
        i += 1

    logging.info("done.\n")

    logging.debug("Parsed {} sequence records.\n".format(i))

    return i


def bin_classes(taxa_abund: dict, positives: set, total_queries: int) -> (int, int, int, int):
    """
    Sorts the classifications into the four classes: true/false positives and true/false negatives
    TP - taxa in positives
    FP - taxa not in false positives
    FN - taxa that were not classified
    TN - 0

    :param taxa_abund: A dictionary mapping an organism name to its number of classifications
    :param positives: A list of expected organisms
    :param total_queries: The total number of query sequences that should have been classified (as positives)
    :return:
    """
    tp = 0
    fp = 0
    for k, v in taxa_abund.items():  # type: (str, int)
        if k in positives:
            tp += v
        else:
            fp += v

    fn = total_queries - (tp + fp)
    tn = 0

    return tp, fp, tn, fn


def margin_calc(chars: int, margin_size=40):
    return " "*(margin_size - chars)


def classification_measures_summary(tp, fp, tn=0, fn=0) -> None:
    tpr = float(tp/(tp + fn))
    ppv = float(tp/(tp + fp))
    tpr_str = "True positive rate (TPR)"
    fnr_str = "False negative rate (FNR)"
    ppv_str = "Positive predictive value (PPV)"
    fdr_str = "False discovery rate (FDR)"
    f1_str = "F1-score"
    summary_str = "{}{}{:.3f}\n".format(tpr_str, margin_calc(len(tpr_str)), tpr) + \
                  "{}{}{:.3f}\n".format(fnr_str, margin_calc(len(fnr_str)), 1-tpr) + \
                  "{}{}{:.3f}\n".format(ppv_str, margin_calc(len(ppv_str)), ppv) + \
                  "{}{}{:.3f}\n".format(fdr_str, margin_calc(len(fdr_str)), 1-ppv) + \
                  "{}{}{:.3f}\n".format(f1_str, margin_calc(len(f1_str)), 2*(ppv*tpr)/(ppv+tpr))
    logging.info(summary_str)
    return


# def format_roc(y_true: np.array, y_counts: np.array) -> (np.array, np.array, np.array):
#     """
#     :type y_true : array, shape = [n_samples]
#     :param y_true: True binary labels. Labels should be either either {-1, 1} or {0, 1}
#     :param y_counts:
#     :return:
#     """
#     # Compute ROC curve and ROC area for each class
#     fpr, tpr, _ = roc_curve(y_test, y_counts)
#     roc_auc = auc(fpr, tpr)
#
#     return tpr, fpr, roc_auc


def plot_roc(tpr: np.array, fpr: np.array, roc_auc: np.array, fig_name: str) -> None:
    plt.figure()
    lw = 2
    plt.plot(fpr, tpr, color='darkorange',
             lw=lw, label='ROC curve (area = %0.2f)' % roc_auc[2])
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic example')
    plt.legend(loc="lower right")
    plt.show()
    return


def metag_roc(sys_args):
    args = get_options(sys_args)

    prep_logging("./roc_log.txt", verbosity=args.verbose)
    # Read the inputs
    num_queries = get_num_fastx_records(args.fastx)
    positive_taxa = read_lines_to_set(args.p_names)
    taxa_map = get_classifications(args.c_tbl, "bracken")

    # Perform ROC analysis and plot
    tp, tn, fp, fn = bin_classes(taxa_map, positive_taxa, num_queries)
    classification_measures_summary(tp=tp, fp=fp, fn=fn)

    return


if __name__ == "__main__":
    metag_roc(sys.argv[1:])
