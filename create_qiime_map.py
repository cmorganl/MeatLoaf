#!/usr/bin/env python

try:
    import argparse
    import sys
    from glob import glob
    import os
    import traceback
    import re
    import subprocess
    import shutil
    from itertools import product
    from random import random
    from time import strftime
    from time import sleep
    from math import ceil
except ImportWarning:
    sys.stderr.write("Could not load some user defined module functions")
    sys.stderr.write(traceback.print_exc(10))
    sys.exit(3)

__author__ = 'Connor Morgan-Lang'


# Heng Li's function for reading a FASTQ file using Generators
def readfq(fp):  # this is a generator function
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last: break
        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs);  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break


def stdprint(obj, channel, cap=""):
    str_obj = str(obj)
    if type(str_obj) is not str:
        sys.exit("TypeError: stdprint only accepts string objects!")
    elif str(channel) == "err":
        sys.stderr.write(str_obj + cap)
        sys.stderr.flush()
    elif str(channel) == "out":
        sys.stdout.write(str_obj + cap)
        sys.stdout.flush()
    else:
        sys.stderr.write("ERROR: Unrecognized input to stdprint:\n")
        sys.stderr.write(str_obj + "\n" + channel + "\n")
        sys.exit()


def launch_write_command(cmd_list):
    proc = subprocess.Popen(' '.join(cmd_list),
                            shell=True,
                            preexec_fn=os.setsid,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    stdout = proc.communicate()[0].decode("utf-8")
    if stdout is None:
        stdout = ""
    proc.wait()
    return stdout, proc.returncode


def get_options():
    parser = argparse.ArgumentParser(description="A script for creating a QIIME-compatible metadata file"
                                                 " from a FASTQ directory and a csv file, mapping the "
                                                 "FASTQ files to other metadata. "
                                                 "ASSUMES unique barcode sequences NEED to be created and prepended.\n")
    parser.add_argument("-m", "--metadata", type=str, required=True,
                        help="The environmental data for all FASTQs. "
                             "This file must contain a header line!")
    parser.add_argument("-i", "--itags_dir", type=str, required=True,
                        help="Path to  directory containing iTags in FASTQ format.")
    parser.add_argument("-d", "--output_dir", type=str, required=False,
                        default="QIIME_ready_iTags",
                        help="The directory for the forward and reverse concatenated FASTQ files."
                             "[ DEFAULT = QIIME_ready_iTags ]")
    parser.add_argument("-o", "--output_map", type=str, required=False,
                        help="Name of the final QIIME-compatible sample mapping file.",
                        default="QIIME_fastq_sample_map.tsv")
    parser.add_argument("-s", "--sep", type=str, required=False,
                        default="\t",
                        help="Field separator used for parsing metadata file. [DEFAULT = '\\t' ]")
    parser.add_argument("--chunk", default=False, action="store_true",
                        required=False,
                        help="Output the FASTQs as concatenated chunks, containing 1Gb of sequence per file."
                             "Default is to output a separate, modified FASTQ file for each sample.")

    args = parser.parse_args()

    if args.itags_dir[-1] != os.sep:
        args.itags_dir += os.sep

    if args.output_dir[-1] != os.sep:
        args.output_dir += os.sep

    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)
    else:
        sys.stdout.write("WARNING: output directory '" + args.output_dir + "' already exists!\n"
                         "You have 5 seconds to kill before these files are overwritten... ")
        sys.stdout.flush()
        sleep(5)
        sys.stdout.write("BOOM.\n")
        sleep(1)
        shutil.rmtree(args.output_dir)
        os.makedirs(args.output_dir)

    return args


def get_fastq_files(data_dir):
    fastq_files = list()
    accepted_extensions = [".fq", ".fastq", ".fastq.gz", ".fq.gz"]
    if not os.path.isdir(data_dir):
        sys.stderr.write("ERROR: --itags_dir directory provided does not exist!\n")
        sys.exit(3)
    else:
        dir_files = glob(data_dir + "*")
        # Remove non-fastq files from the list of files in itags_dir
        for dir_file in dir_files:
            for ext in accepted_extensions:
                if re.search("^.*" + re.escape(ext) + "$", dir_file):
                    fastq_files.append(dir_file)
                    break

    sys.stdout.write(str(len(fastq_files)) + " fastq files found in iTags directory.\n")

    if len(fastq_files) == 0:
        sys.exit(3)

    return fastq_files


def read_metadata_file(metadata_file, sep="\t"):
    num_lines = 0
    sample_dat = dict()
    try:
        metadata_handler = open(metadata_file, 'r')
    except IOError:
        sys.stderr.write("ERROR: Unable to open metadata " + metadata_file + "file for reading!\n")
        sys.exit(5)

    header = metadata_handler.readline().strip()
    header_fields = header.split(sep)
    # TODO: Parse header and ensure the few required fields are present
    if header_fields[0] != "SampleID":
        if re.search("sample", header_fields[0], re.IGNORECASE):
            sys.stderr.write("WARNING: first column in " + metadata_file + " is not 'SampleID' but close enough.\n")
            header_fields[0] = "SampleID"
        else:
            sys.stderr.write("ERROR: first column in " + metadata_file + " does not resemble 'SampleID'.\n")
            sys.stderr.write("The first column must have the name of the sample!\n")
            sys.exit(5)

    line = metadata_handler.readline()
    while line:
        num_lines += 1
        x = 1
        line = line.strip()
        fields = line.split(sep)
        if len(fields) != len(header_fields):
            sys.stderr.write("ERROR: line " + str(num_lines + 1) +
                             " of metadata file does not contain the same number of fields as the header!\n")
            sys.exit(5)
        sample = fields[0]
        if sample in sample_dat:
            sys.stderr.write("ERROR: sample " + sample + " exists in metadata file more than once!\n")
            sys.exit(5)
        else:
            sample_dat[sample] = dict()
        while x < len(header_fields):
            sample_dat[sample][header_fields[x]] = fields[x]
            x += 1

        line = metadata_handler.readline()

    sys.stdout.write(str(num_lines) + " lines read from metadata file " + metadata_file + "\n")
    return sample_dat


def random_subset(iterator, k):
    result = []
    n = 0

    for item in iterator:
        n += 1
        if len(result) < k:
            result.append(item)
        else:
            s = int(random() * n)
            if s < k:
                result[s] = item

    return result


def get_unique_barcodes(num_samples, length=11):
    unique_barcodes = random_subset([''.join(i) for i in product("ACTG", repeat=length)], len(num_samples))
    return unique_barcodes


def assign_barcodes(sample_metadata, barcodes):
    barcode_acc = 0
    for sample in sample_metadata:
        sample_metadata[sample]["BarcodeSequence"] = barcodes[barcode_acc]
        barcode_acc += 1
    return sample_metadata


def match_sample_to_filename(fastq_sample, sample_names, log=None):
    guessed_samples = list()
    # Try and match the real sample name to the one parsed from the FASTQ file:
    for parsed_sample in sample_names:
        if parsed_sample == fastq_sample:
            return parsed_sample
        elif re.search(re.escape(parsed_sample), fastq_sample):
            guessed_samples.append(parsed_sample)
        elif re.search(re.escape(fastq_sample), parsed_sample):
            guessed_samples.append(parsed_sample)
        else:
            pass
    if len(guessed_samples) == 1:
        return guessed_samples[0]
    elif len(guessed_samples) > 1:
        # Find the longest common super-string
        lcs = ""
        for sample in guessed_samples:
            if len(sample) > len(lcs):
                lcs = sample
        if log:
            log.write("WARNING: pairing sample ID '" + lcs + "' with FASTQ prefix '" + fastq_sample + "'\n")
        return lcs
    return ""


def get_fastq_pair(fastq):
    if os.path.isfile(re.sub("_R1", "_R2", fastq)):
        return re.sub("_R1", "_R2", fastq)
    elif os.path.isfile(re.sub("1.fastq", "2.fastq", fastq)):
        return re.sub("1.fastq", "2.fastq", fastq)
    elif get_fastq_pair(fastq + ".gz"):
        return get_fastq_pair(fastq + ".gz")
    else:
        sys.stderr.write("Unable to find paired reverse FASTQ for " + fastq + "\n")
        sys.exit(11)


def read_and_format_fastqs(args, fastq_list, sample_dat, linker_primer="GCTAGCAA"):

    try:
        log_handler = open("QIIME_mapper_log.txt", 'w')
    except IOError:
        sys.stderr.write("ERROR: unable to open log file for writing!\n")
        sys.exit(11)

    sys.stdout.write("Reading and formatting the FASTQ files for QIIME analysis... ")
    sys.stdout.flush()
    chunk = 1
    output_strings = []
    formatted_fastq_file = ""
    for fastq in fastq_list:
        # Continue if this is a reverse read
        if re.search("_R2", fastq) or re.search("_2.", fastq):
            continue

        # Decompress if this is gzipped
        reverse_fq = get_fastq_pair(fastq)
        if re.match(".*.gz$", fastq):
            stdout, retcode = launch_write_command(["gunzip", fastq])
            if retcode != 0:
                sys.stderr.write("ERROR: `gunzip " + fastq + "` did not complete successfully!\n")
                log_handler.write("ERROR: `gunzip " + fastq + "` did not complete successfully!\n")
                sys.exit()
            fastq = fastq.rstrip(".gz")

        sample_guess = '.'.join(os.path.basename(fastq).split('.')[:-1])
        if re.search("^.*_R1", sample_guess):
            sample_guess = sample_guess.split("_R1")[0]
        # Get the matching sample name
        sample_match = match_sample_to_filename(sample_guess, sample_dat.keys(), log_handler)
        if not sample_match:
            log_handler.write("ERROR: Unable to match FASTQ sample ID (" + sample_guess + ") to one in metadata!\n")
            continue
        # Select a random barcode and save the mapping information
        tag = sample_dat[sample_match]["BarcodeSequence"] + linker_primer

        # Open the file to be read
        try:
            fastq_handler = open(fastq, 'r')
        except IOError:
            sys.stderr.write("ERROR: unable to open FASTQ file (" + fastq + ") for reading!\n")
            log_handler.write("ERROR: unable to open FASTQ file (" + fastq + ") for reading!\n")
            sys.exit(7)

        # Decide whether to clear the output strings (if we're writing sample-wise) or not
        if not args.chunk:
            formatted_fastq_file = args.output_dir + sample_match + "R1.fastq"
            output_strings = []
        else:
            formatted_fastq_file = args.output_dir + "chunk_" + str(chunk) + "R1.fastq"

        # Read and format the FASTQ file
        for name, seq, qual in readfq(fastq_handler):
            formatted_line = "\n".join(['@' + name, tag + seq, '+', len(tag)*'5' + qual])
            output_strings.append(formatted_line)
        fastq_handler.close()

        if args.chunk and len(output_strings) >= 1000000:
            write_string_to_file(formatted_fastq_file, "\n".join(output_strings))
            chunk += 1
            output_strings = []
        elif not args.chunk:
            write_string_to_file(formatted_fastq_file, "\n".join(output_strings))

        # Write the reverse FASTQ file to its appropriate destination
        if re.match(".*.gz$", reverse_fq):
            fq_copy_command = ["gunzip", "-c", reverse_fq,
                               ">>" + re.sub("R1.fastq", "R2.fastq", formatted_fastq_file)]
        else:
            fq_copy_command = ["cat", reverse_fq,
                               ">>" + re.sub("R1.fastq", "R2.fastq", formatted_fastq_file)]
        stdout, retcode = launch_write_command(fq_copy_command)
        if retcode != 0:
            sys.stderr.write("ERROR: `" + ' '.join(fq_copy_command) + "` did not complete successfully!\n")
            log_handler.write("ERROR: `" + ' '.join(fq_copy_command) + "` did not complete successfully!\n")
            sys.exit()

    if args.chunk and formatted_fastq_file:
        write_string_to_file(formatted_fastq_file, "\n".join(output_strings))

    sys.stdout.write("done.\n")
    log_handler.close()

    return


def write_string_to_file(filename, output_string):
    try:
        file_handler = open(filename, 'w')
    except IOError:
        sys.stderr.write("ERROR: Unable to open file (" + filename + ") for writing!\n")
        sys.exit(9)

    file_handler.write(output_string + "\n")

    file_handler.close()
    return


def create_qiime_map(args, sample_metadata, linker_primer="GCTAGCAA"):
    qiime_map_strings = list()
    header = list()
    header_value_order = {0: "SampleID",
                          1: "BarcodeSequence",
                          2: "LinkerPrimerSequence"}
    for sample in sample_metadata:
        # Set up the header order with the first sample:
        if not header:
            i = 3
            for field in list(sample_metadata[sample].keys()):
                if field not in header_value_order.values():
                    header_value_order[i] = field
                    i += 1
            for order in sorted(header_value_order, key=int):
                header.append(header_value_order[order])
            # Save the header to the final QIIME output string
            qiime_map_strings.append('#' + "\t".join(header) + "\tDescription")
        barcode = sample_metadata[sample]["BarcodeSequence"]
        sample_metadata[sample]["LinkerPrimerSequence"] = linker_primer
        sample_line = [sample, barcode, linker_primer]
        k = 3
        while k <= len(sample_metadata[sample].keys()):
            sample_line.append(sample_metadata[sample][header_value_order[k]])
            k += 1
        qiime_map_strings.append("\t".join(sample_line) + "\t")

    write_string_to_file(args.output_map, "\n".join(qiime_map_strings))

    return


def main():
    args = get_options()
    # Step 1: Find all the samples from the provided FASTQ directory
    fastq_files = get_fastq_files(args.itags_dir)
    # Step 2: Read metadata file and pull out the descriptors for each sample
    sample_dat = read_metadata_file(args.metadata, args.sep)
    # Step 3: Create unique barcode sequences for each sample
    if len(fastq_files) == len(sample_dat) or len(fastq_files)/2 == len(sample_dat):
        barcodes = get_unique_barcodes(sample_dat, 6)
        sample_dat = assign_barcodes(sample_dat, barcodes)
    else:
        sys.stderr.write("ERROR: Number of samples in the metadata file and the iTags directory are not equal!\n")
        sys.exit(-1)
    # Step 4: Write final FASTQ and QIIME mapping files
    create_qiime_map(args, sample_dat)
    # Step 5: Prepend standard LinkerPrimerSequences and barcodes to each FASTQ sequence with Phred 40 scores
    read_and_format_fastqs(args, fastq_files, sample_dat)


main()
