#!/usr/bin/env python3

import os
import sys
import re
import argparse
import pycurl
import xml.etree.ElementTree as ET

__author__ = 'Connor Morgan-Lang'

c = pycurl.Curl()
# How to install pycurl:
# python3 -m pip uninstall pycurl
# export PYCURL_SSL_LIBRARY=nss
# python3 -m pip install --compile --install-option="--with-nss" --no-cache-dir pycurl

jgi_web = "https://genome.jgi.doe.gov/"


def get_options():
    parser = argparse.ArgumentParser(description="Utility script for programmatically downloading all project files"
                                                 " from the JGI Portal. Details on how this works can be found at"
                                                 " https://genome.jgi.doe.gov/portal/help/download.jsf#/api")
    parser.add_argument("-p", "--portal_name", required=False, dest="name",
                        help="Short portal name in the URL from the JGI's 'Downloads' page for this project")
    parser.add_argument("-c", "--cookies", required=True,
                        help="Path to the cookies file created by curl")
    parser.add_argument("-o", "--output_dir", required=True, dest="output",
                        help="Path to the directory to write all downloaded outputs")
    parser.add_argument("-x", "--extensions", required=False, default=None,
                        help="Only download files with these comma-separated extensions."
                             " Example: pdf,fastq,fasta  [ DEFAULT = download all the things ]")
    parser.add_argument("-d", "--dirs", required=False, default=None, dest="target_dirs",
                        help="Only download files from these comma-separated directories."
                             " Example: QC_Filtered_Raw_Data,Raw_Data  [ DEFAULT = download all the things ]")
    args = parser.parse_args()

    if not os.path.isdir(args.output):
        try:
            os.mkdir(args.output)
        except OSError:
            sys.stderr.write("ERROR: unable to create new directory for output '" + args.output + "'.\n")
            sys.exit(3)

    if args.output[-1] != os.sep:
        args.output = args.output + os.sep

    if args.extensions:
        args.extensions = args.extensions.split(',')
        
    if args.target_dirs:
        args.target_dirs = args.target_dirs.split(',')

    return args


def prep_and_curl(url_address, output_file, cookie_file):
    try:
        output_handler = open(output_file, 'wb')
    except IOError:
        sys.stderr.write("ERROR: Unable to open xml file list " + output_file + " for writing!\n")
        sys.exit(3)

    c.setopt(c.URL, url_address)
    c.setopt(c.COOKIEFILE, cookie_file)
    c.setopt(c.WRITEDATA, output_handler)
    c.perform()

    output_handler.close()
    return


def dl_file_list(portal_name: str, cookies: str, file_xml: str) -> None:
    sys.stdout.write("Downloading XML with listing all files... ")
    sys.stdout.flush()

    if not os.path.isfile(cookies):
        sys.stderr.write("ERROR: cookies file '" + cookies + "' doesn't exist!\n")
        sys.exit(3)

    page_address = jgi_web + "portal/ext-api/downloads/get-directory?organism=" + portal_name
    prep_and_curl(page_address, file_xml, cookies)

    sys.stdout.write("done.\n")

    return


def parse_files_xml(file_xml: str, ext=None) -> dict:
    dirs_n_files = dict()
    try:
        xml_handler = open(file_xml, 'r')
    except IOError:
        sys.stderr.write("ERROR: unable to open xml-file '" + file_xml + "' for reading!\n")
        sys.exit(3)

    tree = ET.parse(xml_handler)
    root = tree.getroot()
    for child in root.iter():  # type: ET.Element
        if child.tag == "folder":
            k, v = list(child.attrib.items())[0]
            # replace all spaces with underscores
            dir_name = re.sub(r' ', '_', v)
            files_addresses = dict()
            # Load all the files for each directory
            for grand_kid in child.iter("file"):
                if grand_kid.tag == "file":
                    filename = ""
                    file_url = ""
                    for k, v in sorted(grand_kid.attrib.items()):
                        if k == "filename":
                            filename = v
                        elif k == "url":
                            file_url = v
                        else:
                            pass
                    # filter file names by extension
                    if ext:
                        file_ext = filename.split('.')[-1]
                        if file_ext in ext:
                            files_addresses[filename] = file_url
                    else:
                        files_addresses[filename] = file_url
            if files_addresses:
                dirs_n_files[dir_name] = files_addresses
    xml_handler.close()
    return dirs_n_files


def dl_files(dirs_n_files: dict, cookies, output: str, ext=None, target_dirs=None) -> None:
    num_dl = 0
    for dir_name in sorted(dirs_n_files):
        if target_dirs and dir_name not in target_dirs:
            continue

        sys.stdout.write("Downloading files in " + dir_name + "... ")
        sys.stdout.flush()
        dir_path = output + dir_name + os.sep
        if not os.path.isdir(dir_path):
            os.mkdir(dir_path)
        for file_name in sorted(dirs_n_files[dir_name]):
            if ext and file_name.split('.')[-1] not in ext:
                continue
            file_path = dir_path + file_name
            # There are some bad characters so we're going to clean them up here
            file_address = jgi_web + re.sub(r"&amp;", '&', dirs_n_files[dir_name][file_name])

            prep_and_curl(file_address, file_path, cookies)
            num_dl += 1
        sys.stdout.write("done.\n")

    if num_dl == 0:
        sys.stderr.write("ERROR: no files were downloaded.\n")
        sys.exit(5)

    return


def main():
    args = get_options()
    file_list = args.output + args.name + "_files.xml"
    dl_file_list(args.name, args.cookies, file_list)
    addresses = parse_files_xml(file_list, args.extensions)
    dl_files(addresses, args.cookies, args.output, args.extensions, args.target_dirs)
    return


main()
