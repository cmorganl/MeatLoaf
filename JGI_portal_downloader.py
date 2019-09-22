#!/usr/bin/env python3

import os
import sys
import re
import io
import argparse
import pycurl
import logging
import xml.etree.ElementTree as ET

__author__ = 'Connor Morgan-Lang'

c = pycurl.Curl()
# How to install pycurl:
# python3 -m pip uninstall pycurl
# export PYCURL_SSL_LIBRARY=nss
# python3 -m pip install --compile --install-option="--with-nss" --no-cache-dir pycurl

website_url = "https://genome.jgi.doe.gov/"
portal_path = "portal/ext-api/downloads/get-directory?organism="


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


def prep_logging(log_file_name=None, verbosity=False):
    """
    Allows for multiple file handlers to be added to the root logger, but only a single stream handler.
    The new file handlers must be removed outside of this function explicitly
    :param log_file_name:
    :param verbosity:
    :return:
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
            logging.error("Unable to make directory '" + output_dir + "'.\n")
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


def progress(total_to_download, total_downloaded, total_to_upload=None, total_uploaded=None):
    try:
        percent_completed = (float(total_downloaded)*100) / total_to_download  # You are calculating amount uploaded
        rate = round(percent_completed, ndigits=2)
        completed_hash = "#" * int(rate)  # Calculate completed percentage
        spaces = " " * (100 - int(percent_completed))  # Calculate remaining completed rate
        buffer = '[%s%s] %s%%\n' % (completed_hash, spaces, rate)  # the pretty progress [####     ] 34%
    except ZeroDivisionError:
        buffer = ""
    if buffer:
        sys.stdout.write(buffer)
        sys.stdout.flush()
    return


def curl_debug(debug_type, msg):
    sys.stderr.write("debug: %s %s\n" % (repr(debug_type), repr(msg)))
    sys.stderr.flush()


def load_curl_opts():
    c.setopt(pycurl.VERBOSE, False)
    c.setopt(c.FOLLOWLOCATION, True)
    # c.setopt(c.NOPROGRESS, False)
    c.setopt(pycurl.DEBUGFUNCTION, curl_debug)
    try:
        c.setopt(c.XFERINFOFUNCTION, progress)
    except AttributeError:
        c.setopt(c.PROGRESSFUNCTION, progress)
    return


def get_options():
    parser = argparse.ArgumentParser(description="Utility script for programmatically downloading all project files"
                                                 " from the JGI Portal. Details on how this works can be found at"
                                                 " https://genome.jgi.doe.gov/portal/help/download.jsf#/api",
                                     add_help=False)
    reqd = parser.add_argument_group("Required arguments")
    user_cred = parser.add_argument_group("User credentials - only if you don't have cookies")
    opts = parser.add_argument_group("Optional arguments")

    reqd.add_argument("-p", "--portal_name", required=True, dest="name",
                      help="Short portal name in the URL from the JGI's 'Downloads' page for this project")
    reqd.add_argument("-o", "--output_dir", required=True, dest="output",
                      help="Path to the directory to write all downloaded outputs")

    opts.add_argument("-c", "--cookies", required=False, default=None,
                      help="Path to the cookies file created by curl")
    opts.add_argument("-x", "--extensions", required=False, default=None,
                      help="Only download files with these comma-separated extensions."
                           " Example: pdf,fastq,fasta  [ DEFAULT = download all the things ]")
    opts.add_argument("-d", "--dirs", required=False, default=None, dest="target_dirs",
                      help="Only download files from these comma-separated directories."
                           " Example: QC_Filtered_Raw_Data,Raw_Data  [ DEFAULT = download all the things ]")
    opts.add_argument("-h", "--help", action="help",
                      help="show this help message and exit")

    user_cred.add_argument("--username", dest="usr", required=False, default=None,
                           help="Username (or email) linked to JGI single sign-on. Don't provide if you have cookies.")
    user_cred.add_argument("--password", dest="pwd", required=False, default=None,
                           help="Password for your JGI single sign-on account. Don't provide if you have cookies.")

    args = parser.parse_args()

    if not os.path.isdir(args.output):
        try:
            os.mkdir(args.output)
        except OSError:
            logging.error("Unable to create new directory for output '" + args.output + "'.\n")
            sys.exit(3)

    if args.output[-1] != os.sep:
        args.output = args.output + os.sep

    if args.extensions:
        args.extensions = args.extensions.split(',')
        
    if args.target_dirs:
        args.target_dirs = args.target_dirs.split(',')

    if not args.cookies:
        logging.error("Function is not yet supported. Must provide a cookie file for now - sorry!\n")
        sys.exit(11)

    return args


def bake_cookies(pwd: str, usr: str) -> str:
    cookies_file = os.getcwd() + os.sep + "cookies"
    logging.info("Creating cookies file '" + cookies_file + "'... ")

    # Basically equivalent to:
    # curl 'https://signon-old.jgi.doe.gov/signon/create' \
    # --data-urlencode 'login=USER_NAME' --data-urlencode 'password=USER_PASSWORD' -c cookies > /dev/null

    # buffer = io.BytesIO()

    c.setopt(c.URL, "https://signon-old.jgi.doe.gov/signon/create")
    c.setopt(c.USERPWD, usr + ':' + pwd)
    c.setopt(c.COOKIEJAR, cookies_file)
    # c.setopt(c.WRITEFUNCTION, buffer.write)
    c.perform()

    # resp = buffer.getvalue()

    logging.info("done.\n")

    return cookies_file


def get_file_deets_curl():
    c.setopt(c.HEADER, 0)
    c.setopt(c.NOBODY, 1)
    c.perform()
    print(c.getinfo(c.CONTENT_LENGTH_DOWNLOAD))
    c.setopt(c.NOBODY, 0)
    print(c.getinfo(c.CONTENT_LENGTH_DOWNLOAD))
    sys.exit()
    return c.getinfo(c.CONTENT_LENGTH_DOWNLOAD)


def prep_and_curl(url_address, output_file, cookie_file):
    try:
        output_handler = open(output_file, 'wb')
    except IOError:
        logging.error("Unable to open file '" + output_file + "' for writing!\n")
        sys.exit(3)

    c.setopt(c.URL, url_address)
    c.setopt(c.COOKIEFILE, cookie_file)

    # filesize = get_file_deets_curl()
    filesize = 100
    if filesize < 0:
        logging.error("Size of " + output_file + " less than 0b. Something went wrong when retrieving it.\n")
        sys.exit(7)
    elif filesize == 0:
        logging.debug("Size of " + output_file + " is 0b. Creating it and moving on.\n")
        open(output_file, 'w').close()
    else:
        c.setopt(c.WRITEDATA, output_handler)
        # c.setopt(c.INFILESIZE, filesize)

        try:
            c.perform()
        except pycurl.error as error:
            errno, errstr = error
            logging.error("An error occurred while running curl: " + errstr + "\n")
            pass

        output_handler.close()
    return


def dl_file_list(portal_name: str, cookies: str, file_xml: str) -> None:
    logging.info("Downloading XML-file listing all target files... ")

    if not os.path.isfile(cookies):
        logging.error("Cookies file '" + cookies + "' doesn't exist!\n")
        sys.exit(3)

    page_address = website_url + portal_path + portal_name
    prep_and_curl(page_address, file_xml, cookies)

    logging.info("done.\n")

    return


def parse_files_xml(file_xml: str, ext=None) -> dict:
    dirs_n_files = dict()
    try:
        xml_handler = open(file_xml, 'r')
    except IOError:
        logging.error("Unable to open xml-file '" + file_xml + "' for reading!\n")
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

        logging.info("Downloading files in " + dir_name + "\n")
        dir_path = output + dir_name + os.sep
        if not os.path.isdir(dir_path):
            os.mkdir(dir_path)
        for file_name in sorted(dirs_n_files[dir_name]):
            if ext and file_name.split('.')[-1] not in ext:
                continue
            file_path = dir_path + file_name
            # There are some bad characters so we're going to clean them up here
            file_address = website_url + re.sub(r"&amp;", '&', dirs_n_files[dir_name][file_name])
            logging.debug("Downloading " + file_name + "... ")
            prep_and_curl(file_address, file_path, cookies)
            logging.debug("done.\n")
            num_dl += 1
        logging.info("Finished downloading all file in " + dir_name + ".\n")

    if num_dl == 0:
        logging.error("No files were downloaded.\n")
        sys.exit(5)

    return


def main():
    args = get_options()
    log_file = args.output + "JGI_downloader_log.txt"
    prep_logging(log_file)
    load_curl_opts()

    if not args.cookies:
        args.cookies = bake_cookies(args.pwd, args.usr)

    file_list = args.output + args.name + "_files.xml"
    dl_file_list(args.name, args.cookies, file_list)
    addresses = parse_files_xml(file_list, args.extensions)
    dl_files(addresses, args.cookies, args.output, args.extensions, args.target_dirs)
    return


main()
