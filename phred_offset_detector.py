#!/usr/bin/env python
# -*- mode: python -*-

import sys
import argparse

# To determine the Phred score of fastq files from their quality line

#1) Save the fourth line of each fastq sequence in a dictionary
#	-open the fastq file with a raw_input call saved as a string, then create a file object (read-only)
#	-extract the fourth line of each fastq sequence by creating a dictionary and each entry would be a new key:value pair. While loop where n+=1 each loop iteration.
#	- make each quality line a list of characters
#	- close the fastq file
#2) Save the Phred score in the dict by finding if the quality line contains any characters that are equal to less than 64 
#3) Return the offset of the Phred score at the end with an if statement (if the value of any key is lower than 64, return 33)

#parse through each list in the dictionary and identify the integer corresponding to each value and identify if there are values lower than 64, if there are then the Phred offset is 33, else it is 64.

# The Key is the decimal value and the Value is the respective glyph

Ascii_dict = {
" ":32,
"!":33,
"\"":34,
"#":35,
"$":36,
"%":37,
"&":38,
"'":39,
"(":40,
")":41,
"*":42,
"+":43,
",":44,
"-":45,
".":46,
"/":47,
"0":48,
"1":49,
"2":50,
"3":51,
"4":52,
"5":53,
"6":54,
"7":55,
"8":56,
"9":57,
":":58,
";":59,
"<":60,
"=":61,
">":62,
"?":63,
"@":64,
"A":65,
"B":66,
"C":67,
"D":68,
"E":69,
"F":70,
"G":71,
"H":72,
"I":73,
"J":74,
"K":75,
"L":76,
"M":77,
"N":78,
"O":79,
"P":80,
"Q":81,
"R":82,
"S":83,
"T":84,
"U":85,
"V":86,
"W":87,
"X":88,
"Y":89,
"Z":90,
"[":91,
"\\":92,
"]":93,
"^":94,
"_":95,
"`":96,
"a":97,
"b":98,
"c":99,
"d":100, 
"e":101,
"f":102,
"g":103,
"h":104,
"i":105,
"j":106,
"k":107,
"l":108,
"m":109,
"n":110,
"o":111,
"p":112,
"q":113,
"r":114,
"s":115,
"t":116,
"u":117,
"v":118,
"w":119,
"x":120,
"y":121,
"z":122,
"{":123,
"|":124,
"}":125,
"~":126,
}

def get_options():
	"""
	Use argparse to get command-line options
	:return: object with file name
	"""
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--file",
						help="FASTQ file to detect Phred offset score.",
						type=str,
						required=True)
	parser.add_argument("-n", "--num_sequences",
						help="The number of sequences to estimate the Phred offset score on [DEFAULT = 10000].",
						default=10000,
						type=int)
	args = parser.parse_args()
	return args

def determine_offset(Scores):
	for x in Scores:
		if Ascii_dict[x] < 59:
			o="-phred33"
			break
		elif Ascii_dict[x] > 74:
			o="-phred64"
			break
	return o

def get_qual_scores(fastq, numLines):
	qual_lines = ""
	with open(fastq, 'r') as input:
		while numLines >= 0:
			line = input.readline().strip()
			if len(line) == 0:
				print "File is fewer than 40000 lines. Basing phred offset estimate on",\
					  (40003 - numLines)/4, "sequences."
				break
			if numLines % 4 == 0:
				qual_lines += line
			else:
				pass
			numLines -= 1
	return qual_lines

def main():
	args = get_options()
	num_lines = (args.num_sequences*4) + 3
	qual_lines = get_qual_scores(args.file, num_lines)
	offset = determine_offset(qual_lines)
	print offset

main()