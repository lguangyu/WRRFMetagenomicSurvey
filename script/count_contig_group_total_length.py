#!/usr/bin/env python3

import Bio.SeqIO
import argparse
import collections
import io
import json
import sys


def get_args():
	ap = argparse.ArgumentParser(description="count total contig length in a "
		"fasta file in each group provided by a group label file")
	ap.add_argument("--fasta", "-f", type=str, required=True,
		metavar="fasta",
		help="input fasta file (required)")
	ap.add_argument("--grouping", "-g", type=str, required=True,
		metavar="txt",
		help="input grouping file, the accepted format is a 2-column table, "
			"with the first column listing sequence names, and the second "
			"column being the label/name of group the corresponding sequence "
			"belongs to (required)")
	ap.add_argument("--delimiter", "-d", type=str, default="\t",
		metavar="char",
		help="the delimiter character used in --grouping file and also affects "
			"the output if --output-format=txt [<tab>]")
	ap.add_argument("--output-format", type=str, default="txt",
		choices=["txt", "json"],
		help="output file format [txt]")
	ap.add_argument("--output", "-o", type=str, default="-",
		metavar="file",
		help="output file [stdout]")

	# parse and refine args
	args = ap.parse_args()
	if args.output == "-":
		args.output = sys.stdout
	return args


def get_fp(f, *ka, factory=open, **kw):
	if isinstance(f, io.IOBase):
		ret = f
	elif isinstance(f, str):
		ret = factory(f, *ka, **kw)
	else:
		raise TypeError("the first argument of get_fp must be str or io.IOBase,"
			" got '%s'" % type(f).__name__)
	return ret


def load_fasta_sequence_length(f, format="fasta") -> dict:
	ret = dict()
	for seq in Bio.SeqIO.parse(f, format="fasta"):
		ret[seq.id] = len(seq)
	return ret


def count_seq_group_total_length(group_file: str, seq_length: dict, *,
		delimiter="\t") -> collections.Counter:
	ret = collections.Counter()
	with get_fp(group_file, "r") as fp:
		for line in fp.read().splitlines():
			fields = line.split(delimiter)
			if len(fields) != 2:
				raise ValueError("wrong number of fields encountered in line: "
					"'%s', expect 2, got %u" % (line, len(fields)))
			seq_id, group = fields
			ret[group] += seq_length.get(seq_id, 0)
	return ret


def save(f, group_length, *, format="txt", delimiter="\t"):
	with get_fp(f, "w") as fp:
		if format == "txt":
			# sort in descending order
			for group, length in group_length.most_common():
				print(group + delimiter + str(length), file=fp)
		elif format == "json":
			json.dump(fp, group_length)
		else:
			raise ValueError("unaccepted output format '%s'" % format)
	return


def main():
	args = get_args()
	seq_length = load_fasta_sequence_length(args.fasta)
	group_length = count_seq_group_total_length(args.grouping, seq_length,
		delimiter=args.delimiter)
	save(args.output, group_length, format=args.output_format,
		delimiter=args.delimiter)
	return


if __name__ == "__main__":
	main()
