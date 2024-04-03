#!/usr/bin/env python3

import argparse
import io
import sys


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("input", type=str, nargs="?", default="-",
		help="input BED-like file, with at least the first 3 columns has the "
			"same definition as a regular BED file")
	ap.add_argument("-o", "--output", type=str, default="-",
		metavar="bed",
		help="output BED-like file [stdout]")
	# parse and refine args
	args = ap.parse_args()
	if args.input == "-":
		args.input = sys.stdin
	if args.output == "-":
		args.output = sys.stdout
	return args


def get_fp(f, *ka, factory=open, **kw):
	if isinstance(f, io.IOBase):
		ret = f
	elif isinstance(f, str):
		ret = factory(f, *ka, **kw)
	else:
		raise TypeError("first arg of get_fp() must be io.IOBase or str, "
			"got '%s'" % type(f).__name__)
	return ret


def rename_bed_contig(ifile, ofile):
	curr_prefix = None
	curr_id = 0
	with get_fp(ifile, "r") as ifp:
		with get_fp(ofile, "w") as ofp:
			for line in ifp:
				cid, data = line.rstrip().split("\t", maxsplit=1)
				if cid != curr_prefix:
					curr_prefix = cid
					curr_id = 0
				curr_id += 1
				print(("\t").join(["%s_%u" % (curr_prefix, curr_id), data]),
					file=ofp)
	return


def main():
	args = get_args()
	rename_bed_contig(args.input, args.output)
	return


if __name__ == "__main__":
	main()
