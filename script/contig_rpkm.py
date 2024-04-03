#!/usr/bin/env python3

import argparse
import io
import sys


class PosInt(int):
	def __new__(cls, *ka, **kw):
		self = super().__new__(cls, *ka, **kw)
		if self <= 0:
			raise ValueError("%s must be positive, got %d"
				% (cls.__name__, self))
		return self


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("input", type=str, nargs="?", default="-",
		help="input contig read mapping stats from samtools idxstats ")
	ap.add_argument("-c", "--total-read-count", type=PosInt, required=True,
		metavar="int",
		help="total read count mapped to contigs (required)")
	ap.add_argument("-o", "--output", type=str, default="-",
		metavar="bed",
		help="output rpkm file [stdout]")
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


def calc_contig_rpkm(ifile, ofile, total_read_count: PosInt):
	with get_fp(ifile, "r") as ifp:
		with get_fp(ofile, "w") as ofp:
			for line in ifp:
				cid, length, n_hit, _ = line.rstrip().split("\t")
				length = int(length)
				if length == 0:
					print(line, end="", file=sys.stderr)
					length = 1
				rpkm = int(n_hit) / length / total_read_count * 1e9
				print("%s\t%.4f" % (cid, rpkm), file=ofp)
	return


def main():
	args = get_args()
	calc_contig_rpkm(args.input, args.output, args.total_read_count)
	return


if __name__ == "__main__":
	main()
