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
		help="input renamed BED-like file, with 1st column as CDS id, 2nd as "
			"start position, 3rd as stop position, 5th as number of reads hit ")
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


def calc_cds_tpm(ifile, ofile, total_read_count: PosInt):
	rpk = dict()
	total_rpk = 0
	with get_fp(ifile, "r") as ifp:
		for line in ifp:
			cid, start, stop, depth, n_hit = line.rstrip().split("\t")
			length = (int(stop) - int(start) + 1)
			assert length > 0
			v = int(n_hit) / length * 1e3
			total_rpk += v
			rpk[cid] = v
	with get_fp(ofile, "w") as ofp:
		for cid in rpk:
			print("%s\t%.4f" % (cid, rpk[cid] / total_rpk * 1e6), file=ofp)
	return


def main():
	args = get_args()
	calc_cds_tpm(args.input, args.output, args.total_read_count)
	return


if __name__ == "__main__":
	main()
