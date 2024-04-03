#!/usr/bin/env python3

import argparse
import collections
import io
import numpy
import sys


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("-c", "--cds-rpkm", type=str, required=True,
		metavar="tsv",
		help="input per CDS rpkm (required)")
	ap.add_argument("-k", "--cds-ko-map", type=str, required=True,
		metavar="tsv",
		help="input CDS to KO map (required)")
	ap.add_argument("-o", "--output", type=str, default="-",
		metavar="tsv",
		help="output per ko rpkm file [stdout]")
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
		raise TypeError("first arg of get_fp() must be io.IOBase or str, "
			"got '%s'" % type(f).__name__)
	return ret


def load_cds_ko_map(fname) -> dict:
	ret = dict()
	with get_fp(fname, "r") as fp:
		for line in fp:
			cid, *ko = line.rstrip().split("\t")
			if not ko:
				continue
			ret[cid] = ko[0]
	return ret


def resolve_ko_rpkm(cds_rpkm_fname, cds_ko_map) -> dict:
	ret = collections.defaultdict(list)
	with get_fp(cds_rpkm_fname, "r") as fp:
		for line in fp:
			cid, rpkm = line.rstrip().split("\t")
			if cid not in cds_ko_map:
				continue
			ret[cds_ko_map[cid]].append(float(rpkm))
	return ret


def save_ko_rpkm(ko_rpkm, fname):
	with get_fp(fname, "w") as fp:
		fields = ["#ko", "median", "mean", "sum"]
		print(("\t").join(fields), file=fp)
		for k in sorted(ko_rpkm.keys()):
			v = ko_rpkm[k]
			fields = [
				k,
				"%.4f" % numpy.median(v),
				"%.4f" % numpy.mean(v),
				"%.4f" % numpy.sum(v),
			]
			print(("\t").join(fields), file=fp)
	return


def main():
	args = get_args()
	cds_ko_map = load_cds_ko_map(args.cds_ko_map)
	ko_rpkm = resolve_ko_rpkm(args.cds_rpkm, cds_ko_map)
	save_ko_rpkm(ko_rpkm, args.output)
	return


if __name__ == "__main__":
	main()
