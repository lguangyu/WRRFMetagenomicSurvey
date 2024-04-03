#!/usr/bin/env python3

import argparse
import collections
import io
import json
import numpy
import sys


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("-k", "--ko-rpkm", type=str, required=True,
		metavar="tsv",
		help="input per KO rpkm (required)")
	ap.add_argument("-m", "--ko-module-map", type=str, required=True,
		metavar="json",
		help="input info for CDS contained in module (required)")
	ap.add_argument("-o", "--output", type=str, default="-",
		metavar="tsv",
		help="output per module rpkm file [stdout]")
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


def load_ko_module_map(fname) -> dict:
	# module -> ko
	with get_fp(fname, "r") as fp:
		raw = json.load(fp)
	ko_module_map = collections.defaultdict(set)
	for m, ko_list in raw.items():
		for ko in ko_list:
			ko_module_map[ko].add(m)
	return ko_module_map


def resolve_module_rpkm(ko_rpkm_fname, ko_module_map) -> dict:
	ret = collections.defaultdict(list)
	with get_fp(ko_rpkm_fname, "r") as fp:
		for line in fp:
			ko, *rpkm_list = line.rstrip().split("\t")
			if ko not in ko_module_map:
				continue
			rpkm = float(rpkm_list[2])  # use sum
			for m in ko_module_map[ko]:
				ret[m].append(float(rpkm))
	return ret


def save_module_rpkm(module_rpkm, fname):
	with get_fp(fname, "w") as fp:
		fields = ["#module", "median", "mean", "sum"]
		print(("\t").join(fields), file=fp)
		for m in sorted(module_rpkm.keys()):
			v = module_rpkm[m]
			fields = [
				m,
				"%.4f" % numpy.median(v),
				"%.4f" % numpy.mean(v),
				"%.4f" % sum(v),
			]
			print(("\t").join(fields), file=fp)
	return


def main():
	args = get_args()
	module_ko_map = load_ko_module_map(args.ko_module_map)
	module_rpkm = resolve_module_rpkm(args.ko_rpkm, module_ko_map)
	save_module_rpkm(module_rpkm, args.output)
	return


if __name__ == "__main__":
	main()
