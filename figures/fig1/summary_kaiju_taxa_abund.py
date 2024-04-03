#!/usr/bin/env python3

import collections
import numpy


# rank_name = ["class", "order", "family", "genus"]
rank_name = ["genus"]


def load_kaiju_abund(fname) -> dict:
	ret = dict()
	skip_rows = 1
	with open(fname, "r") as fp:
		for line in fp:
			if skip_rows:
				skip_rows -= 1
				continue
			fields = line.rstrip().split("\t")
			perc, tax = fields[1], fields[4]
			ret[tax] = float(perc) / 100
	return ret


def load_sample_rank_taxa_rpkm(ranks: list) -> dict:
	ret = dict()
	for sample in ["DU", "HD", "MD", "SC", "UB", "WR"]:
		ret[sample] = dict()
		for rank in ranks:
			ret[sample][rank] = load_kaiju_abund(
				"../../sample/%s/analysis/kaiju.%s.tsv" % (sample, rank))
	return ret


def save_sample_rank_taxa_rpkm(sample_rank_taxa_rpkm: dict):
	sample_list = sorted(sample_rank_taxa_rpkm.keys())
	rank_list = sorted(sample_rank_taxa_rpkm[sample_list[0]].keys())

	for rank in rank_list:
		# find all taxa
		all_taxa = set()
		for sample in sample_list:
			all_taxa.update(sample_rank_taxa_rpkm[sample][rank].keys())
		all_taxa = sorted(all_taxa)
		# save summary
		with open("summary_kaiju_taxa_abund.%s.tsv" % rank, "w") as fp:
			header = ("\t").join([""] + all_taxa)
			print(header, file=fp)
			for sample in sample_list:
				data = [sample]
				values = sample_rank_taxa_rpkm[sample][rank]
				for tax in all_taxa:
					data.append(str(values.get(tax, 0)))
				print(("\t").join(data), file=fp)
	return


def main():
	sample_rank_taxa_rpkm = load_sample_rank_taxa_rpkm(rank_name)
	save_sample_rank_taxa_rpkm(sample_rank_taxa_rpkm)
	return


if __name__ == "__main__":
	main()
