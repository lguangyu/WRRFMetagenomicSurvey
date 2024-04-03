#!/usr/bin/env python3

import collections
import numpy


# rank_name = ["class", "order", "family", "genus"]
# rank_level = [2, 3, 4, 5]
rank_name = ["genus"]
rank_level = [5]


def load_taxa_contig(fname) -> dict:
	ret = collections.defaultdict(lambda: collections.defaultdict(list))
	with open(fname, "r") as fp:
		for line in fp:
			cid, *tax = line.rstrip().split("\t")
			if not tax:
				continue
			tax = tax[0].split(";")
			for i, n in zip(rank_level, rank_name):
				if len(tax) > i:
					ret[n][tax[i]].append(cid)
				else:
					ret[n]["unclassified"].append(cid)
	return ret


def load_contig_rpkm(rpkm_fname) -> dict:
	ret = dict()
	with open(rpkm_fname, "r") as fp:
		for line in fp:
			cid, rpkm = line.rstrip().split("\t")
			ret[cid] = float(rpkm)
	return ret


def calc_taxa_rpkm(taxa_contig, contig_rpkm) -> dict:
	ret = dict()
	for rank, d in taxa_contig.items():
		t = dict()
		for k, v in d.items():
			t[k] = numpy.sum([contig_rpkm[cid] for cid in v])
		ret[rank] = t
	return ret


def load_sample_rank_taxa_rpkm() -> dict:
	ret = dict()
	for sample in ["DU", "HD", "MD", "SC", "UB", "WR"]:
		taxa_contig = load_taxa_contig(
			"../../sample/%s/gene_calling/prodigal.cds.faa.blastp.add_tax.lca" % sample)
		contig_rpkm = load_contig_rpkm(
			"../../sample/%s/analysis/contig_rpkm.tsv" % sample)
		ret[sample] = calc_taxa_rpkm(taxa_contig, contig_rpkm)
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
		with open("summary_blastp_taxa_abund.%s.tsv" % rank, "w") as fp:
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
	sample_rank_taxa_rpkm = load_sample_rank_taxa_rpkm()
	save_sample_rank_taxa_rpkm(sample_rank_taxa_rpkm)
	return


if __name__ == "__main__":
	main()
