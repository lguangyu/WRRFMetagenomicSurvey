#!/usr/bin/env python3

import sys
import Bio.SeqIO


def load_taxdump_nodes(fname) -> dict:
	ret = dict()
	with open(fname, "r") as fp:
		for line in fp:
			fields = line.strip().split("\t")
			ret[fields[0]] = (fields[2], fields[4])  # parent, rank
	return ret


def load_taxdump_names(fname) -> dict:
	ret = dict()
	with open(fname, "r") as fp:
		for line in fp:
			fields = line.strip().split("\t")
			if fields[6] == "scientific name":
				ret[fields[0]] = fields[2]
	return ret


def find_taxid_genus_name(taxid, nodes: dict, names: dict):
	while True:
		if (taxid == "1") or (taxid not in nodes):
			return "unknown"
		parent, rank = nodes[taxid]
		if rank == "genus":
			break
		else:
			taxid = parent
	return names.get(taxid, "unknown")


def add_tax_name_stdio() -> dict:
	nodes = load_taxdump_nodes("../../supp/taxdump/nodes.dmp")
	names = load_taxdump_names("../../supp/taxdump/names.dmp")

	for line in sys.stdin:
		seqid, taxid = line.strip().split("\t")
		tax = find_taxid_genus_name(taxid, nodes, names)
		print("%s\t%s" % (seqid, tax), file=sys.stdout)
	return


def main():
	add_tax_name_stdio()
	return


if __name__ == "__main__":
	main()
