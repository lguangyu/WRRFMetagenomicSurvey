#!/usr/bin/env python3

import argparse
import io
import sys


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("-b", "--blastp-table", type=str, required=True,
		metavar="tsv",
		help="input blastp table to attach taxonomy info, whose 2nd column "
			"must be NCBI accession (required)")
	ap.add_argument("-t", "--tax-table", type=str, required=True,
		metavar="tsv",
		help="input taxonomy table to attach to blastp table, whose 1st "
			"column must be NCBI accession (required)")
	ap.add_argument("-o", "--output", type=str, default="-",
		metavar="tsv",
		help="output blastp table with attached taxonomy info [stdout]")
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


class AccsTaxMap(dict):
	def __init__(self, *ka, delimiter="\t", n_fields=None, **kw):
		super().__init__(*ka, **kw)
		self.delimiter = delimiter
		self.n_fields = n_fields
		return

	def add_record(self, accs: str, tax: str):
		n_fields = len(tax.split(self.delimiter))
		# ensure every line has the same number of fields
		# otherwise raise an error
		if self.n_fields is None:
			self.n_fields = n_fields
		elif self.n_fields != n_fields:
			raise ValueError("each taxonomy record must have the same number "
				"of fields, but '%s' breaks this convention" % tax)
		self[accs] = tax
		return

	def get(self, key, default=None):
		ret = super().get(key, default)
		if default is None:
			# prevent error when self.n_fields is None
			ret = self.delimiter * ((self.n_fields or 1) - 1)
		return ret

	@classmethod
	def from_file(cls, file, *ka, **kw):
		new = cls(*ka, **kw)
		with get_fp(file, "r") as fp:
			for line in fp:
				accs, tax = line.rstrip().split(new.delimiter, 1)
				new.add_record(accs, tax)
		return new


def add_taxonomy_to_blastp(ifile, ofile, tax_table: dict) -> None:
	with get_fp(ifile, "r") as ifp:
		with get_fp(ofile, "w") as ofp:
			for line in ifp:
				line = line.rstrip()
				_, accs, _ = line.split("\t", 2)
				if accs not in tax_table:
					print("'%s' missing annotation (skipped)" % accs,
						file=sys.stderr)
					continue
				print(line + "\t" + tax_table[accs], file=ofp)
	return


def main():
	args = get_args()
	tax_dict = AccsTaxMap.from_file(args.tax_table)
	add_taxonomy_to_blastp(args.blastp_table, args.output, tax_dict)
	return


if __name__ == "__main__":
	main()
