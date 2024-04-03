#!/usr/bin/env python3

import argparse
import collections
import functools
import io
import math
import os
import sys


RANK_ORD = {
	"superkingdom": 0,
	"phylum": 1,
	"class": 2,
	"order": 3,
	"family": 4,
	"genus": 5,
	"species": 6,
}


def get_args():
	ap = argparse.ArgumentParser()
	ap.add_argument("-c", "--contig-gene-taxid", type=str, required=True,
		metavar="tsv",
		help="individual gene taxonomy annotation in each contig; should be "
			"a two-column table input, the first column is gene header "
			"following the format '<contig_id>_<gene_id>', for example: "
			"'c_000000000027_1' means the #1 gene on contig c_000000000027; "
			"the second column is NCBI's numerical taxonomy id (required)")
	ap.add_argument("-n", "--nodes-dump", type=str, required=True,
		metavar="nodes.dmp",
		help="nodes.dump file (required)")
	ap.add_argument("-a", "--names-dump", type=str, required=True,
		metavar="names.dmp",
		help="names.dump file (required)")
	ap.add_argument("-b", "--bootstrap", type=Fraction, default=0.8,
		metavar="float",
		help="bootstrap value to assign contig taxonomy using LCA algorithm "
			"[0.8]")
	ap.add_argument("-o", "--output", type=str, default="-",
		metavar="tsv",
		help="output table [stdout]")
	ap.add_argument("-t", "--tax-tree", type=str,
		metavar="newick",
		help="if set with a file name, write an aggregated taxonomic tree "
			"from all contigs to it [no]")
	# parse and refine args
	args = ap.parse_args()
	if args.output == "-":
		args.output = sys.stdout
	return args


class Fraction(float):
	def __new__(cls, *ka, **kw):
		new = super().__new__(cls, *ka, **kw)
		if (new < 0) or (new > 1.0):
			raise ValueError("fraction should be a decimal from 0.0 to 1.0, "
				"got '%f'" % new)
		return new


def get_fp(f, *ka, factory=open, **kw):
	if isinstance(f, io.IOBase):
		ret = f
	elif isinstance(f, str):
		ret = factory(f, *ka, **kw)
	else:
		raise TypeError("first arg of get_fp() must be io.IOBase or str, "
			"got '%s'" % type(f).__name__)
	return ret


def auto_init_property(init_factory):
	"""
	decorator factory for properties that auto initialize itself at the first
	time access

	ARGUMENTS
	init_factory:
	  factory to be called to generate the initial value to be called at the
	  first time access; also used to wrap the getter method;
	"""
	_attr_name = "_aip__" + init_factory.__name__

	@property
	@functools.wraps(init_factory)
	def getter(self):
		nonlocal _attr_name
		if not hasattr(self, _attr_name):
			setattr(self, _attr_name, init_factory(self))
		return getattr(self, _attr_name)
	assert isinstance(getter, property), type(getter)
	return getter


class Struct(list):
	"""
	list wrapped as struct, use attribute name to access data as well as index
	"""
	class field(property):
		def __init__(self, index, type_cast=str, doc=None):
			super().__init__(doc=doc,
				fget=lambda self: type_cast(self[index]))
			return


class NCBIDumpFormat(object):
	"""
	handles the NCBI taxonomy database dmp file format:
	delimiter: \\t|\\t
	EOL: \\t|\\n
	"""
	@classmethod
	def from_dumped_line(cls, raw_line):
		fields = raw_line.replace("\t|\n", "").split("\t|\t")
		return cls(fields)


class NCBITaxonomyNodeName(Struct, NCBIDumpFormat):
	"""
	name field definition please refer to:
	ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz/readme.txt
	"""
	# properties
	tax_id = Struct.field(0, type_cast=int)
	name_txt = Struct.field(1)
	unique_name = Struct.field(2)
	name_class = Struct.field(3)

	def __init__(self, fields):
		if len(fields) != 4:
			raise ValueError("bad format: " + str(fields))
		super().__init__(fields)
		return

	@property
	def is_scientific_name(self):
		return self.name_class == "scientific name"


class NCBITaxonomyNodeNameList(list):
	@property
	def scientific_name(self) -> str:
		"""
		scientific name of the node;
		"""
		# this property should point to the name object in self whose name class
		# is 'scientific name'
		return ("NULL" if self._sci_name_obj is None
			else self._sci_name_obj.name_txt)

	def add_name_obj(self, name_obj):
		if not isinstance(name_obj, NCBITaxonomyNodeName):
			raise TypeError("name_obj must be NCBITaxonomyNodeName")
		self.append(name_obj)
		# check name class value, if is 'scientific name', update
		# self._sci_name_obj attribute
		if name_obj.is_scientific_name:
			self._sci_name_obj = name_obj
		return

	def __init__(self, *ka, **kw):
		super().__init__(*ka, **kw)
		self._sci_name_obj = None
		return


class NCBITaxonomyNode(Struct, NCBIDumpFormat):
	"""
	node field definition please refer to:
	ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz/readme.txt
	"""
	# properties
	tax_id = Struct.field(0, type_cast=int)
	parent_tax_id = Struct.field(1, type_cast=int)
	rank = Struct.field(2)
	embl_code = Struct.field(3)
	div_id = Struct.field(4, type_cast=int)
	inh_div_flag = Struct.field(5, type_cast=int)  # 0-1 int
	gc_id = Struct.field(6, type_cast=int, doc="genetic code id")
	inh_gc_flag = Struct.field(7, type_cast=int)  # 0-1 int
	mgc_id = Struct.field(8, type_cast=int, doc="mitochondrial g.c. id")
	inh_mgc_flag = Struct.field(9, type_cast=int)  # 0-1 int
	gbk_hid_flag = Struct.field(10, type_cast=int, doc="genbank hidden flag")
	hit_st_flag = Struct.field(11, type_cast=int, doc="hidden subtree flag")
	comments = Struct.field(12)

	@property
	def parent(self):
		return self._parent_node

	@parent.setter
	def parent(self, value):
		if not isinstance(value, NCBITaxonomyNode):
			raise TypeError("parent must be NCBITaxonomyNode")
		self._parent_node = value
		return

	@auto_init_property
	def name_list(self):
		"""
		list of known names
		"""
		return NCBITaxonomyNodeNameList()

	@property
	def scientific_name(self):
		return self.name_list.scientific_name

	def add_name_obj(self, *ka, **kw):
		self.name_list.add_name_obj(*ka, **kw)
		return

	def __init__(self, fields):
		if len(fields) != 13:
			raise ValueError("bad format: " + str(fields))
		super().__init__(fields)
		return


class NCBITaxonomyDB(dict):
	def __getitem__(self, *ka, **kw) -> NCBITaxonomyNode:
		"""
		make database query by tax_id, return a node (NCBITaxonomyNode)
		"""
		return self.query(*ka, **kw)

	def query(self, tax_id: int = None) -> NCBITaxonomyNode:
		"""
		make database query, return a node (NCBITaxonomyNode)
		"""
		return self.get(tax_id, None)

	def add_node(self, node: NCBITaxonomyNode):
		"""
		add a new node to the database
		"""
		if not isinstance(node, NCBITaxonomyNode):
			raise TypeError("must be NCBITaxonomyNode")
		if node.tax_id in self:
			raise ValueError("tax id %d already exist" % node.tax_id)
		self.update({node.tax_id: node})
		return

	############################################################################
	# database I/O: from dump files
	@classmethod
	def from_dumps(cls, nodes_dump, names_dump):
		"""
		load database from .dmp files

		ARGUMENTS:
		nodes_dump:
		  path to the nodes.dmp file
		names_dump:
		  path to the names.dmp file
		"""
		db = cls()
		print("nodes file: %s" % os.path.abspath(nodes_dump), file=sys.stderr)
		print("names file: %s" % os.path.abspath(names_dump), file=sys.stderr)
		print("loading nodes...", file=sys.stderr)
		with get_fp(nodes_dump, "r") as fp:
			for line in fp:
				node = NCBITaxonomyNode.from_dumped_line(line)
				db.add_node(node)
		print("loading names...", file=sys.stderr)
		with get_fp(names_dump, "r") as fp:
			for line in fp:
				name_obj = NCBITaxonomyNodeName.from_dumped_line(line)
				target = db.query(tax_id=name_obj.tax_id)
				target.add_name_obj(name_obj)
		# now polishing data after load has completed
		print("finalizing...", file=sys.stderr)
		db._link_parents()
		print("done", file=sys.stderr)
		return db

	def _link_parents(self):
		# link each nodes to its parent, node.parent will return the parent node
		# instead of parent node's tax_id (node.parent_tax_id)
		for node in self.values():
			assert isinstance(node, NCBITaxonomyNode), type(node).mro()
			node.parent = self.query(tax_id=node.parent_tax_id)
		return


def parse_contig_name(s: str, delimiter="_"):
	*c, _ = s.split(delimiter)
	return delimiter.join(c)


def parse_contig_gene_annotation(file, delimiter="\t") -> dict:
	ret = collections.defaultdict(list)
	with get_fp(file, "r") as fp:
		for line in fp:
			contig_gene, *taxid = line.rstrip().split(delimiter)
			if not taxid:
				print("'%s' has no annotation" % contig_gene, file=sys.stderr)
				continue  # skip unannotated lines
			contig = parse_contig_name(contig_gene)
			ret[contig].append(int(taxid[0]))
	return ret


def contig_taxonomy_by_lca(contig_gene_taxid, tax_db, *,
		bootstrap=0.8, ranks: dict = None) -> dict:
	contig_tax_lca = dict()
	for contig, gene_taxid in contig_gene_taxid.items():
		# this tree deals with all gene taxonomy in one contig
		gene_tax_tree = Tree()
		# add all gene taxoxnomy to the tree
		for i in gene_taxid:
			node = tax_db.query(i)
			if not node:  # tax_db.query() can return None
				continue
			tax_path = _get_node_tax_path(node, ranks=ranks)
			gene_tax_tree.add_path(tax_path)  # tax_path must be descent
		# now find the desired path using lca and bootstrap threshold
		if not gene_tax_tree:  # when no valid nodes are in the Tree()
			continue
		min_count = math.ceil(bootstrap * len(gene_taxid))
		contig_tax_lca[contig] = gene_tax_tree.get_longest_path(min_count)
	return contig_tax_lca


def _get_node_tax_path(node: NCBITaxonomyNode, ranks: dict = None) -> list:
	if ranks is None:
		path = list()
		while node.parent != node:
			path.append(node.name_list.scientific_name)
			node = node.parent
		path = path[::-1]
	else:
		path = [None] * len(ranks)
		for k, v in ranks.items():
			path[v] = "unknown_" + k
		while node.parent != node:
			if node.rank in ranks:
				path[ranks[node.rank]] = node.name_list.scientific_name
			node = node.parent
	return path


class Tree(object):
	class TreeNode(object):
		def __init__(self, parent, name, *ka, count=0, **kw):
			super().__init__(*ka, **kw)
			self.parent = parent
			self.name = name
			self.children = list()
			self.count = count
			return

		@property
		def name_newick(self) -> str:
			return "\"" + self.name + "\""

		def __str__(self) -> str:
			ret = ""
			if self.children:
				ret += "(" + (",").join([str(c) for c in self.children]) + ")"
			ret += self.name_newick
			return ret

		def get_longest_path(self, min_count=1) -> list:
			"""
			find the longest path from this node that has all nodes with the 
			count attribute larger than the min_cnt value
			"""
			ret = list()
			# add the longest child path
			if self.count >= min_count:
				ret += [self.name]
				if self.children:
					sp = [c.get_longest_path(min_count) for c in self.children]
					ret += sorted(sp, key=lambda x: len(x), reverse=True)[0]
			return ret

	def __init__(self, *ka, **kw):
		super().__init__(*ka, **kw)
		self._node_dict = dict()
		return

	def __str__(self) -> str:
		return self.newick_str

	def __bool__(self) -> bool:
		return bool(self._node_dict)

	@property
	def newick_str(self) -> str:
		"""
		dump the whole tree as a newick string
		"""
		return str(self.root) + ";"  # for a tree, need this trailing ';'

	@property
	def root(self):
		for c in self._node_dict.values():
			if c.parent == None:
				ret = c
				break
		else:
			ret = None
		return ret

	def get_node(self, name):
		return self._node_dict[name]

	def add_node(self, name, parent_name=None):
		# if node already exists, don't do anything
		# can happen when adding two paths sharing (some) ancester(s)
		if name in self._node_dict:
			node = self.get_node(name)
		else:
			parent = None if parent_name is None else self.get_node(parent_name)
			node = self.TreeNode(parent, name)
			self._node_dict[name] = node
			if parent:
				parent.children.append(node)
		node.count += 1
		return node

	def add_path(self, path):
		# path must be in descendent order
		for i, name in enumerate(path):
			self.add_node(name, parent_name=(path[i - 1] if i else None))
		return

	def get_longest_path(self, min_count: int = 0):
		"""
		find the longest path from root that have all nodes have count attribute
		larger than the min_count value
		"""
		return self.root.get_longest_path(min_count)


def save_contig_tax(file, contig_tax: dict):
	with get_fp(file, "w") as fp:
		for c in sorted(contig_tax.keys()):
			print("%s\t%s" % (c, (";").join(contig_tax[c])), file=fp)
	return


def save_aggregated_tax_tree(file, contig_tax: dict):
	tree = Tree()
	for v in contig_tax.values():
		tree.add_path(v)
	with get_fp(file, "w") as fp:
		fp.write(str(tree))
	return


def main():
	args = get_args()
	tax_db = NCBITaxonomyDB.from_dumps(args.nodes_dump, args.names_dump)
	contig_gene_taxid = parse_contig_gene_annotation(args.contig_gene_taxid)
	contig_tax = contig_taxonomy_by_lca(contig_gene_taxid, tax_db,
		bootstrap=args.bootstrap, ranks=RANK_ORD)
	save_contig_tax(args.output, contig_tax)
	if args.tax_tree:
		save_aggregated_tax_tree(args.tax_tree, contig_tax)
	return


if __name__ == "__main__":
	main()
