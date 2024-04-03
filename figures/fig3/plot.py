#!/usr/bin/env python3

import collections
import dataclasses
import json
import re
import sys
import typing

import matplotlib
import matplotlib.patches
import matplotlib.pyplot
import mpllayout
import numpy
import scipy.cluster
import scipy.spatial
import scipy.stats
import sklearn.decomposition
import Bio.SeqIO


def create_layout():
	lc = mpllayout.LayoutCreator(
		left_margin=0.6,
		right_margin=1.5,
		top_margin=0.2,
		bottom_margin=0.5,
	)

	ppk_j = lc.add_frame("ppk_j")
	ppk_j.set_anchor("bottomleft")
	ppk_j.set_size(1.5, 3.0)

	ppk_b = lc.add_frame("ppk_b")
	ppk_b.set_anchor("bottomleft", ppk_j, "bottomright", offsets=(0, 0))
	ppk_b.set_size(3.6, 3.0)

	ppk = lc.add_frame("ppk")
	ppk.set_anchor("bottomleft", ppk_j, "topleft", offsets=(0, 0.5))
	ppk.set_size(1.8, 2.0)

	ppx = lc.add_frame("ppx")
	ppx.set_anchor("bottomleft", ppk, "bottomright", offsets=(2.0, 0))
	ppx.set_size(1.8, 2.0)

	layout = lc.create_figure_layout()

	return layout


class NonNegFloat(float):
	def __new__(cls, value):
		self = super().__new__(cls, value)
		if self < 0:
			raise ValueError("value must be non-negative")
		return self


class JTree(object):
	unknown_taxonomy = 0x7fffffff

	@dataclasses.dataclass
	class JTreeNode(object):
		name: typing.Optional[str] = None
		dist: NonNegFloat = 0
		edge_id: typing.Optional[int] = None
		children: list = dataclasses.field(default_factory=lambda: list())
		_parent = None

		def __hash__(self) -> int:
			return id(self)

		@property
		def parent(self):
			return self._parent

		@property
		def is_root(self) -> bool:
			return self._parent is None

		def add_child(self, node):
			self.children.append(node)
			node._parent = self
			return

		@classmethod
		def parse(cls, s: str):
			new = cls()
			last_comma = 0
			stack = 0
			for i, c in enumerate(s):
				if c == "(":
					stack += 1
				elif ((c == ",") or c == ")") and (stack == 1):
					new.add_child(cls.parse(s[last_comma + 1: i]))
					last_comma = i

				if c == ")":
					stack -= 1
				if stack == 0:
					break

			if c == ")":
				i += 1

			m = re.match(r"(.*):(.*)\{(.*)\}(\[\d+\])?;?$", s[i:])
			name, dist, edge_id, _ = m.groups()

			new.name = name if name else None
			new.dist = float(dist)
			new.edge_id = int(edge_id)

			return new

		def traverse(self):
			yield self
			for c in self.children:
				yield from c.traverse()
			return

		@property
		def n_leaves(self) -> int:
			return sum([node.is_leaf for node in self.traverse()])

		@property
		def max_leaf_depth(self) -> float:
			if self.is_leaf:
				return 0
			return max([n.dist + n.max_leaf_depth for n in self.children])

		@property
		def min_leaf_depth(self) -> float:
			if self.is_leaf:
				return 0
			return min([n.dist + n.min_leaf_depth for n in self.children])

		@property
		def is_leaf(self):
			return not (self.children)

		def path_to_child(self, node) -> typing.Optional[list]:
			path = [node]
			if node is self:
				return path
			while node.parent is not None:
				path.append(node.parent)
				if node.parent is self:
					return path[::-1]
				node = node.parent
			else:
				return None
			return

		def dist_to_child(self, node) -> typing.Optional[float]:
			path = self.path_to_child(node)
			if path is None:
				return None
			dist = 0.
			for n in path[1:]:
				dist += n.dist
			return dist

		def reroot(self):
			pre_parent = self.parent
			if pre_parent is None:
				return
			if pre_parent.parent is not None:
				pre_parent.reroot()
			pre_parent.children = [i for i in self.parent.children
				if i is not self]
			pre_parent._parent = self
			pre_parent.dist = self.dist
			self._parent = None
			self.children.append(pre_parent)
			self.dist = 0
			return

	def __init__(self, root: JTreeNode, *ka, **kw):
		super().__init__(*ka, **kw)
		self.root = root
		self.renew_map_name_to_node()
		self.renew_map_edge_to_node()
		return

	def renew_map_name_to_node(self):
		d = dict()
		for node in self.traverse():
			if node.name:
				d[node.name] = node
		self._name_to_node = d
		return

	def renew_map_edge_to_node(self):
		d = dict()
		for node in self.traverse():
			d[node.edge_id] = node
		self._edge_to_node = d
		return

	def get_node_by_name(self, name: str) -> JTreeNode:
		return self._name_to_node[name]

	def get_node_by_edge_id(self, edge_id: int) -> JTreeNode:
		return self._edge_to_node[edge_id]

	@classmethod
	def parse(cls, s: str, *, node_cls=JTreeNode):
		new = cls(root=node_cls.parse(s))
		return new

	def traverse(self):
		yield from self.root.traverse()

	def get_path_to_node(self, node: JTreeNode) -> typing.Optional[list]:
		return self.root.path_to_child(node)

	def lowest_common_ancestor(self, *nodes) -> JTreeNode:
		paths = [self.get_path_to_node(i) for i in nodes]
		if any([i is None for i in paths]):
			raise ValueError("all nodes must be from the same tree")
		ret = self.root
		for i in zip(*paths):
			uniques = set(i)
			if len(uniques) != 1:
				break
			ret = uniques.pop()
		return ret

	def reroot(self, node):
		node.reroot()
		self.root = node
		return


def load_seqid_taxonomy(samples: list) -> dict:
	ret = dict()
	for sample in samples:
		v = dict()
		for gene in ("ppk", "ppx"):
			with open("data/%s.%s.tax" % (sample, gene), "r") as fp:
				for line in fp:
					fields = line.strip().split("\t")
					seqid = fields[0]
					tax = fields[1].replace("Candidatus ", "Ca_")
					v[seqid] = tax
		ret[sample] = v
	return ret


def load_filtered_cds_tpm_data(seqid_tax) -> dict:
	ret = dict()
	for sample in seqid_tax.keys():
		v = dict()
		with open("../../sample/%s/analysis/cds_tpm.tsv" % sample, "r") as fp:
			for line in fp:
				seqid, tpm = line.strip().split("\t")
				if seqid in seqid_tax[sample]:
					v[seqid] = float(tpm)
		ret[sample] = v
	return ret


def calc_tpm_table(samples: list, seqid_tax: dict, seqid_tpm: dict, gene: str):
	# get targeted seqids from fasta, for each sample
	tax_tpm = dict()
	for sample in samples:
		v = collections.defaultdict(float)
		fasta_file = "reannotate/%s.faa" % gene
		seqio = Bio.SeqIO.parse(fasta_file, format="fasta")
		for seq in seqio:
			seqid = seq.id
			if seqid in seqid_tpm[sample]:
				v[seqid_tax[sample][seqid]] += seqid_tpm[sample][seqid]
		tax_tpm[sample] = v
	# make into a table
	taxa = sorted(set(t for v in tax_tpm.values() for t in v.keys()))
	taxa = [t for t in taxa if t[0].isupper()]
	data = numpy.zeros((len(samples), len(taxa)), dtype=float)
	for i, sample in enumerate(samples):
		for j, tax in enumerate(taxa):
			data[i, j] = tax_tpm[sample].get(tax, 0)
	return taxa, data


class ColorDict(dict):
	def __init__(self, *ka, proto_list=None, **kw):
		super().__init__(*ka, **kw)
		self._colors = self.get_color_list(proto_list)
		return

	def get_color_list(self, proto_list=None) -> list:
		if proto_list is None:
			proto_list = ["tab10", "Accent", "Dark2"]
		colors = list()
		for i in proto_list:
			colors.extend(matplotlib.colormaps.get_cmap(i).colors)
		return [matplotlib.colors.to_hex(c) for c in colors]

	def get_color(self, tax: str):
		if tax not in self:
			self[tax] = self._colors[len(self)]
		return self[tax]


def plot_bar(layout: dict, gene: str, samples: list, taxa: list,
		tax_tpm: numpy.ndarray, sample_cfg: dict, bar_colors: ColorDict,
		label: str = None):
	print(tax_tpm.sum(axis=1))
	# find the top n taxa, sort by mean
	top_n = 10
	top_idx = numpy.mean(tax_tpm, axis=0).argsort()[::-1][:top_n]

	plot_taxa = [taxa[i] for i in top_idx]
	plot_tpm = tax_tpm[:, top_idx]

	# plot stacked bars
	ubar = layout[gene]

	bottom = numpy.zeros(len(samples), dtype=float)
	handles = []
	for i, tax in enumerate(plot_taxa):
		c = bar_colors.get_color(tax)
		height = plot_tpm[:, i]
		bar = ubar.bar(samples, height, bottom=bottom, width=0.9,
			edgecolor="none", facecolor=c, label=tax
		)
		handles.append(bar)
		bottom += height

	# add legend
	ubar.legend(handles=handles, loc=2, bbox_to_anchor=(0.98, 1.02),
		handlelength=0.8, fontsize=8, frameon=False
	)

	# misc
	ubar.spines["top"].set_visible(False)
	ubar.spines["right"].set_visible(False)
	ubar.tick_params(bottom=False)
	ubar.set_ylabel("TPM (%s)" % gene)
	for l in ubar.get_xticklabels():
		l.set(color=sample_cfg[l.get_text()]["color"],
		fontweight="bold")

	# label
	if label:
		ubar.text(0.02, 0.97, label, fontsize=12, fontweight="bold",
			transform=ubar.transAxes,
			horizontalalignment="left", verticalalignment="top")
	return


def load_seqid_clade(fname: str) -> dict:
	ret = dict()
	with open(fname, "r") as fp:
		for line in fp:
			fields = line.strip().split("\t")
			ret[fields[0].replace(".", "_")] = fields[1]
	return ret


def _recurse_assign_tax(node: JTree.JTreeNode):
	if hasattr(node, "taxonomy"):
		return getattr(node, "taxonomy", None)

	children_tax = set(_recurse_assign_tax(child) for child in node.children)

	if None in children_tax:
		children_tax.remove(None)

	if not children_tax:
		return None

	if JTree.unknown_taxonomy in children_tax:
		tax = JTree.unknown_taxonomy
	elif len(children_tax) == 1:
		node.taxonomy = children_tax.pop()
		tax = node.taxonomy
	else:
		tax = JTree.unknown_taxonomy
	node.taxonomy = tax
	return tax


def _recurse_tax_collapse(collapse_dict: dict, node: JTree.JTreeNode):
	if not hasattr(node, "taxonomy"):
		return
	if node.taxonomy != JTree.unknown_taxonomy:
		collapse_dict[node.edge_id] = node.taxonomy
	else:
		for child in node.children:
			_recurse_tax_collapse(collapse_dict, child)
	return


def calc_jtree_node_collapse(jtree: JTree, seqid_clade: dict) -> list:
	for seqid, tax in seqid_clade.items():
		node = jtree.get_node_by_name(seqid)
		node.taxonomy = tax
	_recurse_assign_tax(jtree.root)
	collapse_dict = dict()
	_recurse_tax_collapse(collapse_dict, jtree.root)
	return collapse_dict


def _recurse_place_nodes(node: JTree.JTreeNode, d: float, h0: float,
		h1: float):
	node.d = d + node.dist
	node.h0 = h0
	node.h1 = h1
	node.h = (h0 + h1) / 2
	if node.is_leaf:
		return

	span = h0 - h1
	# distribute span to children
	c_leaves = [c.n_leaves for c in node.children]
	unit = span / sum(c_leaves)
	h = h1
	for child in node.children:
		c_span = unit * child.n_leaves
		_recurse_place_nodes(child, node.d, h, h + c_span)
		h += c_span
	return


def _recurse_draw_node(axes: matplotlib.axes.Axes, node: JTree.JTreeNode,
		collapse_dict: dict, color_dict: ColorDict):
	# draw the link to parent

	line = matplotlib.lines.Line2D([node.d, node.d - node.dist],
		[node.h, node.h],
		clip_on=False, color="#606060", linewidth=1.5, zorder=2)
	axes.add_line(line)

	if node.edge_id in collapse_dict:
		# draw a collapsed node
		xy = [
			(node.d, node.h),
			(node.d + min(node.max_leaf_depth, 0.2), node.h0),
			(node.d + min(node.min_leaf_depth, 0.2), node.h1),
		]
		color = color_dict.get_color(collapse_dict[node.edge_id])
		patch = matplotlib.patches.Polygon(xy, closed=True, clip_on=False,
			edgecolor=color, facecolor=color, zorder=2)
		axes.add_patch(patch)
	elif not node.is_leaf:
		c_hs = [c.h for c in node.children]
		line = matplotlib.lines.Line2D([node.d, node.d],
			[numpy.min(c_hs), numpy.max(c_hs)],
			clip_on=False, color="#606060", linewidth=1.5, zorder=2)
		axes.add_line(line)
		for child in node.children:
			_recurse_draw_node(axes, child, collapse_dict, color_dict)

	return


def draw_tree(axes: matplotlib.axes.Axes, jtree: JTree, collapse_dict: dict,
		color_dict: ColorDict):
	_recurse_place_nodes(jtree.root, 0, 0, 1)
	_recurse_draw_node(axes, jtree.root, collapse_dict, color_dict)

	# misc
	for sp in axes.spines.values():
		sp.set_visible(False)
	axes.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
	axes.set_xlim(0, 0.38)
	axes.set_ylim(0, 1.0)
	return


def calc_clade_tpm(seqid_tpm: dict, jplace: dict, jtree: JTree) -> tuple:
	clade_tpm = collections.defaultdict(lambda: collections.defaultdict(float))
	for p in jplace["placements"]:
		edge_id = p["p"][0][1]
		tax = getattr(jtree.get_node_by_edge_id(edge_id), "taxonomy", None)
		for seqid, w in p["nm"]:
			sample = seqid[:2]
			clade_tpm[sample][tax] += seqid_tpm[sample][seqid]
	return clade_tpm


def _get_clade_draw_order(jtree: JTree, collapse_dict: dict) -> list:
	edge_ids = sorted(collapse_dict.keys())
	hs = [jtree.get_node_by_edge_id(i).h for i in edge_ids]
	sort_idx = numpy.argsort(hs)
	order = [collapse_dict[edge_ids[i]] for i in sort_idx]
	return [i for i in order if i != "Unknown"] + ["Unknown"]  # unknown last


def draw_bars(axes: matplotlib.axes.Axes, samples: list, clade_order: list,
		clade_tpm: dict, sample_cfg: dict, color_dict: ColorDict):
	# plot bars
	bottom = numpy.zeros(len(samples), dtype=float)
	handles = list()
	for clade in clade_order:
		height = numpy.array([clade_tpm[sample][clade] for sample in samples])
		color = color_dict.get_color(clade)
		bar = axes.bar(samples, height, bottom=bottom, width=0.9, clip_on=False,
			edgecolor="#000000", facecolor=color, label=clade)
		handles.append(bar)
		bottom += height

	# add legend
	axes.legend(handles=handles, loc=2, bbox_to_anchor=(0.98, 1.02),
		handlelength=0.8, fontsize=10, frameon=False,
		title="Ca. Accumulibacter clades")

	# misc
	for sp in axes.spines.values():
		sp.set_visible(False)
	axes.tick_params(left=False, bottom=False, right=False, top=False,
		labelleft=False)
	for l in axes.get_xticklabels():
		l.set(color=sample_cfg[l.get_text()]["color"], fontweight="bold")

	# label
	axes.text(0.02, 0.97, "(c)", fontsize=12, fontweight="bold",
		transform=axes.transAxes,
		horizontalalignment="left", verticalalignment="top")
	return


def plot_jplace_tree(layout: dict, jplace_fname: str, seqid_tpm: dict,
		samples: list, sample_cfg: dict):
	# load data
	jplace = json.load(open(jplace_fname, "r"))
	jtree = JTree.parse(jplace["tree"])
	seqid_clade = load_seqid_clade("../../supp/ppk1_cap.ref_tree/clade.tsv")
	collapse_dict = calc_jtree_node_collapse(jtree, seqid_clade)

	# draw tree
	ppk_j = layout["ppk_j"]
	color_dict = ColorDict(proto_list=["Set2", "Dark2"])
	draw_tree(ppk_j, jtree, collapse_dict, color_dict)

	# draw bar
	ppk_b = layout["ppk_b"]
	clade_tpm = calc_clade_tpm(seqid_tpm, jplace, jtree)
	clade_order = _get_clade_draw_order(jtree, collapse_dict)
	draw_bars(ppk_b, samples, clade_order, clade_tpm, sample_cfg, color_dict)
	return


def plot(png, dpi=300):
	samples = ["DU", "MD", "UB", "HD", "SC", "WR"]

	# load data
	sample_cfg = json.load(open("../fig1/sample_cfg.json", "r"))
	seqid_tax = load_seqid_taxonomy(samples)
	seqid_tpm = load_filtered_cds_tpm_data(seqid_tax)
	ppk_tax, ppk_tpm = calc_tpm_table(samples, seqid_tax, seqid_tpm, "ppk")
	ppx_tax, ppx_tpm = calc_tpm_table(samples, seqid_tax, seqid_tpm, "ppx")

	layout = create_layout()
	figure = layout["figure"]

	bar_colors = ColorDict()
	plot_bar(layout, "ppk", samples, ppk_tax, ppk_tpm, sample_cfg, bar_colors,
		label="(a)")
	plot_bar(layout, "ppx", samples, ppx_tax, ppx_tpm, sample_cfg, bar_colors,
		label="(b)")

	plot_jplace_tree(layout, "data/acc.ppk.fna.mafft.fasta.jplace",
		seqid_tpm, samples, sample_cfg)

	# save fig and clean up
	figure.savefig(png, dpi=dpi)
	matplotlib.pyplot.close(figure)
	return


def main():
	plot("fig3.png")
	return


if __name__ == "__main__":
	main()
