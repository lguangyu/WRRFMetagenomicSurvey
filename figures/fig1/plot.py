#!/usr/bin/env python3

import collections
import json
import sys
import matplotlib
import matplotlib.patches
import matplotlib.pyplot
import mpllayout
import numpy
import scipy.cluster
import scipy.spatial
import scipy.stats
import sklearn.decomposition


def create_layout():
	lc = mpllayout.LayoutCreator(
		left_margin=1.0,
		right_margin=0.2,
		top_margin=0.5,
		bottom_margin=0.8,
	)

	pca = lc.add_frame("pca")
	pca.set_anchor("bottomleft", offsets=(0, 0))
	pca.set_size(3, 3)

	cbar = lc.add_frame("cbar")
	cbar.set_anchor("topleft", pca, "topright", offsets=(0.5, -0.1))
	cbar.set_size(3, 0.1)

	abund = lc.add_frame("abund")
	abund.set_anchor("topleft", cbar, "bottomleft", offsets=(0, -0.25))
	abund.set_size(3, 1.2)

	dendro = lc.add_frame("dendro")
	dendro.set_anchor("bottomleft", abund, "bottomright", offsets=(0, 0))
	dendro.set_anchor("topright", abund, "topright", offsets=(0.6, 0))

	mash = lc.add_frame("mash")
	mash.set_anchor("bottomleft", pca, "topleft", offsets=(0, 1.0))
	mash.set_size(3, 3)

	mcbar = lc.add_frame("mcbar")
	mcbar.set_anchor("bottomleft", mash, "topleft", offsets=(0, 0.25))
	mcbar.set_anchor("topright", mash, "topright", offsets=(0, 0.35))

	mdendro = lc.add_frame("mdendro")
	mdendro.set_anchor("bottomleft", mash, "bottomright", offsets=(0, 0))
	mdendro.set_anchor("topright", mash, "topright", offsets=(1.0, 0))

	alphadiv = lc.add_frame("alphadiv")
	alphadiv.set_anchor("topleft", mdendro, "topright", offsets=(1.0, 0))
	alphadiv.set_anchor("bottomright", cbar, "topright", offsets=(0, 1.0))

	layout = lc.create_figure_layout()

	return layout


def load_data(fname) -> tuple:
	# load data
	raw = numpy.loadtxt(fname, delimiter="\t", dtype=object)
	taxa = raw[0, 1:]
	samples = raw[1:, 0]
	data = raw[1:, 1:].astype(float)
	data /= numpy.sum(data, axis=1, keepdims=True)

	# clean up unclassified
	drop = [(i.startswith("unclassified") or i.startswith("cannot be assigned")
		or i.startswith("unknown"))
		for i in taxa
	]
	drop = numpy.array(drop)
	keep = drop ^ True

	# calculated unclassified taxa abundance
	for sample, row in zip(samples, data):
		print("%s\t%.4f" % (sample, numpy.sum(row[drop])), file=sys.stderr)

	taxa = taxa[keep]
	data = data[:, keep]
	# re-normalize
	data /= numpy.sum(data, axis=1, keepdims=True)

	return taxa, samples, data


def plot_pca(layout: dict, taxa, samples, data, sample_cfg):
	n_samples = len(samples)
	assert n_samples == len(sample_cfg)
	n_taxa = len(taxa)

	# pca
	pca = sklearn.decomposition.PCA(n_components=2, whiten=True)
	transformed = pca.fit_transform(data)
	assert transformed.shape == (n_samples, 2)

	# pca scatter
	axes = layout["pca"]
	edgecolors = [sample_cfg[s]["color"] for s in samples]
	facecolors = ["#ffffff40"] * n_samples
	axes.scatter(transformed[:, 0], transformed[:, 1], marker="o", s=40,
		edgecolors=edgecolors, facecolors=facecolors, zorder=2)

	# add text
	for xy, sample in zip(transformed, samples):
		axes.text(xy[0], xy[1] + 0.001, sample, fontsize=10,
			color=sample_cfg[sample]["color"], zorder=3,
			horizontalalignment="center", verticalalignment="bottom")

	# add convex hull
	groups = dict()
	for cfg in sample_cfg.values():
		groups[cfg["category"]] = cfg["color"]
	for cate, color in groups.items():
		indices = [i for i, sample in enumerate(samples)
			if sample_cfg[sample]["category"] == cate]
		vertices = transformed[indices]
		# find convex hull
		hull = scipy.spatial.ConvexHull(vertices)
		hull_xys = vertices[hull.vertices]
		# plot convex hull
		patch = matplotlib.patches.Polygon(hull_xys, closed=True,
			linestyle="--", edgecolor=color + "80", facecolor=color + "20",
			zorder=1)
		axes.add_patch(patch)
		# add convec hull text
		hull_center = hull_xys.mean(axis=0)
		axes.text(*hull_center, cate, fontsize=10, fontweight="bold",
			color=color, zorder=4,
			horizontalalignment="center", verticalalignment="center")

	# pca misc
	xmin, ymin = numpy.min(transformed, axis=0)
	xmax, ymax = numpy.max(transformed, axis=0)
	xspan = xmax - xmin
	yspan = ymax - ymin

	padding = 0.1
	xlim = (xmin - padding * xspan, xmax + padding * xspan)
	ylim = (ymin - padding * yspan, ymax + padding * yspan)

	axes.axhline(0, linestyle="--", linewidth=1.0, color="#a0a0a0", zorder=0)
	axes.axvline(0, linestyle="--", linewidth=1.0, color="#a0a0a0", zorder=0)
	axes.set_xlabel("PC1 (%.2f%%)" % (pca.explained_variance_ratio_[0] * 100),
		fontweight="bold")
	axes.set_ylabel("PC2 (%.2f%%)" % (pca.explained_variance_ratio_[1] * 100),
		fontweight="bold")
	axes.set_xlim(xlim)
	axes.set_ylim(ylim)

	# add arrow of top 5 norms
	norms = numpy.linalg.norm(pca.components_, axis=0)
	arrwo_indices = numpy.argsort(norms)[::-1][:5]

	# constraint arrow in 40% of xlim and ylim
	constraint = 0.6
	arrow_scale = min(
		[xlim[0] / pca.components_[0, i] * constraint
			for i in arrwo_indices if pca.components_[0, i] < 0] +
		[xlim[1] / pca.components_[0, i] * constraint
			for i in arrwo_indices if pca.components_[0, i] > 0] +
		[ylim[0] / pca.components_[1, i] * constraint
			for i in arrwo_indices if pca.components_[1, i] < 0] +
		[ylim[1] / pca.components_[1, i] * constraint
			for i in arrwo_indices if pca.components_[1, i] > 0]
	)

	# plot arrows
	for i in arrwo_indices:
		xy = pca.components_[:, i] * arrow_scale
		axes.annotate("", xy=xy, xytext=(0, 0), color="#4040ff",
			arrowprops=dict(arrowstyle="->", color="#4040ff", linewidth=1.0))
		axes.text(*xy, taxa[i].replace(" ", "\n"), fontsize=8, color="#4040ff",
			horizontalalignment="center", verticalalignment="center")

	for label in (axes.get_xticklabels() + axes.get_yticklabels()):
		label.set_fontsize(8)

	# label
	axes.text(0.95, 0.95, "(c)", fontsize=12, fontweight="bold",
		transform=axes.transAxes,
		horizontalalignment="right", verticalalignment="top")
	return


def plot_abund(layout: dict, taxa, samples, data, sample_cfg):
	n_samples = len(samples)
	assert n_samples == len(sample_cfg)
	n_taxa = len(taxa)

	# find the top 10 taxa
	filter_size = 15
	mean_abund = numpy.mean(data, axis=0)
	sort_indices = numpy.argsort(mean_abund)[::-1][:filter_size]
	plot_data = data[:, sort_indices]
	plot_taxa = taxa[sort_indices]
	plot_taxa = [i.replace("Candidatus", "Ca.") for i in plot_taxa]

	# add hca and dendro
	linkage = scipy.cluster.hierarchy.linkage(plot_data, method="average",
		optimal_ordering=True)
	dendro = scipy.cluster.hierarchy.dendrogram(linkage,
		orientation="right", no_plot=True)
	leaves = numpy.asarray(dendro["leaves"])

	# plot dendrogram
	dax = layout["dendro"]
	for x, y in zip(dendro["dcoord"], dendro["icoord"]):
		dax.plot(x, y, clip_on=False, linestyle="-", linewidth=1.0,
			color="#4080ff")
	imax = 10 * len(leaves)
	dmax = numpy.max(dendro["dcoord"])
	dax.set_xlim(0, dmax)
	dax.set_ylim(0, imax)
	dax.tick_params(
		left=False, right=False, bottom=False, top=False,
		labelleft=False, labelright=False, labelbottom=False, labeltop=False,
	)
	for sp in dax.spines.values():
		sp.set_visible(False)

	# adjust row order
	plot_samples = samples[leaves]

	# heatmap
	abund = layout["abund"]
	cmap = matplotlib.colormaps.get_cmap("jet")
	p = abund.pcolor(plot_data[leaves], cmap=cmap, vmin=0, vmax=0.2)

	# add colorbar
	cax = layout["cbar"]
	cbar = abund.figure.colorbar(p, cax=cax, orientation="horizontal",
		label="genus relative abundance")
	cax.xaxis.set_ticks_position("top")

	# add text
	abund.tick_params(left=False, right=False, bottom=False, top=False)
	abund.set_ylim(n_samples, 0)
	abund.set_xticks(numpy.arange(len(plot_taxa)) + 0.5)
	abund.set_xticklabels(plot_taxa, rotation=90, fontsize=8)
	abund.set_yticks(numpy.arange(n_samples) + 0.5)
	abund.set_yticklabels(plot_samples, fontsize=8)
	for s, label in zip(plot_samples, abund.get_yticklabels()):
		label.set(color=sample_cfg[s]["color"], fontweight="bold")

	# label
	abund.text(0.97, 0.05, "(d)", fontsize=12, fontweight="bold",
		color="#ffffff", transform=abund.transAxes,
		horizontalalignment="right", verticalalignment="bottom")
	return


def plot_alphadiv(layout: dict, taxa, samples, data, sample_cfg):
	alphadiv = scipy.stats.entropy(data, axis=1)
	assert len(alphadiv) == len(samples)

	# group samples
	cate_color = dict()
	cate_dict = collections.defaultdict(list)
	for i, s in enumerate(samples):
		cate_dict[sample_cfg[s]["category"]].append(i)
		cate_color[sample_cfg[s]["category"]] = sample_cfg[s]["color"]
	cate_order = sorted(cate_dict.keys())

	# add alpha diversity
	divax = layout["alphadiv"]
	for x, c in enumerate(cate_order):
		indices = cate_dict[c]
		xs = [x + 0.35] * len(indices)
		ys = alphadiv[indices]
		print(c, ys, file=sys.stderr)
		color = cate_color[c]
		divax.scatter(xs, ys, marker="o", s=60, linewidths=2.0,
			edgecolor=color, facecolor="#ffffff", zorder=3)
		# add text
		for tx, ty, text in zip(xs, ys, samples[indices]):
			divax.text(tx, ty, "  " + text, fontsize=10, color=color,
				horizontalalignment="left", verticalalignment="center",
				zorder=3)

	# misc
	divax.tick_params(bottom=False)
	divax.set_xlim(0, len(cate_dict))
	ymin = numpy.floor(numpy.min(alphadiv) / 0.5) * 0.5
	ymax = numpy.ceil(numpy.max(alphadiv) / 0.5) * 0.5
	divax.set_ylim(ymin, ymax)
	yticks = numpy.arange(ymin, ymax + 0.01, 0.25)
	divax.set_xticks(numpy.arange(len(cate_order)) + 0.5)
	divax.set_xticklabels(cate_order, rotation=45, fontsize=10)
	for label in divax.get_xticklabels():
		label.set(color=cate_color[label.get_text()], fontweight="bold")
	divax.set_yticks(yticks)
	divax.set_ylabel("Alpha diversity (Shannon)", fontweight="bold")

	# add background
	for i, c in enumerate(cate_order):
		patch = matplotlib.patches.Rectangle((i, ymin), 1, ymax - ymin,
			linewidth=0, edgecolor="none", facecolor=cate_color[c] + "20",
			zorder=2)
		divax.add_patch(patch)

	# label
	divax.text(0.95, 0.02, "(b)", fontsize=12, fontweight="bold",
		transform=divax.transAxes,
		horizontalalignment="right", verticalalignment="bottom")
	return


def load_mash_data(fname: str) -> tuple:
	data = collections.defaultdict(dict)
	with open(fname, "r") as fp:
		for line in fp:
			fields = line.strip().split("\t")
			f_str = "%s (%s)" % tuple(fields[0].split(".")[:2])
			t_str = "%s (%s)" % tuple(fields[1].split(".")[:2])
			data[f_str][t_str] = float(fields[2])
	tags = numpy.asarray(sorted(data.keys()), dtype=object)
	tag_to_id = {tag: i for i, tag in enumerate(tags)}

	arr = numpy.full((len(tags), len(tags)), numpy.nan)
	for f, d in data.items():
		for t, v in d.items():
			arr[tag_to_id[f], tag_to_id[t]] = v
	return tags, arr


def plot_mash(layout: dict, fname: str, sample_cfg: dict):
	# load data
	tags, data = load_mash_data(fname)

	# hca to rearrange order
	sq_dist = scipy.spatial.distance.squareform(data, checks=False)
	linkage = scipy.cluster.hierarchy.linkage(sq_dist, method="average",
		optimal_ordering=True)
	dendro = scipy.cluster.hierarchy.dendrogram(linkage,
		orientation="right", no_plot=True)
	leaves = numpy.asarray(dendro["leaves"])
	plot_tags = tags[leaves]
	plot_data = data[numpy.ix_(leaves, leaves)]

	# plot dendrogram
	dax = layout["mdendro"]
	for x, y in zip(dendro["dcoord"], dendro["icoord"]):
		dax.plot(x, y, clip_on=False, linestyle="-", linewidth=1.0,
			color="#4080ff")
	imax = 10 * len(leaves)
	dmax = numpy.max(dendro["dcoord"])
	dax.set_xlim(0, dmax)
	dax.set_ylim(0, imax)
	dax.tick_params(
		left=False, right=False, bottom=False, top=False,
		labelleft=False, labelright=False, labelbottom=False, labeltop=False,
	)
	for sp in dax.spines.values():
		sp.set_visible(False)

	# plot heatmap
	mash = layout["mash"]
	cmap = matplotlib.colormaps.get_cmap("viridis_r")
	vmax = numpy.ceil(numpy.nanmax(plot_data) / 0.1) * 0.1
	p = mash.pcolor(plot_data, cmap=cmap, vmin=0, vmax=vmax)

	# add colorbar
	cax = layout["mcbar"]
	cbar = mash.figure.colorbar(p, cax=cax, orientation="horizontal",
		label="mash distance")
	cax.xaxis.set_ticks_position("top")

	# add labels
	mash.tick_params(left=False, right=False, bottom=False, top=False)
	ticks = numpy.arange(len(plot_tags)) + 0.5
	mash.set_xticks(ticks)
	mash.set_xticklabels(plot_tags, rotation=90, fontsize=8)
	mash.set_yticks(ticks)
	mash.set_yticklabels(plot_tags, fontsize=8)

	for label in (mash.get_xticklabels() + mash.get_yticklabels()):
		sample = label.get_text().split(" ")[0]
		label.set(color=sample_cfg[sample]["color"], fontweight="bold")

	# label
	mash.text(0.95, 0.05, "(a)", fontsize=12, fontweight="bold",
		color="#ffffff", transform=mash.transAxes,
		horizontalalignment="right", verticalalignment="bottom")
	return


def plot(png, dpi=300):
	# load data
	taxa, samples, data = load_data("summary_blastp_taxa_abund.genus.tsv")
	kjtaxa, kjsamples, kjdata = load_data("summary_kaiju_taxa_abund.genus.tsv")
	sample_cfg = json.load(open("sample_cfg.json"))

	layout = create_layout()
	figure = layout["figure"]

	plot_pca(layout, taxa, samples, data, sample_cfg)
	plot_abund(layout, taxa, samples, data, sample_cfg)
	plot_alphadiv(layout, taxa, samples, data, sample_cfg)
	plot_mash(layout, "../../mash/dist.txt", sample_cfg)

	# save fig and clean up
	figure.savefig(png, dpi=dpi)
	matplotlib.pyplot.close(figure)
	return


def main():
	plot("fig1.png")
	return


if __name__ == "__main__":
	main()
