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


def create_layout(nrow: int, ncol: int):
	lc = mpllayout.LayoutCreator(
		left_margin=0.6,
		right_margin=0.5,
		top_margin=0.2,
		bottom_margin=0.5,
	)

	bartb_w = 0.60 * ncol
	bartb_h = 0.20 * nrow

	bartb = lc.add_frame("bartb")
	bartb.set_anchor("bottomleft", offsets=(0.5, 0))
	bartb.set_size(bartb_w, bartb_h)

	diag = lc.add_frame("diag")
	diag.set_anchor("bottomleft", bartb, "topleft", offsets=(-0.3, 0.6))
	diag.set_size(2.5, 2.5)

	pca = lc.add_frame("pca")
	pca.set_anchor("bottomleft", diag, "bottomright", offsets=(1.0, 0))
	pca.set_size(2.5, 2.5)

	layout = lc.create_figure_layout()
	layout["bartb_w"] = bartb_w
	layout["bartb_h"] = bartb_h
	layout["nrow"] = nrow
	layout["ncol"] = ncol

	return layout


def load_ko_tpm_data(samples) -> dict:
	ret = dict()
	for s in samples:
		v = dict()
		with open("../../sample/%s/analysis/ko_tpm.tsv" % s, "r") as fp:
			for line in fp:
				if line.startswith("#"):
					continue
				ko, *tpm_list = line.strip().split("\t")
				tpm = float(tpm_list[2])  # use sum
				v[ko] = tpm
		ret[s] = v
	return ret


def load_module_tpm_data(samples) -> dict:
	ret = dict()
	for s in samples:
		v = dict()
		with open("../../sample/%s/analysis/module_tpm.tsv" % s, "r") as fp:
			for line in fp:
				if line.startswith("#"):
					continue
				module, *tpm_list = line.strip().split("\t")
				tpm = float(tpm_list[0])  # use median
				v[module] = tpm
		ret[s] = v
	return ret


def calc_pathway_rows(pathway_cfg: dict) -> int:
	ret = 0
	for cfg in pathway_cfg:
		ret += len(cfg.get("modules", [])) + len(cfg.get("genes", []))
	return ret


def calc_multi_ko_tpm(samples: list, kos: list, ko_tpm: dict):
	return [numpy.mean([ko_tpm[s].get(ko, 0) for ko in kos]) for s in samples]


def add_bars(axes: matplotlib.axes.Axes, y: float, tpm: list, tpm_ref: float,
		color: list, show_edge=False, fill_alpha="20"):
	gap = 0.1
	h = 1.0 - gap * 2
	for i, v in enumerate(tpm):
		w = min(v / tpm_ref, 1.0)
		facecolor = color[i] + fill_alpha
		edgecolor = color[i] if show_edge else "none"
		p = matplotlib.patches.Rectangle((i, y + gap), w, h,
			linewidth=0.5, facecolor=facecolor, edgecolor=edgecolor,
			zorder=2)
		axes.add_patch(p)
		axes.text(i + 0.5, y + 1 - 0.5, "%.1f" % v,
			fontsize=6, color="#000000", zorder=3,
			horizontalalignment="center", verticalalignment="center")
	return


def plot_diag(layout: dict, sample_cfg: dict, ko_tpm: dict) -> list:
	cate = collections.defaultdict(list)
	for s, v in sample_cfg.items():
		cate[v["category"]].append(s)
	cate_list = list(cate.keys())

	all_kos = set()
	for v in ko_tpm.values():
		all_kos.update(v.keys())

	tpm_max = max([max(v.values()) for v in ko_tpm.values()])
	tpm_max = numpy.ceil(tpm_max / 200) * 200

	# plot boxes
	sig_kos = dict()
	diag = layout["diag"]
	diag.set_xscale("symlog")
	diag.set_yscale("symlog")
	for ko in all_kos:
		vals = list()
		for s in cate_list:
			vals.append([ko_tpm[s].get(ko, 0) for s in cate[s]])
		vals = numpy.array(vals)
		# t-test
		test_res = scipy.stats.ttest_ind(vals[0], vals[1], equal_var=False)
		sig = test_res.pvalue < 0.01
		if sig:
			sig_kos[ko] = {k: v for k, v in zip(cate_list, vals.flatten())}
		# plot box
		xy = numpy.min(vals, axis=1)
		w, h = numpy.ptp(vals, axis=1)
		p = matplotlib.patches.Rectangle(xy, w, h, linewidth=0.5,
			edgecolor="#ff0000" if sig else "none",
			facecolor="#ffc0c0" if sig else "#c0c0ff",
			alpha=0.5 if sig else 0.2,
			zorder=3 if sig else 2)
		diag.add_patch(p)

	# add diagonal line
	diag.plot([0, tpm_max], [0, tpm_max], linestyle="--", linewidth=1.0,
		color="#a0a0a0", zorder=3)

	# misc
	diag.spines["top"].set_visible(False)
	diag.spines["right"].set_visible(False)
	diag.set_xlim(0, tpm_max)
	diag.set_ylim(0, tpm_max)
	diag.set_xlabel("KO TPM (EBPR)", fontsize=10, fontweight="bold",
		color=sample_cfg["DU"]["color"])
	diag.set_ylabel("KO TPM (S2EBPR)", fontsize=10, fontweight="bold",
		color=sample_cfg["HD"]["color"])

	# label
	diag.text(0.03, 0.97, "(a)", fontsize=12, fontweight="bold",
		transform=diag.transAxes,
		horizontalalignment="left", verticalalignment="top")

	# save sig kos
	json.dump(sig_kos, open("sig_kos.json", "w"), sort_keys=True, indent="\t")
	with open("sig_kos.tsv", "w") as ofp:
		print(("\t").join(["#ko"] + cate_list), file=ofp)
		for k in sorted(sig_kos.keys()):
			print(("\t").join([k] + [str(sig_kos[k][c]) for c in cate_list]),
				file=ofp)
	return


def plot_pca(layout: dict, sample_cfg: dict, ko_tpm: dict):
	cate = collections.defaultdict(list)
	cate_color = dict()
	for s, v in sample_cfg.items():
		cate[v["category"]].append(s)
		cate_color[v["category"]] = v["color"]
	cate_list = list(cate.keys())
	samples = sorted(sample_cfg.keys())
	n_samples = len(samples)

	all_kos = set()
	for v in ko_tpm.values():
		all_kos.update(v.keys())
	all_kos = sorted(all_kos)

	# prepare data
	data = numpy.empty((n_samples, len(all_kos)), dtype=float)
	for i, s in enumerate(samples):
		for j, k in enumerate(all_kos):
			data[i, j] = ko_tpm[s].get(k, 0)

	# pca
	pca = sklearn.decomposition.PCA(n_components=2)
	transformed = pca.fit_transform(data)

	# plot
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

	# misc
	axes.axhline(0, linestyle="--", linewidth=1.0, color="#a0a0a0", zorder=0)
	axes.axvline(0, linestyle="--", linewidth=1.0, color="#a0a0a0", zorder=0)
	axes.set_xlabel("PC1 (%.2f%%)" % (pca.explained_variance_ratio_[0] * 100),
		fontweight="bold")
	axes.set_ylabel("PC2 (%.2f%%)" % (pca.explained_variance_ratio_[1] * 100),
		fontweight="bold")

	# label
	axes.text(0.03, 0.97, "(b)", fontsize=12, fontweight="bold",
		transform=axes.transAxes,
		horizontalalignment="left", verticalalignment="top")
	return


def plot_bartb(layout: dict, sample_cfg: dict, pathway_cfg: list,
		samples: list, module_tpm: dict, ko_tpm: dict):
	tpm_ref = 800.0

	bartb = layout["bartb"]
	bartb_w = layout["bartb_w"]
	bartb_h = layout["bartb_h"]
	nrow = layout["nrow"]
	ncol = layout["ncol"]

	for sp in bartb.spines.values():
		sp.set_visible(False)

	# add bars by pathway cfg
	sec_y: float = 0
	curr_y: float = 0
	color = [sample_cfg[s]["color"] for s in samples]
	for p_cfg in pathway_cfg:
		curr_y = sec_y
		# add horizontal line
		bartb.plot([-0.8, ncol], [nrow - curr_y] * 2, clip_on=False,
			linestyle="-", linewidth=1.0, color="#000000", zorder=5)

		# add modules
		for m in p_cfg.get("modules", []):
			# add bar
			tpm = [module_tpm[s].get(m["mo"], 0) for s in samples]
			add_bars(bartb, nrow - curr_y - 1, tpm, tpm_ref, color,
				show_edge=False, fill_alpha="d0")
			# add text
			test_res = scipy.stats.ttest_ind(tpm[:3], tpm[3:], equal_var=False)
			text_color = "#ff0000" if test_res.pvalue < 0.05 else "#000000"
			bartb.text(-0.05, nrow - curr_y - 0.5, m.get("display", m["name"]),
				fontsize=7, fontweight="bold", color=text_color,
				horizontalalignment="right", verticalalignment="center"
			)
			bartb.text(ncol + 0.05, nrow - curr_y - 0.5, m.get("description", ""),
				fontsize=7, color=text_color,
				horizontalalignment="left", verticalalignment="center"
			)
			# increment y
			curr_y += 1
		# add genes
		for g in p_cfg.get("genes", []):
			# add bar
			tpm = calc_multi_ko_tpm(samples, g["kos"], ko_tpm)
			add_bars(bartb, nrow - curr_y - 1, tpm, tpm_ref, color,
				show_edge=True, fill_alpha="60")
			# add text
			test_res = scipy.stats.ttest_ind(tpm[:3], tpm[3:])
			text_color = "#ff4040" if test_res.pvalue < 0.05 else "#404040"
			bartb.text(-0.05, nrow - curr_y - 0.5, g.get("display", g["name"]),
				fontsize=7, color=text_color,
				horizontalalignment="right", verticalalignment="center"
			)
			bartb.text(ncol + 0.05, nrow - curr_y - 0.5, g.get("description", ""),
				fontsize=7, color=text_color,
				horizontalalignment="left", verticalalignment="center"
			)
			# increment y
			curr_y += 1

		# add section label
		bartb.text(-1.1, nrow - (curr_y + sec_y) / 2, p_cfg["name"].upper(),
			rotation=90, fontsize=8, fontweight="bold", color="#000000",
			horizontalalignment="center", verticalalignment="center")

		sec_y = curr_y
	# add last horizontal line
	bartb.plot([-0.8, ncol], [0, 0], clip_on=False,
		linestyle="-", linewidth=1.0, color="#000000", zorder=5)

	# other lines and background
	for i, s in enumerate(samples):
		color = sample_cfg[s]["color"]
		bartb.axvline(i, linestyle="-", linewidth=1.0, clip_on=False,
			color=color, zorder=4)
		p = matplotlib.patches.Rectangle((i, 0), 1, nrow,
			edgecolor="none", facecolor=color + "20", zorder=1)
		bartb.add_patch(p)

	# misc
	bartb.tick_params(left=False, right=False, bottom=False, top=False,
		labelleft=False,
	)
	bartb.set_xlim(0, ncol)
	bartb.set_ylim(0, nrow)
	bartb.set_xticks(numpy.arange(ncol) + 0.5)
	bartb.set_xticklabels(samples, fontsize=10, fontweight="bold")
	for label in bartb.get_xticklabels():
		label.set_color(sample_cfg[label.get_text()]["color"])

	# label
	bartb.text(-0.15, 0.07, "(c)", fontsize=12, fontweight="bold",
		transform=bartb.transAxes,
		horizontalalignment="right", verticalalignment="top")
	return


def plot(png, dpi=300):
	samples = ["DU", "MD", "UB", "HD", "SC", "WR"]

	# load data
	sample_cfg = json.load(open("../fig1/sample_cfg.json", "r"))
	pathway_cfg = json.load(open("pathway.json", "r"))
	ko_tpm = load_ko_tpm_data(samples)
	module_tpm = load_module_tpm_data(samples)

	nrow = calc_pathway_rows(pathway_cfg)
	ncol = len(samples)

	layout = create_layout(nrow, ncol)
	figure = layout["figure"]

	plot_diag(layout, sample_cfg, ko_tpm)
	plot_pca(layout, sample_cfg, ko_tpm)
	plot_bartb(layout, sample_cfg, pathway_cfg, samples, module_tpm, ko_tpm)

	# save fig and clean up
	figure.savefig(png, dpi=dpi)
	matplotlib.pyplot.close(figure)
	return


def main():
	plot("fig2.png")
	return


if __name__ == "__main__":
	main()
