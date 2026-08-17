#!/usr/bin/env python3
"""
Community figure: distance tree | ANI heatmap | CheckM2 quality | plasmid map.

All four panels share one row order, taken from the tree, so a row can be read
straight across.

A note on the tree: it is a UPGMA dendrogram over mash distances, not a
marker-gene phylogeny. Mash distance is defined across the full divergence
range present here, whereas skani reports nothing below ~80% ANI, which would
leave most inter-genus pairs without a value. The dendrogram is used for
ordering and to show community structure; it should be labelled as an ANI/mash
distance tree rather than a phylogeny.

Usage:
    python3 13_plot_community.py --data subset_char --out community_figure
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform


def read_full_matrix(path):
    names, rows = [], []
    with open(path) as fh:
        first = fh.readline().strip()
        try:
            int(first)
        except ValueError:
            fh.seek(0)
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            p = line.split("\t")
            names.append(os.path.basename(p[0]).replace(".fna", ""))
            rows.append([float(x) for x in p[1:] if x.strip() != ""])
    return names, np.array(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", default="community_figure")
    ap.add_argument("--min-hosts", type=int, default=2,
                    help="plasmid groups shown as columns must occur in at "
                         "least this many genomes")
    ap.add_argument("--label-len", type=int, default=42)
    ap.add_argument("--heatmap", choices=["mash", "skani"], default="mash",
                    help="mash: ANI estimated as 100*(1-mash distance), "
                         "defined across the full divergence range. "
                         "skani: alignment-based ANI, but it reports nothing "
                         "below ~80%% so most inter-genus cells are empty.")
    ap.add_argument("--vmin", type=float, default=None)
    args = ap.parse_args()

    D = args.data
    P = lambda f: os.path.join(D, f)

    # ------------------------------------------------------------- inputs
    names, ani = read_full_matrix(P("ani_matrix.tsv"))
    idx = {n: i for i, n in enumerate(names)}
    n = len(names)
    print(f"ANI matrix: {n} genomes")

    # mash distances -> full matrix
    mash = np.zeros((n, n))
    with open(P("mash_dist.tsv")) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            a = os.path.basename(f[0]).replace(".fna", "")
            b = os.path.basename(f[1]).replace(".fna", "")
            if a in idx and b in idx:
                mash[idx[a], idx[b]] = float(f[2])
    mash = (mash + mash.T) / 2
    np.fill_diagonal(mash, 0.0)

    org, size, npl = {}, {}, {}
    with open(P("organisms.tsv")) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            org[f[0]], size[f[0]], npl[f[0]] = f[1], float(f[2]), int(f[4])

    comp, cont = {}, {}
    qp = P("checkm2/quality_report.tsv")
    if os.path.exists(qp):
        with open(qp) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ci = hdr.index("Completeness") if "Completeness" in hdr else 1
            xi = hdr.index("Contamination") if "Contamination" in hdr else 2
            for line in fh:
                f = line.rstrip("\n").split("\t")
                comp[f[0]] = float(f[ci])
                cont[f[0]] = float(f[xi])
        print(f"CheckM2: {len(comp)} genomes")
    else:
        print("CheckM2 report absent -- quality panel will be blank",
              file=sys.stderr)

    pg = {}
    hosts_of = {}
    gp = P("plasmid_groups.tsv")
    if os.path.exists(gp):
        with open(gp) as fh:
            next(fh)
            for line in fh:
                g, _c, grp, gs, _l = line.rstrip("\n").split("\t")
                pg.setdefault(g, set()).add(grp)
                hosts_of.setdefault(grp, set()).add(g)
    shared = sorted((g for g, h in hosts_of.items()
                     if len(h) >= args.min_hosts),
                    key=lambda g: (-len(hosts_of[g]), g))
    print(f"plasmid groups: {len(hosts_of)} total, {len(shared)} in "
          f">={args.min_hosts} genomes")

    # ------------------------------------------------------------- ordering
    Z = linkage(squareform(mash, checks=False), method="average")
    dend = dendrogram(Z, no_plot=True, labels=names)
    order = [names.index(l) for l in dend["ivl"]]
    onames = [names[i] for i in order]

    # ------------------------------------------------------------- layout
    ncol_pl = max(len(shared), 1)
    fig_w = 14.0 + max(3, ncol_pl) * 0.26
    fig = plt.figure(figsize=(fig_w, 0.34 * n + 2.4), dpi=200)
    left = 0.035
    w_tree, w_lab, w_ani, w_q = 0.10, 0.20, 0.30, 0.055
    w_pl = max(0.055, 0.013 * ncol_pl)
    gap = 0.012
    bottom, height = 0.085, 0.845

    x = left
    ax_t = fig.add_axes([x, bottom, w_tree, height]); x += w_tree
    ax_l = fig.add_axes([x, bottom, w_lab, height]);  x += w_lab + gap
    ax_a = fig.add_axes([x, bottom, w_ani, height]);  x += w_ani + gap
    ax_c = fig.add_axes([x, bottom, w_q, height]);    x += w_q + 0.006
    ax_x = fig.add_axes([x, bottom, w_q, height]);    x += w_q + gap
    ax_p = fig.add_axes([x, bottom, w_pl, height])

    # ------------------------------------------------------------- tree
    dendrogram(Z, orientation="left", ax=ax_t, no_labels=True,
               color_threshold=0, above_threshold_color="#43505e",
               link_color_func=lambda k: "#43505e")
    ax_t.invert_yaxis()
    ax_t.set_xticks([]); ax_t.set_yticks([])
    for sp in ax_t.spines.values():
        sp.set_visible(False)
    ax_t.set_title("mash distance", fontsize=8, color="#5c6b7a", pad=6)

    # ------------------------------------------------------------- labels
    ax_l.set_xlim(0, 1); ax_l.set_ylim(n, 0)
    ax_l.axis("off")
    ecoli_col, other_col = "#2f6f9f", "#1b2430"
    for r, g in enumerate(onames):
        o = org.get(g, g)
        if len(o) > args.label_len:
            o = o[:args.label_len - 1] + "\u2026"
        is_ec = o.lower().startswith("escherichia coli")
        ax_l.text(0.01, r + 0.5, o, va="center", ha="left", fontsize=7.4,
                  color=ecoli_col if is_ec else other_col,
                  fontweight="bold" if is_ec else "normal", style="italic")
        ax_l.text(0.995, r + 0.5, f"{size.get(g,0):.2f} Mb",
                  va="center", ha="right", fontsize=6.6, color="#8d99a6")

    # ------------------------------------------------------------- ANI
    # skani is an alignment-based ANI but stops reporting below ~80%, which
    # leaves most inter-genus cells empty in a community spanning 16 genera.
    # Mash distance is defined across the whole range, so 100*(1-d) fills the
    # matrix and makes genus-level structure visible. The E. coli block is
    # identical either way.
    if args.heatmap == "mash":
        mat = 100.0 * (1.0 - mash)
        vmin = args.vmin if args.vmin is not None else float(
            np.percentile(mat[mat < 99.9], 2))
        title = "estimated ANI (%), mash"
        note = None
    else:
        mat = ani
        vmin = args.vmin if args.vmin is not None else 80.0
        title = "pairwise ANI (%), skani"
        note = "grey = below detection"
    sub = mat[np.ix_(order, order)]
    masked = np.ma.masked_where(sub <= 0, sub)
    cmap_a = plt.get_cmap("YlGnBu").copy()
    cmap_a.set_bad("#f2f4f6")
    im = ax_a.imshow(masked, cmap=cmap_a, vmin=vmin, vmax=100,
                     aspect="auto", interpolation="nearest")
    ax_a.set_xticks([]); ax_a.set_yticks([])
    ax_a.set_title(title, fontsize=8, color="#5c6b7a", pad=6)
    for sp in ax_a.spines.values():
        sp.set_color("#c3ccd6")
    cax = fig.add_axes([ax_a.get_position().x0,
                        bottom - 0.055, w_ani * 0.45, 0.014])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    ticks = [t for t in (70, 75, 80, 85, 90, 95, 100) if t >= vmin - 1]
    cb.set_ticks(ticks[:: max(1, len(ticks) // 4)])
    cb.ax.tick_params(labelsize=6.5, length=2)
    cb.outline.set_visible(False)
    if note:
        cax.text(1.05, 0.5, note, transform=cax.transAxes,
                 fontsize=6.5, va="center", color="#8d99a6")

    # ------------------------------------------------------------- CheckM2
    def qpanel(ax, vals, cmap, vmin, vmax, title, fmt="{:.0f}"):
        arr = np.array([[vals.get(g, np.nan)] for g in onames])
        m = np.ma.masked_invalid(arr)
        c = plt.get_cmap(cmap).copy(); c.set_bad("#f2f4f6")
        ax.imshow(m, cmap=c, vmin=vmin, vmax=vmax, aspect="auto",
                  interpolation="nearest")
        for r, g in enumerate(onames):
            v = vals.get(g)
            if v is None:
                continue
            frac = (v - vmin) / (vmax - vmin + 1e-9)
            ax.text(0, r, fmt.format(v), ha="center", va="center", fontsize=6,
                    color="white" if frac > 0.55 else "#1b2430")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(title, fontsize=7.4, color="#5c6b7a", pad=6)
        for sp in ax.spines.values():
            sp.set_color("#c3ccd6")

    qpanel(ax_c, comp, "Greens", 90, 100, "compl.\n(%)", "{:.1f}")
    qpanel(ax_x, cont, "OrRd", 0, 5, "contam.\n(%)", "{:.1f}")

    # ------------------------------------------------------------- plasmids
    if shared:
        mat = np.array([[1.0 if pgrp in pg.get(g, ()) else 0.0
                         for pgrp in shared] for g in onames])
        cmap_p = LinearSegmentedColormap.from_list(
            "pl", ["#f2f4f6", "#c8801f"])
        ax_p.imshow(mat, cmap=cmap_p, vmin=0, vmax=1, aspect="auto",
                    interpolation="nearest")
        ax_p.set_xticks(range(len(shared)))
        ax_p.set_xticklabels(
            [f"{g} ({len(hosts_of[g])})" for g in shared],
            rotation=90, fontsize=5.8, color="#5c6b7a")
        ax_p.set_yticks([])
        ax_p.set_title("shared\nplasmid groups", fontsize=7.4,
                       color="#5c6b7a", pad=6)
        for sp in ax_p.spines.values():
            sp.set_color("#c3ccd6")
        for r in range(len(onames) + 1):
            ax_p.axhline(r - 0.5, color="white", lw=0.5)
        for c in range(len(shared) + 1):
            ax_p.axvline(c - 0.5, color="white", lw=0.5)
    else:
        ax_p.axis("off")

    # total plasmid count per genome, at the far right
    xr = ax_p.get_position().x1 + 0.016
    ax_n = fig.add_axes([xr, bottom, 0.022, height])
    ax_n.set_xlim(0, 1); ax_n.set_ylim(n, 0); ax_n.axis("off")
    ax_n.set_title("all\nplasmids", fontsize=7.4, color="#5c6b7a", pad=6)
    for r, g in enumerate(onames):
        ax_n.text(0.5, r + 0.5, str(npl.get(g, 0)), ha="center", va="center",
                  fontsize=6.8, color="#1b2430")

    fig.suptitle("Hy-SCALE in silico validation community "
                 f"(n = {n} Enterobacterales genomes)",
                 fontsize=12, fontweight="bold", color="#1b2430",
                 x=left, ha="left", y=0.975)

    for ext in ("png", "svg", "pdf"):
        fig.savefig(f"{args.out}.{ext}", bbox_inches="tight",
                    facecolor="white")
    print(f"wrote {args.out}.png / .svg / .pdf")


if __name__ == "__main__":
    main()
