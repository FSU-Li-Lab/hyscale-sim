#!/usr/bin/env python3
"""
Stage 8: the experiments.

E1   precision / recall / yield vs coverage tier          (clean barcodes)
E2   error decomposition                                  (clean barcodes)
E2b  doublet detection sensitivity, stratified by ANI     (injected doublets)
E2c  ambient carryover vs spike fraction                  (injected ambient)
E3   error rate vs genome similarity

Two corrections relative to naive scoring, both material:

1. Truth is the STRAIN GROUP, not the assembly accession. Genomes collapsed by
   05_strain_groups.py are the same strain deposited twice; links between them
   are correct, and counting them as errors puts an unreachable floor under FDR.

2. Every false link is classified before it counts. A link to a genome carrying
   a 100%-identical plasmid is a genuine biological ambiguity, not an
   algorithmic mistake. Strict and adjusted precision are both reported;
   adjusted is the number M5.1 should be assessed against, with the ambiguous
   fraction stated alongside it.

Usage:
    python 04_analyse.py --tiers 0.01x 0.1x 0.2x 0.4x 0.8x 1x \
                         --linkdir links --truth truth --out results
"""
import argparse
import gzip
import os
from collections import defaultdict

import numpy as np
from scipy.stats import spearmanr


def maybe(path):
    return path if os.path.exists(path) else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tiers", nargs="+", required=True)
    ap.add_argument("--linkdir", default="links")
    ap.add_argument("--truth", default="truth")
    ap.add_argument("--out", default="results")
    ap.add_argument("--near-clone-ani", type=float, default=99.9,
                    help="below the strain-group collapse threshold: pairs in "
                         "this band are separable but hard, so links between "
                         "them are ambiguous rather than plainly wrong")
    ap.add_argument("--doublet-spec", type=float, default=0.99,
                    help="clean-barcode specificity at which doublet "
                         "sensitivity is reported")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    T = lambda f: os.path.join(args.truth, f)

    # ------------------------------------------------------------ truth ----
    g2sg = {}
    if maybe(T("strain_groups.tsv")):
        with open(T("strain_groups.tsv")) as fh:
            next(fh)
            for line in fh:
                g, sg, _n = line.rstrip("\n").split("\t")
                g2sg[g] = sg
    else:
        print("WARNING: no strain_groups.tsv; scoring on raw accessions. "
              "Indistinguishable genome pairs will inflate FDR.")
    sg = lambda g: g2sg.get(g, g)

    bc_tier, bc_genome = {}, {}
    with open(T("barcode_truth.tsv")) as fh:
        for line in fh:
            tier, bc, g = line.rstrip("\n").split("\t")
            bc_tier[bc] = tier
            bc_genome[bc] = g

    lr_genome = {}
    with open(T("longread_truth.tsv")) as fh:
        for line in fh:
            lr, g = line.rstrip("\n").split("\t")
            lr_genome[lr] = g

    ani = defaultdict(dict)
    with open(T("genome_pairs.tsv")) as fh:
        next(fh)
        for line in fh:
            g1, g2, a = line.rstrip("\n").split("\t")
            a = float(a)
            ani[g1][g2] = a
            ani[g2][g1] = a

    nn_ani = {}
    nn_file = maybe(T("strain_nn.tsv")) or maybe(T("genome_nn.tsv"))
    if nn_file:
        with open(nn_file) as fh:
            next(fh)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                nn_ani[f[0]] = float(f[2])

    shared_plasmid = defaultdict(set)
    if maybe(T("plasmid_pairs.tsv")):
        with open(T("plasmid_pairs.tsv")) as fh:
            next(fh)
            for line in fh:
                _p1, g1, _p2, g2, _a = line.rstrip("\n").split("\t")
                shared_plasmid[g1].add(g2)
                shared_plasmid[g2].add(g1)

    def classify(gb, gl):
        """gb = barcode's genome, gl = long read's genome."""
        if sg(gb) == sg(gl):
            return "true_link"
        if gl in shared_plasmid.get(gb, ()):
            return "shared_plasmid"
        if ani.get(gb, {}).get(gl, 0.0) >= args.near_clone_ani:
            return "near_clone"
        return "true_fp"

    # -------------------------------------------------------- per tier ----
    e1_rows, e3_rows, dbl_rows, amb_rows = [], [], [], []

    for tier in args.tiers:
        path = os.path.join(args.linkdir, tier, "retained.tsv.gz")
        if not os.path.exists(path):
            print(f"  skip {tier}: {path} missing")
            continue

        inj = {}
        ipath = T(f"injected_{tier}.tsv")
        if os.path.exists(ipath):
            with open(ipath) as fh:
                next(fh)
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    inj[f[0]] = dict(kind=f[1], major=f[4],
                                     minor=f[5].split(","), stratum=f[6],
                                     ani=float(f[7]), frac=float(f[8]))

        counts = defaultdict(int)
        per_genome = defaultdict(lambda: [0, 0])
        per_bc_kept = defaultdict(int)
        candidates = defaultdict(int)
        bc_sg_hits = defaultdict(lambda: defaultdict(int))

        with gzip.open(path, "rt") as fh:
            next(fh)
            for line in fh:
                bc, lr, _h, _k, _p, _q, keep = line.rstrip("\n").split("\t")
                gl = lr_genome.get(lr)
                if gl is None:
                    continue
                gb = inj[bc]["major"] if bc in inj else bc_genome.get(bc)
                if gb is None:
                    continue
                candidates[bc] += 1
                if keep != "1":
                    continue
                per_bc_kept[bc] += 1
                bc_sg_hits[bc][sg(gl)] += 1
                if bc in inj:
                    continue                      # negatives scored separately
                cls = classify(gb, gl)
                counts[cls] += 1
                if cls == "true_link":
                    per_genome[gb][0] += 1
                elif cls == "true_fp":
                    per_genome[gb][1] += 1

        clean_cand = sum(v for b, v in candidates.items() if b not in inj)
        kept = sum(counts.values())
        tp, fp = counts["true_link"], counts["true_fp"]
        amb = counts["shared_plasmid"] + counts["near_clone"]

        strict = tp / kept if kept else float("nan")
        adjusted = tp / (tp + fp) if (tp + fp) else float("nan")
        recall = tp / clean_cand if clean_cand else float("nan")
        clean_kept = [v for b, v in per_bc_kept.items() if b not in inj]

        e1_rows.append(dict(
            tier=tier, kept=kept, true_link=tp, true_fp=fp,
            shared_plasmid=counts["shared_plasmid"],
            near_clone=counts["near_clone"],
            precision_strict=strict, precision_adjusted=adjusted,
            fdr_adjusted=(1 - adjusted) if adjusted == adjusted else float("nan"),
            recall_reachable=recall,
            mean_links_per_bc=float(np.mean(clean_kept)) if clean_kept else 0.0))

        for g, (t, f) in per_genome.items():
            if t + f and g in nn_ani:
                e3_rows.append((tier, g, nn_ani[g], f / (t + f), t + f))

        # ---- minor-strain fraction: the doublet / ambient statistic -------
        def minor_fraction(bc):
            h = bc_sg_hits.get(bc)
            if not h:
                return None
            tot = sum(h.values())
            return (tot - max(h.values())) / tot if tot else None

        clean_mf = [m for b in per_bc_kept if b not in inj
                    for m in [minor_fraction(b)] if m is not None]
        if clean_mf:
            thr = float(np.quantile(clean_mf, args.doublet_spec))
        else:
            thr = float("nan")

        by_stratum = defaultdict(lambda: [0, 0])
        for b, meta in inj.items():
            if meta["kind"] != "doublet":
                continue
            m = minor_fraction(b)
            if m is None:
                continue
            by_stratum[meta["stratum"]][1] += 1
            if m > thr:
                by_stratum[meta["stratum"]][0] += 1
        for stratum, (det, tot) in sorted(by_stratum.items()):
            dbl_rows.append((tier, stratum, det, tot,
                             det / tot if tot else float("nan"), thr))

        by_frac = defaultdict(lambda: [0.0, 0])
        for b, meta in inj.items():
            if meta["kind"] != "ambient":
                continue
            h = bc_sg_hits.get(b)
            if not h:
                continue
            major_sg = sg(meta["major"])
            tot = sum(h.values())
            carry = (tot - h.get(major_sg, 0)) / tot if tot else 0.0
            by_frac[meta["frac"]][0] += carry
            by_frac[meta["frac"]][1] += 1
        for frac, (s, n) in sorted(by_frac.items()):
            amb_rows.append((tier, frac, n, s / n if n else float("nan")))

    # ------------------------------------------------------------ write ----
    cols = ["tier", "kept", "true_link", "true_fp", "shared_plasmid",
            "near_clone", "precision_strict", "precision_adjusted",
            "fdr_adjusted", "recall_reachable", "mean_links_per_bc"]
    p1 = os.path.join(args.out, "E1_E2_coverage_and_error.tsv")
    with open(p1, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in e1_rows:
            fh.write("\t".join(f"{r[c]:.4f}" if isinstance(r[c], float)
                               else str(r[c]) for c in cols) + "\n")

    p2 = os.path.join(args.out, "E2b_doublet_detection.tsv")
    with open(p2, "w") as fh:
        fh.write("tier\tstratum\tdetected\ttotal\tsensitivity\tthreshold\n")
        for r in dbl_rows:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]:.4f}\t{r[5]:.4f}\n")

    p3 = os.path.join(args.out, "E2c_ambient_carryover.tsv")
    with open(p3, "w") as fh:
        fh.write("tier\tspike_fraction\tn_barcodes\tmean_carryover\n")
        for r in amb_rows:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.4f}\n")

    p4 = os.path.join(args.out, "E3_error_vs_ani.tsv")
    with open(p4, "w") as fh:
        fh.write("tier\tgenome\tnn_ani\tfdr\tn_links\n")
        for tier, g, a, f, n in e3_rows:
            fh.write(f"{tier}\t{g}\t{a:.2f}\t{f:.6f}\t{n}\n")

    # ----------------------------------------------------------- report ----
    print("\n=== E1/E2  precision & error by coverage tier ===")
    print(f"{'tier':>7} {'kept':>11} {'prec_adj':>9} {'FDR':>8} "
          f"{'recall':>8} {'ambig':>9} {'links/bc':>9}")
    for r in e1_rows:
        print(f"{r['tier']:>7} {r['kept']:>11} {r['precision_adjusted']:>9.4f} "
              f"{r['fdr_adjusted']:>8.4f} {r['recall_reachable']:>8.4f} "
              f"{r['shared_plasmid'] + r['near_clone']:>9} "
              f"{r['mean_links_per_bc']:>9.1f}")

    if dbl_rows:
        print(f"\n=== E2b  doublet detection "
              f"(threshold at {args.doublet_spec:.0%} clean specificity) ===")
        for tier, stratum, det, tot, sens, thr in dbl_rows:
            print(f"  {tier:>7}  {stratum:<10} {det:>5}/{tot:<5} "
                  f"sensitivity={sens:.3f}  (minor-frac cutoff {thr:.4f})")

    if amb_rows:
        print("\n=== E2c  ambient carryover into retained links ===")
        for tier, frac, n, mean in amb_rows:
            print(f"  {tier:>7}  spike={frac:<6} n={n:<5} "
                  f"mean carryover={mean:.4f}")

    print("\n=== E3  does error scale with genome similarity? ===")
    bins = [(0, 96), (96, 97), (97, 98), (98, 99),
            (99, 99.5), (99.5, 99.9), (99.9, 100.01)]
    for tier in args.tiers:
        rows = [(a, f) for t, _g, a, f, _n in e3_rows if t == tier]
        if len(rows) < 10:
            continue
        a = np.array([r[0] for r in rows])
        f = np.array([r[1] for r in rows])
        rho, p = spearmanr(a, f)
        print(f"\n  tier {tier}:  Spearman rho={rho:+.3f}  p={p:.2e}  "
              f"(n={len(rows)} genomes)")
        for lo, hi in bins:
            m = (a >= lo) & (a < hi)
            if m.sum():
                print(f"    NN-ANI {lo:5.1f}-{hi:<6.2f} n={m.sum():>4}  "
                      f"mean FDR={f[m].mean():.4f}")

    for p in (p1, p2, p3, p4):
        print(f"\nwrote {p}")


if __name__ == "__main__":
    main()
