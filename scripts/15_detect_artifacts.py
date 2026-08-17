#!/usr/bin/env python3
"""
Artifact detection on the UNPRUNED association profile.

Why this is separate from 04_analyse.py
---------------------------------------
Linkage and artifact detection are different tasks and want different inputs.

Linkage asks whether an individual barcode-read association is real, so it uses
the pruned set: associations that survive hypergeometric testing.

Detection asks whether a barcode group contains one cell or two. That is a
property of the barcode's whole association profile, not of any single link. A
compartment holding two cells produces two peaks across strain groups whether
or not each individual association clears significance. Scoring detection on
the pruned set therefore throws away most of the evidence, and at strict
thresholds leaves injected barcodes with too few associations to profile at all
-- which is exactly what was observed when detection was run on retained.tsv.gz.

This script reads the raw per-batch link tables instead, so detection is scored
on everything minimap2 found.

Statistic
---------
For each barcode, associations are aggregated by the strain group of the long
read, weighted by supporting short reads. The minor fraction is

    minor = 1 - (reads to the largest strain group) / (reads to all groups)

A clean barcode has minor near 0. A doublet approaches 0.5 when its two cells
contribute equally. An ambient-contaminated barcode sits near its spike level.

Thresholds are reported two ways: sensitivity at a threshold set to a chosen
specificity against clean barcodes, and AUC, which is threshold-free and
therefore comparable across coverage tiers.

Usage
-----
    python 15_detect_artifacts.py --tiers 0.01x 0.1x 0.2x 0.4x 0.8x 1x \\
        --linkdir links --truth truth --out results_detect
"""
import argparse
import glob
import gzip
import os
from collections import defaultdict

import numpy as np


def auc(pos, neg):
    """Rank-based AUC (Mann-Whitney U); ties contribute 0.5."""
    if not len(pos) or not len(neg):
        return float("nan")
    allv = np.concatenate([pos, neg])
    order = allv.argsort()
    ranks = np.empty(len(allv), dtype=float)
    ranks[order] = np.arange(1, len(allv) + 1)
    # average ranks within ties
    _, inv, cnt = np.unique(allv, return_inverse=True, return_counts=True)
    sums = np.zeros(len(cnt))
    np.add.at(sums, inv, ranks)
    ranks = (sums / cnt)[inv]
    n1, n2 = len(pos), len(neg)
    return (ranks[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tiers", nargs="+", required=True)
    ap.add_argument("--linkdir", default="links")
    ap.add_argument("--truth", default="truth")
    ap.add_argument("--out", default="results_detect")
    ap.add_argument("--specificity", type=float, default=0.99,
                    help="clean-barcode specificity at which sensitivity is "
                         "reported")
    ap.add_argument("--min-assoc", type=int, default=10,
                    help="barcodes with fewer associations are not profiled")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    T = lambda f: os.path.join(args.truth, f)

    # ---------------------------------------------------------------- truth
    g2sg = {}
    if os.path.exists(T("strain_groups.tsv")):
        with open(T("strain_groups.tsv")) as fh:
            next(fh)
            for line in fh:
                g, sg, _n = line.rstrip("\n").split("\t")
                g2sg[g] = sg
    sg = lambda g: g2sg.get(g, g)

    lr_sg = {}
    with open(T("longread_truth.tsv")) as fh:
        for line in fh:
            lr, g = line.rstrip("\n").split("\t")
            lr_sg[lr] = sg(g)
    print(f"long reads: {len(lr_sg)}")

    rows_out = []
    auc_out = []

    for tier in args.tiers:
        files = sorted(glob.glob(os.path.join(args.linkdir, tier,
                                              "*.links.tsv.gz")))
        if not files:
            print(f"  skip {tier}: no raw link tables")
            continue

        clean = {}
        with open(T("barcode_truth.tsv")) as fh:
            for line in fh:
                t, bc, g = line.rstrip("\n").split("\t")
                if t == tier:
                    clean[bc] = sg(g)

        inj = {}
        ip = T(f"injected_{tier}.tsv")
        if os.path.exists(ip):
            with open(ip) as fh:
                next(fh)
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    inj[f[0]] = dict(kind=f[1], major=f[4], stratum=f[6],
                                     ani=float(f[7]), frac=float(f[8]))

        # ---- aggregate the unpruned profile -----------------------------
        prof = defaultdict(lambda: defaultdict(float))
        for path in files:
            with gzip.open(path, "rt") as fh:
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 4:
                        continue
                    bc, lr, _h, reads = p[0], p[1], p[2], p[3]
                    s = lr_sg.get(lr)
                    if s is None:
                        continue
                    prof[bc][s] += float(reads)

        def minor_fraction(bc):
            d = prof.get(bc)
            if not d:
                return None
            tot = sum(d.values())
            if tot < args.min_assoc:
                return None
            return (tot - max(d.values())) / tot

        clean_mf = [m for bc in clean
                    for m in [minor_fraction(bc)] if m is not None]
        if not clean_mf:
            print(f"  skip {tier}: no clean barcodes profiled")
            continue
        thr = float(np.quantile(clean_mf, args.specificity))

        by_stratum = defaultdict(list)
        amb_by_frac = defaultdict(list)
        for bc, meta in inj.items():
            m = minor_fraction(bc)
            if m is None:
                continue
            if meta["kind"] == "doublet":
                by_stratum[meta["stratum"]].append(m)
            else:
                amb_by_frac[meta["frac"]].append(m)

        print(f"\n=== {tier} ===")
        print(f"  clean barcodes profiled: {len(clean_mf)}  "
              f"median minor fraction {np.median(clean_mf):.4f}")
        print(f"  threshold at {args.specificity:.0%} specificity: {thr:.4f}")

        cl = np.array(clean_mf)
        for stratum in ("far", "mid", "near"):
            vals = np.array(by_stratum.get(stratum, []))
            if not len(vals):
                continue
            det = int((vals > thr).sum())
            a = auc(vals, cl)
            print(f"  doublet {stratum:<5} {det:>3}/{len(vals):<3} "
                  f"sens={det/len(vals):.3f}  AUC={a:.3f}  "
                  f"median minor={np.median(vals):.3f}")
            rows_out.append((tier, "doublet", stratum, det, len(vals),
                             det / len(vals), a, float(np.median(vals)), thr))

        for frac in sorted(amb_by_frac):
            vals = np.array(amb_by_frac[frac])
            det = int((vals > thr).sum())
            a = auc(vals, cl)
            print(f"  ambient {frac:<5} {det:>3}/{len(vals):<3} "
                  f"sens={det/len(vals):.3f}  AUC={a:.3f}  "
                  f"median minor={np.median(vals):.3f}")
            rows_out.append((tier, "ambient", str(frac), det, len(vals),
                             det / len(vals), a, float(np.median(vals)), thr))

        alld = np.concatenate([np.array(v) for v in by_stratum.values()]) \
            if by_stratum else np.array([])
        if len(alld):
            auc_out.append((tier, "doublet_all", auc(alld, cl), len(alld),
                            len(cl)))

    # ---------------------------------------------------------------- write
    p = os.path.join(args.out, "artifact_detection.tsv")
    with open(p, "w") as fh:
        fh.write("tier\tkind\tstratum\tdetected\ttotal\tsensitivity\tauc\t"
                 "median_minor_fraction\tthreshold\n")
        for r in rows_out:
            fh.write("\t".join(str(x) if not isinstance(x, float)
                               else f"{x:.4f}" for x in r) + "\n")
    print(f"\nwrote {p}")

    if auc_out:
        p2 = os.path.join(args.out, "artifact_auc_summary.tsv")
        with open(p2, "w") as fh:
            fh.write("tier\tclass\tauc\tn_positive\tn_negative\n")
            for t, k, a, n1, n2 in auc_out:
                fh.write(f"{t}\t{k}\t{a:.4f}\t{n1}\t{n2}\n")
        print(f"wrote {p2}")


if __name__ == "__main__":
    main()
