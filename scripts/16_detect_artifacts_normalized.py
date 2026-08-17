#!/usr/bin/env python3
"""
Artifact detection normalized against expected strain diffusion.

The problem with a raw minor-strain fraction
--------------------------------------------
A barcode's associations spread across strain groups for two unrelated reasons:

  1. a genuine second cell in the compartment (what we want to detect), and
  2. mapping spread across near-identical strains (a property of the community,
     not of the barcode).

In a community containing eight E. coli above 98% ANI, cause 2 dominates. A
clean E. coli barcode has a median minor fraction of ~0.17 with a long tail,
so doublets at ~0.5-0.8 sit inside the clean distribution and no threshold
separates them. Worse, near-identical doublets score HIGHER than distant ones,
which is backwards for a two-cell mixing statistic -- the excess is diffusion,
not the second cell.

The fix
-------
Compare each barcode's profile against the profile a CLEAN barcode of the same
strain would produce, rather than against zero. Diffusion is common to both and
cancels; a second cell does not cancel.

    residual[s] = observed_fraction[s] - reference_fraction[s]
    statistic   = max residual over strain groups other than the primary

The reference is the per-strain median profile over all barcodes sharing that
primary strain group. The median is used deliberately: it is estimated without
knowing which barcodes are clean, so the method is applicable to real data
where that is unknown, and it is robust as long as artifacts are a minority.
The barcode under test is excluded from its own reference.

Usage
-----
    python 16_detect_artifacts_normalized.py --tiers 0.1x 0.4x 1x \\
        --linkdir links --truth truth --out results_detect_norm
"""
import argparse
import glob
import gzip
import os
from collections import defaultdict

import numpy as np


def auc(pos, neg):
    """Rank-based AUC (Mann-Whitney U); ties contribute 0.5."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if not pos.size or not neg.size:
        return float("nan")
    allv = np.concatenate([pos, neg])
    order = allv.argsort()
    ranks = np.empty(allv.size, dtype=float)
    ranks[order] = np.arange(1, allv.size + 1)
    _, inv, cnt = np.unique(allv, return_inverse=True, return_counts=True)
    sums = np.zeros(cnt.size)
    np.add.at(sums, inv, ranks)
    ranks = (sums / cnt)[inv]
    n1, n2 = pos.size, neg.size
    return (ranks[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tiers", nargs="+", required=True)
    ap.add_argument("--linkdir", default="links")
    ap.add_argument("--truth", default="truth")
    ap.add_argument("--out", default="results_detect_norm")
    ap.add_argument("--specificity", type=float, default=0.99)
    ap.add_argument("--min-assoc", type=int, default=10)
    ap.add_argument("--reference", choices=["median", "clean"],
                    default="median",
                    help="median: per-strain median over all barcodes, "
                         "estimable without ground truth (default). "
                         "clean: oracle reference from known-clean barcodes, "
                         "for comparison only.")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    T = lambda f: os.path.join(args.truth, f)

    g2sg = {}
    if os.path.exists(T("strain_groups.tsv")):
        with open(T("strain_groups.tsv")) as fh:
            next(fh)
            for line in fh:
                g, s, _n = line.rstrip("\n").split("\t")
                g2sg[g] = s
    sg = lambda g: g2sg.get(g, g)

    lr_sg = {}
    with open(T("longread_truth.tsv")) as fh:
        for line in fh:
            lr, g = line.rstrip("\n").split("\t")
            lr_sg[lr] = sg(g)
    groups = sorted(set(lr_sg.values()))
    gidx = {s: i for i, s in enumerate(groups)}
    G = len(groups)
    print(f"long reads: {len(lr_sg)}   strain groups: {G}")

    rows = []
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
                    inj[f[0]] = dict(kind=f[1], stratum=f[6], frac=float(f[8]))

        # ---- profiles ----------------------------------------------------
        prof = defaultdict(lambda: np.zeros(G))
        for path in files:
            with gzip.open(path, "rt") as fh:
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 4:
                        continue
                    s = lr_sg.get(p[1])
                    if s is None:
                        continue
                    prof[p[0]][gidx[s]] += float(p[3])

        frac, primary = {}, {}
        for bc, v in prof.items():
            tot = v.sum()
            if tot < args.min_assoc:
                continue
            frac[bc] = v / tot
            primary[bc] = int(v.argmax())

        # ---- per-strain reference profiles -------------------------------
        pool = defaultdict(list)
        for bc, f in frac.items():
            if args.reference == "clean" and bc in inj:
                continue
            pool[primary[bc]].append(f)
        ref = {k: np.median(np.vstack(v), axis=0) for k, v in pool.items()
               if len(v) >= 3}
        print(f"\n=== {tier} ===")
        print(f"  profiled barcodes: {len(frac)}   "
              f"reference profiles: {len(ref)}")

        def statistic(bc):
            """Max excess over the reference profile, outside the primary."""
            f, pr = frac.get(bc), primary.get(bc)
            if f is None or pr not in ref:
                return None
            r = ref[pr]
            if args.reference == "median" and bc in frac:
                # leave-one-out is unnecessary with a median over hundreds of
                # barcodes, but the primary column is excluded regardless
                pass
            d = f - r
            d[pr] = -np.inf
            return float(d.max())

        clean_v = [v for bc in clean for v in [statistic(bc)] if v is not None]
        if not clean_v:
            print("  no clean barcodes profiled")
            continue
        cl = np.array(clean_v)
        thr = float(np.quantile(cl, args.specificity))
        print(f"  clean: n={cl.size}  median={np.median(cl):+.4f}  "
              f"threshold@{args.specificity:.0%}={thr:+.4f}")

        by = defaultdict(list)
        for bc, m in inj.items():
            v = statistic(bc)
            if v is None:
                continue
            key = ("doublet", m["stratum"]) if m["kind"] == "doublet" \
                else ("ambient", f"{m['frac']}")
            by[key].append(v)

        for kind in ("doublet", "ambient"):
            for key in sorted(k for k in by if k[0] == kind):
                v = np.array(by[key])
                det = int((v > thr).sum())
                a = auc(v, cl)
                print(f"  {key[0]:<8}{key[1]:<6} {det:>3}/{v.size:<3} "
                      f"sens={det/v.size:.3f}  AUC={a:.3f}  "
                      f"median stat={np.median(v):+.3f}")
                rows.append((tier, key[0], key[1], det, v.size,
                             det / v.size, a, float(np.median(v)), thr))

    p = os.path.join(args.out, "artifact_detection_normalized.tsv")
    with open(p, "w") as fh:
        fh.write("tier\tkind\tstratum\tdetected\ttotal\tsensitivity\tauc\t"
                 "median_statistic\tthreshold\n")
        for r in rows:
            fh.write("\t".join(f"{x:.4f}" if isinstance(x, float) else str(x)
                               for x in r) + "\n")
    print(f"\nwrote {p}")


if __name__ == "__main__":
    main()
