#!/usr/bin/env python3
"""
Stage 4: prune spurious barcode-longread links with a hypergeometric test
and Benjamini-Hochberg correction.

Model, per tier:
    N = total short-read alignments retained in the tier
    K = alignments from barcode b
    n = alignments landing on long read r (any barcode)
    k = alignments from b landing on r
    p = P(X >= k | hypergeom(N, K, n))

The null assumes uniform alignment opportunity, which repeats violate. Run
repeat masking (Stage 2) first or rRNA operons alone will inflate the FDR.

Usage:
    python 03_score.py --links links/1x --out links/1x/retained.tsv.gz --fdr 0.01
"""
import argparse
import glob
import gzip
import os
import sys

import numpy as np
from scipy.stats import hypergeom


def open_maybe_gz(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--links", required=True,
                    help="directory of *.links.tsv.gz from 02_link.sh")
    ap.add_argument("--out", required=True)
    ap.add_argument("--fdr", type=float, default=0.01)
    ap.add_argument("--min-reads", type=int, default=2,
                    help="minimum distinct short reads before a pair is tested")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.links, "*.links.tsv.gz")))
    if not files:
        sys.exit(f"no *.links.tsv.gz under {args.links}")

    # ---- pass 1: marginals -------------------------------------------------
    K = {}   # barcode -> total reads
    n_ = {}  # long read -> total reads
    N = 0
    n_pairs = 0
    for path in files:
        with open_maybe_gz(path) as fh:
            for line in fh:
                bc, lr, _hits, reads = line.rstrip("\n").split("\t")
                r = int(reads)
                K[bc] = K.get(bc, 0) + r
                n_[lr] = n_.get(lr, 0) + r
                N += r
                n_pairs += 1
    print(f"[score] barcodes={len(K)} longreads={len(n_)} "
          f"pairs={n_pairs} total_reads={N}", file=sys.stderr)

    # ---- pass 2: p-values --------------------------------------------------
    recs = []
    pvals = []
    for path in files:
        with open_maybe_gz(path) as fh:
            for line in fh:
                bc, lr, hits, reads = line.rstrip("\n").split("\t")
                k = int(reads)
                if k < args.min_reads:
                    continue
                p = hypergeom.sf(k - 1, N, K[bc], n_[lr])
                recs.append((bc, lr, int(hits), k))
                pvals.append(p)

    if not recs:
        sys.exit("no pairs passed --min-reads; lower it or check the linkage step")

    pvals = np.asarray(pvals, dtype=float)
    m = pvals.size
    print(f"[score] tested {m} pairs", file=sys.stderr)

    # ---- Benjamini-Hochberg ------------------------------------------------
    order = np.argsort(pvals)
    ranked = pvals[order]
    qvals = np.empty(m, dtype=float)
    qvals[order] = np.minimum.accumulate(
        (ranked * m / np.arange(1, m + 1))[::-1])[::-1]
    np.clip(qvals, 0, 1, out=qvals)

    keep = qvals <= args.fdr
    print(f"[score] retained {keep.sum()} / {m} at q<={args.fdr} "
          f"({100*keep.mean():.1f}%)", file=sys.stderr)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with gzip.open(args.out, "wt") as fh:
        fh.write("barcode\tlong_read\tn_hits\tn_reads\tpval\tqval\tretained\n")
        for (bc, lr, h, k), p, q, kp in zip(recs, pvals, qvals, keep):
            fh.write(f"{bc}\t{lr}\t{h}\t{k}\t{p:.6g}\t{q:.6g}\t{int(kp)}\n")
    print(f"[score] wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
