#!/usr/bin/env python3
"""
M5.1: plasmid-to-host linkage accuracy and FDR.

For every long read originating from a plasmid contig, the pipeline's implied
host is the strain group contributing the most retained barcode links. That
assignment is compared against the read's true source genome.

Three outcomes, and keeping them apart is the whole point:

    correct     assigned host == true source strain group
    ambiguous   assigned host carries a plasmid >= --plasmid-ani identical to
                the source plasmid. The read really did come from another cell,
                so this is biologically wrong -- but no method operating on
                sequence can know that, because the molecules are identical.
                Reported separately, never silently folded into either bucket.
    wrong       neither: a genuine linkage error

    accuracy_strict   correct / all assigned
    accuracy_adjusted correct / (correct + wrong)
    FDR               wrong / (correct + wrong)

M5.1 should be assessed on accuracy_adjusted and FDR, with the ambiguous count
quoted alongside. Chromosomal reads are scored identically as a control: if
chromosome accuracy is not clearly higher than plasmid accuracy, the problem is
in the linkage stage rather than in plasmid biology.

Also emits an empirical repeat-prone read list. With no pbsim .maf files there
are no true read coordinates, so repeats cannot be masked by annotation. Long
reads drawing retained links from an unusually large number of distinct strain
groups are the observable signature of rRNA operons, IS elements and other
conserved sequence. Feed the list back as an exclusion set and compare FDR with
and without it.

Usage
-----
    python 08_plasmid_linkage.py --tier 1x --links links/1x/retained.tsv.gz \
        --truth truth --out results
"""
import argparse
import gzip
import os
from collections import defaultdict

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", required=True)
    ap.add_argument("--links", required=True)
    ap.add_argument("--truth", default="truth")
    ap.add_argument("--out", default="results")
    ap.add_argument("--repeat-pct", type=float, default=99.0,
                    help="percentile of distinct-strain-group count above "
                         "which a long read is flagged repeat-prone")
    ap.add_argument("--min-links", type=int, default=1,
                    help="minimum retained links before a long read is assigned")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    T = lambda f: os.path.join(args.truth, f)

    # ---------------------------------------------------------------- truth
    g2sg = {}
    if os.path.exists(T("strain_groups.tsv")):
        with open(T("strain_groups.tsv")) as fh:
            next(fh)
            for line in fh:
                g, sg_, _n = line.rstrip("\n").split("\t")
                g2sg[g] = sg_
    sg = lambda g: g2sg.get(g, g)

    lr_info = {}
    p = T("longread_contig_truth.tsv")
    if not os.path.exists(p):
        raise SystemExit(f"missing {p} -- run 07_longread_contigs.py first")
    with open(p) as fh:
        next(fh)
        for line in fh:
            lr, genome, contig, ctype, _clen = line.rstrip("\n").split("\t")
            lr_info[lr] = (genome, contig, ctype)
    print(f"long reads with contig truth: {len(lr_info)}")

    bc_genome = {}
    with open(T("barcode_truth.tsv")) as fh:
        for line in fh:
            tier, bc, g = line.rstrip("\n").split("\t")
            if tier == args.tier:
                bc_genome[bc] = g
    ipath = T(f"injected_{args.tier}.tsv")
    injected = set()
    if os.path.exists(ipath):
        with open(ipath) as fh:
            next(fh)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                bc_genome[f[0]] = f[4]
                injected.add(f[0])

    # plasmid contig -> genomes carrying a near-identical copy
    plasmid_hosts = defaultdict(set)
    if os.path.exists(T("plasmid_pairs.tsv")):
        with open(T("plasmid_pairs.tsv")) as fh:
            next(fh)
            for line in fh:
                p1, g1, p2, g2, _a = line.rstrip("\n").split("\t")
                plasmid_hosts[p1].add(g2)
                plasmid_hosts[p2].add(g1)

    # ---------------------------------------------------------------- links
    votes = defaultdict(lambda: defaultdict(int))
    with gzip.open(args.links, "rt") as fh:
        next(fh)
        for line in fh:
            bc, lr, _h, k, _pv, _q, keep = line.rstrip("\n").split("\t")
            if keep != "1":
                continue
            gb = bc_genome.get(bc)
            if gb is None or lr not in lr_info:
                continue
            votes[lr][sg(gb)] += int(k)

    if not votes:
        raise SystemExit("no retained links matched the contig truth table")

    # ------------------------------------------------- repeat-prone reads
    breadth = {lr: len(v) for lr, v in votes.items()}
    cutoff = float(np.percentile(list(breadth.values()), args.repeat_pct))
    repeat_reads = {lr for lr, b in breadth.items() if b > cutoff}
    rpath = os.path.join(args.out, f"repeat_prone_{args.tier}.txt")
    with open(rpath, "w") as fh:
        for lr in sorted(repeat_reads):
            fh.write(f"{lr}\t{breadth[lr]}\t{lr_info[lr][2]}\n")

    # ----------------------------------------------------------- scoring
    def score(exclude_repeats):
        res = defaultdict(lambda: defaultdict(int))
        for lr, v in votes.items():
            if exclude_repeats and lr in repeat_reads:
                continue
            if sum(v.values()) < args.min_links:
                continue
            genome, contig, ctype = lr_info[lr]
            assigned = max(v, key=v.get)
            true_sg = sg(genome)
            if assigned == true_sg:
                out = "correct"
            elif any(sg(h) == assigned for h in plasmid_hosts.get(contig, ())):
                out = "ambiguous"
            else:
                out = "wrong"
            res[ctype][out] += 1
            res[ctype]["assigned"] += 1
        return res

    rows = []
    for label, excl in (("all_reads", False), ("repeats_excluded", True)):
        res = score(excl)
        for ctype in ("plasmid", "chromosome"):
            d = res.get(ctype)
            if not d:
                continue
            c, w, a = d["correct"], d["wrong"], d["ambiguous"]
            denom = c + w
            rows.append(dict(
                tier=args.tier, subset=label, contig_type=ctype,
                assigned=d["assigned"], correct=c, wrong=w, ambiguous=a,
                accuracy_strict=c / d["assigned"] if d["assigned"] else float("nan"),
                accuracy_adjusted=c / denom if denom else float("nan"),
                fdr=w / denom if denom else float("nan")))

    total_plasmid = sum(1 for lr, (_g, _c, t) in lr_info.items() if t == "plasmid")
    assigned_plasmid = sum(r["assigned"] for r in rows
                           if r["subset"] == "all_reads" and r["contig_type"] == "plasmid")

    # ------------------------------------------------------------- output
    cols = ["tier", "subset", "contig_type", "assigned", "correct", "wrong",
            "ambiguous", "accuracy_strict", "accuracy_adjusted", "fdr"]
    opath = os.path.join(args.out, f"M5.1_plasmid_linkage_{args.tier}.tsv")
    with open(opath, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[c]:.4f}" if isinstance(r[c], float)
                               else str(r[c]) for c in cols) + "\n")

    print(f"\n=== M5.1  plasmid-to-host linkage, tier {args.tier} ===")
    print(f"{'subset':<18}{'type':<12}{'assigned':>9}{'acc_adj':>9}"
          f"{'FDR':>8}{'ambig':>8}")
    for r in rows:
        print(f"{r['subset']:<18}{r['contig_type']:<12}{r['assigned']:>9}"
              f"{r['accuracy_adjusted']:>9.4f}{r['fdr']:>8.4f}"
              f"{r['ambiguous']:>8}")

    print(f"\nplasmid long reads in truth   {total_plasmid}")
    print(f"plasmid long reads assigned   {assigned_plasmid}"
          f" ({100*assigned_plasmid/total_plasmid:.1f}%)" if total_plasmid else "")
    print(f"repeat-prone reads flagged    {len(repeat_reads)} "
          f"(>{cutoff:.0f} distinct strain groups)")
    print(f"\nwrote {opath}")
    print(f"wrote {rpath}")


if __name__ == "__main__":
    main()
