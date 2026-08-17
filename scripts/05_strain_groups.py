#!/usr/bin/env python3
"""
Collapse indistinguishable genomes into strain groups.

Genomes at ~100% ANI are separate assembly accessions for what is effectively
the same strain. Simulated separately, they produce barcodes and long reads
that no method can tell apart, so every cross-link between them scores as a
false positive that is impossible to avoid. Left uncorrected this puts a hard
floor under FDR that has nothing to do with the pipeline.

Each connected component at --threshold becomes one strain for truth purposes.
A link between two members of a group is a true link, not an error.

Choosing --threshold:
    100.00  collapse only what skani cannot separate at all
    99.99   collapse pairs differing by <~500 SNPs over 5 Mb (recommended)
    99.9    aggressive; also collapses genuinely phaseable pairs, i.e. deletes
            the cases the allelic phasing layer exists to solve

Outputs
-------
    strain_groups.tsv     genome, strain_group, group_size
    strain_group_sizes.tsv strain_group, n_genomes, members
    strain_nn.tsv         genome, nn_genome, nn_ani, nn_strain_group
                          nearest neighbour OUTSIDE the genome's own group --
                          this is the covariate E3 regresses error against,
                          because within-group ANI is uninformative once
                          collapsed

Usage
-----
    python 05_strain_groups.py --pairs truth/genome_pairs.tsv \
        --all-genomes truth/genome_list.txt --outdir truth --threshold 99.99
"""
import argparse
import os
from collections import defaultdict


class UnionFind:
    def __init__(self):
        self.parent = {}

    def add(self, x):
        self.parent.setdefault(x, x)

    def find(self, x):
        self.add(x)
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != root:      # path compression
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default="truth/genome_pairs.tsv")
    ap.add_argument("--all-genomes", default=None,
                    help="list of every genome, so singletons absent from the "
                         "pair table still get a group")
    ap.add_argument("--outdir", default="truth")
    ap.add_argument("--threshold", type=float, default=99.99)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    uf = UnionFind()
    all_g = set()
    ani = defaultdict(dict)

    n_merge = 0
    with open(args.pairs) as fh:
        next(fh)
        for line in fh:
            g1, g2, a = line.rstrip("\n").split("\t")
            a = float(a)
            all_g.update((g1, g2))
            uf.add(g1)
            uf.add(g2)
            ani[g1][g2] = a
            ani[g2][g1] = a
            if a >= args.threshold:
                uf.union(g1, g2)
                n_merge += 1

    if args.all_genomes and os.path.exists(args.all_genomes):
        with open(args.all_genomes) as fh:
            for line in fh:
                g = line.strip()
                if g:
                    all_g.add(g)
                    uf.add(g)

    groups = defaultdict(list)
    for g in sorted(all_g):
        groups[uf.find(g)].append(g)

    ordered = sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    name = {root: f"SG{i+1:04d}" for i, (root, _) in enumerate(ordered)}
    g2sg = {g: name[root] for root, members in ordered for g in members}

    with open(os.path.join(args.outdir, "strain_groups.tsv"), "w") as fh:
        fh.write("genome\tstrain_group\tgroup_size\n")
        for root, members in ordered:
            for g in members:
                fh.write(f"{g}\t{name[root]}\t{len(members)}\n")

    with open(os.path.join(args.outdir, "strain_group_sizes.tsv"), "w") as fh:
        fh.write("strain_group\tn_genomes\tmembers\n")
        for root, members in ordered:
            fh.write(f"{name[root]}\t{len(members)}\t{','.join(members)}\n")

    # cross-group nearest neighbour
    with open(os.path.join(args.outdir, "strain_nn.tsv"), "w") as fh:
        fh.write("genome\tnn_genome\tnn_ani\tnn_strain_group\n")
        for g in sorted(all_g):
            best, best_a = None, 0.0
            for h, a in ani.get(g, {}).items():
                if g2sg.get(h) == g2sg.get(g):
                    continue
                if a > best_a:
                    best, best_a = h, a
            fh.write(f"{g}\t{best or 'NA'}\t{best_a:.2f}\t"
                     f"{g2sg.get(best, 'NA')}\n")

    multi = [(name[r], m) for r, m in ordered if len(m) > 1]
    print(f"genomes                  {len(all_g)}")
    print(f"strain groups            {len(ordered)}")
    print(f"collapsed groups (>1)    {len(multi)}")
    print(f"genomes inside them      {sum(len(m) for _, m in multi)}")
    print(f"merge events at >={args.threshold}  {n_merge}")
    if multi:
        print("\ncollapsed groups:")
        for gname, members in multi:
            print(f"  {gname}  n={len(members)}")
            for g in members:
                print(f"      {g}")
    print(f"\nwrote {args.outdir}/strain_groups.tsv")
    print(f"wrote {args.outdir}/strain_nn.tsv")


if __name__ == "__main__":
    main()
