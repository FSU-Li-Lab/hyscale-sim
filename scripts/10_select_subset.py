#!/usr/bin/env python3
"""
Select a subset of genomes whose ANI structure resembles a real microbiome.

Why subset at all
-----------------
473 Enterobacterales at 96-99% ANI, overwhelmingly E. coli, means a core-genome
150 bp read has ~473 x 13 = 6,000 valid targets in a 13x long-read reference.
minimap2's -N cap then keeps a near-random 50 of them and own-genome hits are
crowded out. That is ~100x harder than the fecal community Task 5 targets,
where a handful of E. coli strains compete against hundreds of distant taxa and
a read has ~65 targets.

The fix is not fewer genomes, it is the right ANI SHAPE:

    a small CLOSE cluster   (>= --close-min-ani)  -- the actual strain problem
    a DIVERSE background    (<= --max-background-ani between any two)

Reads from the close cluster stay hard, which is the point. Reads from
everything else become discriminative, which is what makes the run tractable
and what a real sample looks like.

Selection priorities, in order:
  1. a near-clone pair below the strain-group collapse threshold (phasing test)
  2. hosts of cross-host plasmids whose CHROMOSOMES are distant -- the ideal
     M5.1 case: plasmid identical, host unambiguous, so a mis-assignment is a
     real error rather than an ambiguity
  3. the close cluster, grown from the densest high-ANI neighbourhood
  4. diverse background, greedily maximising distance

Usage
-----
    python 10_select_subset.py --truth truth --out truth/subset_genomes.txt \
        --n-close 5 --n-diverse 25
"""
import argparse
import os
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", default="truth")
    ap.add_argument("--out", default="truth/subset_genomes.txt")
    ap.add_argument("--n-close", type=int, default=5,
                    help="size of the closely-related strain cluster")
    ap.add_argument("--close-min-ani", type=float, default=98.0)
    ap.add_argument("--n-diverse", type=int, default=25,
                    help="background genomes, mutually distant")
    ap.add_argument("--max-background-ani", type=float, default=92.0,
                    help="no two background genomes may exceed this, and none "
                         "may exceed it against the close cluster")
    ap.add_argument("--near-clone-ani", type=float, default=99.9)
    ap.add_argument("--collapse-ani", type=float, default=99.99,
                    help="pairs at or above this are one strain group and "
                         "test nothing")
    args = ap.parse_args()

    T = lambda f: os.path.join(args.truth, f)

    # ---------------------------------------------------------------- load
    ani = defaultdict(dict)
    with open(T("genome_pairs.tsv")) as fh:
        next(fh)
        for line in fh:
            g1, g2, a = line.rstrip("\n").split("\t")
            a = float(a)
            ani[g1][g2] = a
            ani[g2][g1] = a

    genomes = set()
    with open(T("contig_truth.tsv")) as fh:
        for line in fh:
            genomes.add(line.split("\t")[1])
    genomes = sorted(genomes)

    sg = {}
    if os.path.exists(T("strain_groups.tsv")):
        with open(T("strain_groups.tsv")) as fh:
            next(fh)
            for line in fh:
                g, grp, _n = line.rstrip("\n").split("\t")
                sg[g] = grp

    plasmid_partners = defaultdict(set)
    if os.path.exists(T("plasmid_pairs.tsv")):
        with open(T("plasmid_pairs.tsv")) as fh:
            next(fh)
            for line in fh:
                _p1, g1, _p2, g2, _a = line.rstrip("\n").split("\t")
                plasmid_partners[g1].add(g2)
                plasmid_partners[g2].add(g1)

    A = lambda x, y: ani.get(x, {}).get(y, 0.0)

    selected, why = [], {}

    def add(g, reason):
        if g in genomes and g not in selected:
            selected.append(g)
            why[g] = reason
            return True
        return False

    # ------------------------------------------- 1. near-clone pair --------
    clones = sorted(
        ((A(x, y), x, y) for i, x in enumerate(genomes) for y in genomes[i+1:]
         if args.near_clone_ani <= A(x, y) < args.collapse_ani),
        reverse=True)
    if clones:
        a, x, y = clones[0]
        add(x, f"near-clone {a:.2f}%")
        add(y, f"near-clone {a:.2f}%")

    # ------------------------- 2. distant hosts sharing a plasmid ----------
    gold = []
    for g1, partners in plasmid_partners.items():
        for g2 in partners:
            if g1 < g2 and A(g1, g2) < args.max_background_ani:
                gold.append((A(g1, g2), g1, g2))
    gold.sort()
    for a, g1, g2 in gold[:3]:
        add(g1, f"shared plasmid, host ANI {a:.1f}%")
        add(g2, f"shared plasmid, host ANI {a:.1f}%")

    # ------------------------------------------- 3. close cluster ----------
    close_nb = {g: [h for h in genomes
                    if h != g and A(g, h) >= args.close_min_ani
                    and A(g, h) < args.collapse_ani]
                for g in genomes}
    seeds = sorted(genomes, key=lambda g: -len(close_nb[g]))
    cluster = []
    for seed in seeds:
        if not close_nb[seed]:
            continue
        cand = [seed]
        for h in sorted(close_nb[seed], key=lambda h: -A(seed, h)):
            if all(A(h, c) >= args.close_min_ani for c in cand):
                cand.append(h)
            if len(cand) >= args.n_close:
                break
        if len(cand) >= min(args.n_close, 3):
            cluster = cand
            break
    for g in cluster:
        add(g, f"close cluster (>={args.close_min_ani}%)")

    anchor = selected[:]                     # everything chosen so far

    # ------------------------------------------- 4. diverse background -----
    def max_ani_to_selected(g):
        return max((A(g, s) for s in selected), default=0.0)

    pool = [g for g in genomes if g not in selected]
    # prefer plasmid hosts, then genomes far from everything already picked
    pool.sort(key=lambda g: (-len(plasmid_partners.get(g, ())),
                             max_ani_to_selected(g)))
    added = 0
    for g in pool:
        if added >= args.n_diverse:
            break
        if max_ani_to_selected(g) < args.max_background_ani:
            tag = ("background + plasmid host"
                   if plasmid_partners.get(g) else "background")
            if add(g, f"{tag} (max ANI {max_ani_to_selected(g):.1f}%)"):
                added += 1

    with open(args.out, "w") as fh:
        fh.write("\n".join(selected) + "\n")

    # ---------------------------------------------------------------- report
    print(f"selected {len(selected)} genomes -> {args.out}\n")
    for g in selected:
        print(f"  {g:<52} {why[g]}")

    print("\n--- ANI structure of the subset ---")
    hi = [(A(x, y), x, y) for i, x in enumerate(selected)
          for y in selected[i+1:] if A(x, y) > 0]
    hi.sort(reverse=True)
    bands = [(99.9, 100.1), (99, 99.9), (98, 99), (95, 98),
             (90, 95), (0.01, 90)]
    for lo, up in bands:
        n = sum(1 for a, _x, _y in hi if lo <= a < up)
        print(f"  pairs {lo:5.1f}-{up:<5.1f}  {n}")
    print(f"  pairs below skani screen (very distant): "
          f"{len(selected)*(len(selected)-1)//2 - len(hi)}")

    print("\n--- expected mapping difficulty ---")
    print("  confusable = genomes at >=95% ANI (a 150bp core read cannot")
    print("  distinguish these); targets ~ confusable x 13x long-read depth")
    worst = 0
    for g in selected:
        c = 1 + sum(1 for h in selected if h != g and A(g, h) >= 95.0)
        worst = max(worst, c)
        if c > 1:
            print(f"  {g:<52} confusable={c}  ~{c*13} targets/read")
    print(f"  worst case ~{worst*13} targets/read "
          f"(was ~6000 with all 473; -N 100 covers this)")

    print("\n--- M5.1 test cases inside the subset ---")
    n = 0
    for g1 in selected:
        for g2 in plasmid_partners.get(g1, ()):
            if g2 in selected and g1 < g2:
                n += 1
                if n <= 8:
                    print(f"  {g1[:34]:<36} <-> {g2[:34]:<36} "
                          f"host ANI {A(g1,g2):.1f}%")
    print(f"  total cross-host plasmid pairs retained: {n}")
    if n == 0:
        print("  !! none -- M5.1 cannot be scored on this subset")

    print("\n--- strain groups represented ---")
    grps = {sg.get(g, g) for g in selected}
    print(f"  {len(grps)} groups across {len(selected)} genomes")


if __name__ == "__main__":
    main()
