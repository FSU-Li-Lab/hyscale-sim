#!/usr/bin/env python3
"""
Parse skani --full-matrix output into long-form pair tables.

skani full-matrix format:
    line 1        : N
    lines 2..N+1  : <name>\t<v1>\t<v2>\t...\t<vN>
Row order == column order. Zero means "below skani's screening threshold",
NOT "0% ANI" -- treated as missing, not as a low value.

Outputs
-------
    genome_pairs.tsv    g1, g2, ani                      (upper triangle)
    genome_nn.tsv       genome, nearest_neighbour, ani   (raw, any partner)
    genome_list.txt     every genome seen
    plasmid_pairs.tsv   p1, g1, p2, g2, ani              (cross-host only)
    pilot_genomes.txt   suggested pilot subset

Pilot selection deliberately avoids pairs at or above --pilot-clone-max
(default 99.995). Those are re-deposits of the same strain: no method can
separate them, so they test nothing. The pilot wants the hardest *separable*
pair instead, plus a plasmid-sharing hub and a distant outgroup.

Usage
-----
    python parse_skani.py --genome-matrix genome_ani.txt \
                          --plasmid-matrix plasmid_ani.txt \
                          --outdir truth
"""
import argparse
import os
import re
import sys


def read_full_matrix(path):
    """Return (names, rows) where rows[i] is a list of floats."""
    names, rows = [], []
    with open(path) as fh:
        first = fh.readline().strip()
        try:
            n_expected = int(first)
        except ValueError:
            sys.exit(f"{path}: first line should be the matrix size, got {first!r}")
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            names.append(parts[0].strip())
            rows.append([float(x) for x in parts[1:] if x.strip() != ""])
    if len(names) != n_expected:
        print(f"WARN {path}: header says {n_expected} rows, read {len(names)}",
              file=sys.stderr)
    bad = [i for i, r in enumerate(rows) if len(r) != len(names)]
    if bad:
        print(f"WARN {path}: {len(bad)} rows have unexpected width "
              f"(first: row {bad[0]} has {len(rows[bad[0]])} of {len(names)})",
              file=sys.stderr)
    return names, rows


def genome_id(name):
    """genome/GCA_000167875.2_ASM16787v2_genomic.fna -> GCA_000167875.2_ASM16787v2_genomic"""
    base = os.path.basename(name)
    return re.sub(r"\.(fna|fa|fasta)(\.gz)?$", "", base)


def plasmid_fields(name):
    """'GCA_x|CP035754.1 E. coli ... plasmid p33' -> ('GCA_x', 'CP035754.1')"""
    if "|" not in name:
        return None, name.split()[0]
    genome, rest = name.split("|", 1)
    return genome.strip(), rest.split()[0].strip()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome-matrix", required=True)
    ap.add_argument("--plasmid-matrix", required=True)
    ap.add_argument("--outdir", default="truth")
    ap.add_argument("--plasmid-ani", type=float, default=99.0,
                    help="min ANI to call two plasmids shared")
    ap.add_argument("--near-clone-ani", type=float, default=99.9)
    ap.add_argument("--pilot-clone-max", type=float, default=99.995,
                    help="pilot avoids pairs at or above this: they are the "
                         "same strain deposited twice and test nothing")
    ap.add_argument("--pilot-size", type=int, default=10)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    out = lambda f: os.path.join(args.outdir, f)

    # ---------------- genomes ----------------
    gnames, grows = read_full_matrix(args.genome_matrix)
    gids = [genome_id(n) for n in gnames]

    with open(out("genome_list.txt"), "w") as fh:
        fh.write("\n".join(gids) + "\n")

    n_pairs = 0
    nn = {}
    undirected = []          # (ani, g1, g2) for i<j
    with open(out("genome_pairs.tsv"), "w") as fh:
        fh.write("g1\tg2\tani\n")
        for i, gi in enumerate(gids):
            row = grows[i]
            best, best_ani = None, 0.0
            for j, gj in enumerate(gids):
                if i == j:
                    continue
                v = row[j]
                if v <= 0:
                    continue
                if v > best_ani:
                    best, best_ani = gj, v
                if i < j:
                    fh.write(f"{gi}\t{gj}\t{v:.2f}\n")
                    undirected.append((v, gi, gj))
                    n_pairs += 1
            nn[gi] = (best, best_ani)

    with open(out("genome_nn.tsv"), "w") as fh:
        fh.write("genome\tnearest_neighbour\tani\n")
        for g in gids:
            p, a = nn[g]
            fh.write(f"{g}\t{p or 'NA'}\t{a:.2f}\n")

    undirected.sort(reverse=True)
    degenerate = [(a, x, y) for a, x, y in undirected if a >= args.pilot_clone_max]
    separable = [(a, x, y) for a, x, y in undirected
                 if args.near_clone_ani <= a < args.pilot_clone_max]

    # ---------------- plasmids ----------------
    pnames, prows = read_full_matrix(args.plasmid_matrix)
    pmeta = [plasmid_fields(n) for n in pnames]

    shared_hosts = {}
    n_cross = 0
    with open(out("plasmid_pairs.tsv"), "w") as fh:
        fh.write("p1\tg1\tp2\tg2\tani\n")
        for i, (g1, p1) in enumerate(pmeta):
            row = prows[i]
            for j in range(i + 1, len(pmeta)):
                v = row[j]
                if v < args.plasmid_ani:
                    continue
                g2, p2 = pmeta[j]
                if g1 == g2:
                    continue
                fh.write(f"{p1}\t{g1}\t{p2}\t{g2}\t{v:.2f}\n")
                n_cross += 1
                shared_hosts.setdefault(g1, set()).add(g2)
                shared_hosts.setdefault(g2, set()).add(g1)

    # ---------------- pilot selection ----------------
    pilot, why = [], {}

    def add(g, reason):
        if g and g not in pilot and len(pilot) < args.pilot_size:
            pilot.append(g)
            why[g] = reason

    if separable:                       # hardest *separable* pair
        a, x, y = separable[0]
        add(x, f"near-clone {a:.2f}%")
        add(y, f"near-clone {a:.2f}%")
    if shared_hosts:                    # plasmid-sharing hub + partners
        hub = max(shared_hosts, key=lambda g: len(shared_hosts[g]))
        add(hub, f"plasmid hub, shares with {len(shared_hosts[hub])}")
        for p in sorted(shared_hosts[hub])[:2]:
            add(p, "plasmid partner of hub")
    if nn:                              # most distant genome as outgroup
        far = min(nn, key=lambda g: nn[g][1] if nn[g][1] > 0 else 999)
        add(far, f"outgroup, NN ANI {nn[far][1]:.2f}%")
    ranked = sorted(nn.items(), key=lambda kv: -kv[1][1])
    for g, (_p, a) in ranked[::max(1, len(ranked) // max(1, args.pilot_size))]:
        add(g, f"ANI-range filler, NN {a:.2f}%")

    with open(out("pilot_genomes.txt"), "w") as fh:
        fh.write("\n".join(pilot) + "\n")

    # ---------------- report ----------------
    print(f"genomes                       {len(gids)}")
    print(f"genome pairs above cutoff     {n_pairs}")
    print(f"pairs >= {args.pilot_clone_max} (indistinguishable)  {len(degenerate)}")
    for a, x, y in degenerate[:15]:
        print(f"    {a:6.2f}  {x}  ==  {y}")
    if degenerate:
        print("    -> collapse these with 05_strain_groups.py before scoring")
    print(f"pairs {args.near_clone_ani}-{args.pilot_clone_max} (hard but separable)  {len(separable)}")
    for a, x, y in separable[:10]:
        print(f"    {a:6.2f}  {x}  <->  {y}")
    print(f"plasmid contigs               {len(pmeta)}")
    print(f"cross-host plasmid pairs >={args.plasmid_ani}  {n_cross}")
    print(f"genomes sharing >=1 plasmid   {len(shared_hosts)}")
    if not n_cross:
        print("\n  !! No cross-host plasmids: M5.1 has no test case as built.")
    print(f"\npilot set ({len(pilot)}) -> {out('pilot_genomes.txt')}")
    for g in pilot:
        print(f"    {g:<55} {why[g]}")


if __name__ == "__main__":
    main()
