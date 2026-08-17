#!/usr/bin/env python3
"""
Resolve each long read to its source CONTIG, and thus to plasmid vs chromosome.

pbsim was run per contig, and the middle field of the long-read header is the
1-based contig index in the order contigs appear in the genome's .fna:

    GCA_000167875.2_ASM16787v2_genomic_0003_142
    |__________ genome ______________| |__| |_|
                                    contig  read

Verified on this dataset: for all 473 genomes, max(chunk index) == number of
contigs in the .fna, with no genome exceeding its contig count.

This is what makes M5.1 measurable. Without it, long reads resolve only to a
genome and there is no way to tell a plasmid read from a chromosomal one, so
plasmid-to-host linkage accuracy cannot be computed at all.

Contig order is taken from truth/contig_truth.tsv, which 01_build_truth.sh
writes in .fna order per genome. Do not sort that file.

Output
------
    longread_contig_truth.tsv   long_read, genome, contig, type, contig_len

Usage
-----
    python 07_longread_contigs.py --truth truth
"""
import argparse
import os
import re
import sys
from collections import defaultdict

HEADER = re.compile(r"^(?P<genome>.+)_(?P<chunk>\d+)_(?P<read>\d+)$")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", default="truth")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    T = lambda f: os.path.join(args.truth, f)
    out_path = args.out or T("longread_contig_truth.tsv")

    # contig order per genome, as written by 01_build_truth.sh
    order = defaultdict(list)
    with open(T("contig_truth.tsv")) as fh:
        for line in fh:
            acc, genome, ctype, clen = line.rstrip("\n").split("\t")
            order[genome].append((acc, ctype, int(clen)))
    print(f"contigs loaded for {len(order)} genomes")

    n = ok = bad_parse = bad_index = 0
    type_counts = defaultdict(int)
    unresolved = []

    with open(T("longread_truth.tsv")) as fh, open(out_path, "w") as out:
        out.write("long_read\tgenome\tcontig\ttype\tcontig_len\n")
        for line in fh:
            lr, genome = line.rstrip("\n").split("\t")
            n += 1
            m = HEADER.match(lr)
            if not m:
                bad_parse += 1
                if len(unresolved) < 5:
                    unresolved.append(lr)
                continue
            idx = int(m.group("chunk"))
            contigs = order.get(genome)
            if not contigs or not (1 <= idx <= len(contigs)):
                bad_index += 1
                if len(unresolved) < 5:
                    unresolved.append(lr)
                continue
            acc, ctype, clen = contigs[idx - 1]
            out.write(f"{lr}\t{genome}\t{acc}\t{ctype}\t{clen}\n")
            type_counts[ctype] += 1
            ok += 1

    print(f"long reads          {n}")
    print(f"resolved            {ok} ({100*ok/n:.2f}%)")
    for t, c in sorted(type_counts.items()):
        print(f"    {t:<12} {c:>10} ({100*c/ok:.2f}%)")
    if bad_parse:
        print(f"unparseable headers {bad_parse}")
    if bad_index:
        print(f"chunk out of range  {bad_index}")
    if unresolved:
        print("examples:")
        for u in unresolved:
            print(f"    {u}")
    if ok < n:
        print("\n  NOTE: unresolved reads are dropped from plasmid scoring.",
              file=sys.stderr)
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()
