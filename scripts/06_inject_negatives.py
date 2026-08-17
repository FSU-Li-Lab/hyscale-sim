#!/usr/bin/env python3
"""
Inject doublets and ambient DNA to create the negative class.

Every barcode in the simulation came from exactly one genome, so FDR has no
false linkages to reject and the dual-layer validation framework cannot be
demonstrated. This builds the missing negatives.

Doublets are stratified by the ANI between the two source genomes. A doublet of
two distant genomes is trivially detectable; a doublet of two near-identical
genomes is the case that defeats barcode clustering. Strata are graded rather
than binary so doublet detection can be reported as a curve against similarity,
which is the same covariate E3 uses.

Default bands (edit with --bands):
    far   ANI < 98
    mid   98 <= ANI < 99
    near  ANI >= 99
Recipients cycle through the bands, so counts come out roughly balanced
instead of collapsing into whichever band happens to be commonest.

Ambient contamination spikes a fraction of reads from several random donor
genomes into an otherwise clean barcode, mimicking free DNA in the emulsion.

Outputs
-------
    injected/<TIER>/{DBL,AMB}xxxxxx{1,2}.fq.gz
    truth/injected_<TIER>.tsv

Usage
-----
    # whole tier
    python 06_inject_negatives.py --root . --tier 1x

    # pilot only, high rate so E2b has enough events to measure
    python 06_inject_negatives.py --root . --tier 1x \
        --genomes truth/pilot_genomes.txt --doublet-rate 0.30 \
        --ambient-rate 0.30 --prefix P
"""
import argparse
import gzip
import os
import random
import shutil
import sys
from collections import defaultdict


def read_fastq(path):
    with gzip.open(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            yield h + fh.readline() + fh.readline() + fh.readline()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--tier", required=True)
    ap.add_argument("--genomes", default=None,
                    help="restrict recipients to genomes listed in this file "
                         "(donors are still drawn from the whole set)")
    ap.add_argument("--doublet-rate", type=float, default=0.05)
    ap.add_argument("--ambient-rate", type=float, default=0.05)
    ap.add_argument("--ambient-fracs", type=float, nargs="+",
                    default=[0.01, 0.05, 0.10])
    ap.add_argument("--bands", type=float, nargs="+", default=[98.0, 99.0],
                    help="ANI edges separating far/mid/near doublet strata")
    ap.add_argument("--prefix", default="",
                    help="extra tag in injected barcode IDs, e.g. P for pilot, "
                         "so pilot and full-run injections do not collide")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    random.seed(args.seed)
    root = args.root
    sc = os.path.join(root, "single_cell")
    outdir = os.path.join(root, "injected", args.tier)
    os.makedirs(outdir, exist_ok=True)

    # ---- load ---------------------------------------------------------------
    bcs, bc_genome = [], {}
    with open(os.path.join(root, "truth", "barcode_truth.tsv")) as fh:
        for line in fh:
            tier, bc, g = line.rstrip("\n").split("\t")
            if tier == args.tier:
                bcs.append(bc)
                bc_genome[bc] = g
    if not bcs:
        sys.exit(f"no barcodes for tier {args.tier}")

    by_genome = defaultdict(list)
    for bc in bcs:
        by_genome[bc_genome[bc]].append(bc)
    genomes = sorted(by_genome)

    keep = None
    if args.genomes:
        with open(args.genomes) as fh:
            keep = {l.strip() for l in fh if l.strip()}
        recips_pool = [b for b in bcs if bc_genome[b] in keep]
        if not recips_pool:
            sys.exit(f"no barcodes for tier {args.tier} in {args.genomes}")
    else:
        recips_pool = bcs[:]

    ani = defaultdict(dict)
    with open(os.path.join(root, "truth", "genome_pairs.tsv")) as fh:
        next(fh)
        for line in fh:
            g1, g2, a = line.rstrip("\n").split("\t")
            a = float(a)
            ani[g1][g2] = a
            ani[g2][g1] = a

    group = {}
    sg = os.path.join(root, "truth", "strain_groups.tsv")
    if os.path.exists(sg):
        with open(sg) as fh:
            next(fh)
            for line in fh:
                g, grp, _n = line.rstrip("\n").split("\t")
                group[g] = grp

    def same_group(a, b):
        return group.get(a, a) == group.get(b, b)

    lo, hi = (args.bands + [98.0, 99.0])[:2]
    band_names = ["far", "mid", "near"]

    def band_of(a):
        return "far" if a < lo else ("mid" if a < hi else "near")

    # candidate partners per genome, bucketed by band
    partners = {}
    for g in genomes:
        buckets = {b: [] for b in band_names}
        for h in genomes:
            if h == g or same_group(g, h):
                continue
            buckets[band_of(ani.get(g, {}).get(h, 0.0))].append(h)
        partners[g] = buckets

    print(f"tier {args.tier}: {len(bcs)} barcodes, {len(genomes)} genomes")
    if keep:
        print(f"  recipients restricted to {len(keep)} genomes "
              f"-> {len(recips_pool)} barcodes")

    # ---- doublets -----------------------------------------------------------
    # Recipients are chosen per band from those genomes that actually HAVE a
    # partner in that band, rather than cycling bands blindly and falling back.
    # Blind rotation collapses into whichever band is commonest -- with these
    # genomes that is 'far', which is the least informative one. Conditioning on
    # availability keeps 'near' at a genuinely hard ANI while still filling the
    # quota.
    n_dbl = int(round(len(recips_pool) * args.doublet_rate))
    quota = {b: n_dbl // 3 for b in band_names}
    for b in band_names[:n_dbl - sum(quota.values())]:
        quota[b] += 1

    pool = recips_pool[:]
    random.shuffle(pool)
    used = set()
    picks = []
    for band in ("near", "mid", "far"):          # scarcest band first
        want = quota[band]
        for bc in pool:
            if want <= 0:
                break
            if bc in used or not partners[bc_genome[bc]][band]:
                continue
            used.add(bc)
            picks.append((bc, band))
            want -= 1
        if want > 0:
            print(f"  note: only {quota[band] - want}/{quota[band]} recipients "
                  f"have a '{band}' partner available")

    plan = []
    for i, (bc, band) in enumerate(picks):
        g1 = bc_genome[bc]
        g2 = random.choice(partners[g1][band])
        donor = random.choice(by_genome[g2])
        plan.append(dict(kind="doublet", new=f"DBL{args.prefix}{i+1:06d}",
                         recipient=bc, donor=donor, g1=g1, g2=g2,
                         stratum=band, ani=ani.get(g1, {}).get(g2, 0.0),
                         frac=0.5))

    remaining = [b for b in pool if b not in used]

    # ---- ambient ------------------------------------------------------------
    n_amb = int(round(len(recips_pool) * args.ambient_rate))
    amb_recipients = random.sample(remaining, min(n_amb, len(remaining)))
    for i, bc in enumerate(amb_recipients):
        g1 = bc_genome[bc]
        frac = args.ambient_fracs[i % len(args.ambient_fracs)]
        donors = []
        for _ in range(3):
            cand = [g for g in genomes if not same_group(g, g1)]
            donors.append(random.choice(by_genome[random.choice(cand)]))
        plan.append(dict(kind="ambient", new=f"AMB{args.prefix}{i+1:06d}",
                         recipient=bc, donor=",".join(donors), g1=g1,
                         g2=",".join(sorted({bc_genome[d] for d in donors})),
                         stratum=f"frac{frac}", ani=0.0, frac=frac))

    # ---- truth table --------------------------------------------------------
    tpath = os.path.join(root, "truth", f"injected_{args.tier}.tsv")
    mode = "a" if (os.path.exists(tpath) and args.prefix) else "w"
    with open(tpath, mode) as fh:
        if mode == "w":
            fh.write("new_bc\tkind\trecipient_bc\tdonor_bc\tgenome_major\t"
                     "genome_minor\tstratum\tani\tfraction\n")
        for p in plan:
            fh.write(f"{p['new']}\t{p['kind']}\t{p['recipient']}\t{p['donor']}\t"
                     f"{p['g1']}\t{p['g2']}\t{p['stratum']}\t{p['ani']:.2f}\t"
                     f"{p['frac']}\n")

    counts = defaultdict(int)
    for p in plan:
        if p["kind"] == "doublet":
            counts[p["stratum"]] += 1
    print(f"  doublets {sum(counts.values())}  " +
          "  ".join(f"{b}={counts[b]}" for b in band_names))
    print(f"  ambient  {len(amb_recipients)}  fracs {args.ambient_fracs}")
    print(f"  plan -> {tpath} ({'appended' if mode == 'a' else 'written'})")
    if args.dry_run:
        print("  --dry-run: no FASTQ written")
        return

    # ---- write reads --------------------------------------------------------
    # Concatenated gzip members form a valid gzip stream, so a doublet is a
    # byte copy of two .fq.gz files -- no decompression, no recompression.
    # Only the ambient spike needs the donor decompressed, and only the first
    # few percent of it. This is roughly two orders of magnitude faster than
    # reading and rewriting every record.
    def raw_copy(src, out):
        with open(src, "rb") as fh:
            shutil.copyfileobj(fh, out, length=1 << 20)

    def count_reads(path):
        n = 0
        with gzip.open(path, "rb") as fh:
            for _ in fh:
                n += 1
        return n // 4

    n_reads_cache = {}
    for j, p in enumerate(plan, 1):
        if j % 100 == 0:
            print(f"  writing {j}/{len(plan)}", flush=True)
        for mate in (1, 2):
            dest = os.path.join(outdir, f"{p['new']}{mate}.fq.gz")
            if os.path.exists(dest):
                continue
            rec = os.path.join(sc, f"{p['recipient']}{mate}.fq.gz")
            with open(dest + ".tmp", "wb") as out:
                raw_copy(rec, out)
                if p["kind"] == "doublet":
                    raw_copy(os.path.join(sc, f"{p['donor']}{mate}.fq.gz"), out)
                else:
                    key = (p["recipient"], mate)
                    if key not in n_reads_cache:
                        n_reads_cache[key] = count_reads(rec)
                    donors = p["donor"].split(",")
                    take = max(1, int(n_reads_cache[key] * p["frac"]
                                      / len(donors)))
                    for d in donors:
                        gz = gzip.GzipFile(fileobj=out, mode="wb",
                                           compresslevel=1)
                        for i, r in enumerate(read_fastq(
                                os.path.join(sc, f"{d}{mate}.fq.gz"))):
                            if i >= take:
                                break
                            gz.write(r.encode())
                        gz.close()
            os.replace(dest + ".tmp", dest)

    print(f"  wrote {len(plan) * 2} files to {outdir}")


if __name__ == "__main__":
    main()
