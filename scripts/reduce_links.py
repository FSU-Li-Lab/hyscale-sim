#!/usr/bin/env python3
"""
Reduce a stream of minimap2 PAF records to (barcode, long_read, n_hits, n_reads).

Reads PAF on stdin, writes TSV on stdout:
    barcode  long_read  n_hits  n_reads

  n_hits  = alignment records, secondaries included
  n_reads = distinct short reads (the quantity the hypergeometric test uses)

Query names are expected as <BARCODE>:<mate>:<n>, produced by 02_link.sh.

Memory note: minimap2 emits all records for a given query contiguously, so
distinct reads are counted by remembering only the last query name seen per
pair. The obvious implementation -- a set of query names per pair -- costs
roughly 100x more. At 1x coverage a 200-barcode batch produces ~84M alignment
records across ~800k pairs; sets would need ~10 GB, this needs ~150 MB. That
headroom is what lets --batch go high enough to amortise loading a 90 GB index.

If a future minimap2 ever interleaves queries, n_reads would over-count; n_hits
is unaffected. Check with:
    head -100000 file.paf | cut -f1 | uniq | sort | uniq -d | head
Empty output means queries are contiguous.
"""
import argparse
import sys
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-hits", type=int, default=1,
                    help="drop (barcode, long_read) pairs with fewer hits")
    ap.add_argument("--min-alen", type=int, default=80,
                    help="minimum alignment block length")
    ap.add_argument("--min-ident", type=float, default=0.95,
                    help="minimum nmatch/alen; simulated short reads are"
                         " near error-free, so this is safe")
    ap.add_argument("--mask", default=None,
                    help="optional BED of masked long-read intervals to ignore")
    ap.add_argument("--progress", type=int, default=0,
                    help="log every N million PAF records (0 = off)")
    args = ap.parse_args()

    mask = defaultdict(list)
    if args.mask:
        with open(args.mask) as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track")):
                    continue
                c, s, e = line.split()[:3]
                mask[c].append((int(s), int(e)))

    def masked(target, start, end):
        for s, e in mask.get(target, ()):
            if start < e and end > s:
                return True
        return False

    # key -> [n_hits, n_reads, last_qname]
    acc = {}
    n_in = n_kept = 0
    step = args.progress * 1_000_000

    for line in sys.stdin:
        f = line.split("\t", 12)
        if len(f) < 12:
            continue
        n_in += 1
        if step and n_in % step == 0:
            print(f"[reduce] {n_in//1_000_000}M records, {len(acc)} pairs",
                  file=sys.stderr, flush=True)
        alen = int(f[10])
        if alen < args.min_alen or int(f[9]) / alen < args.min_ident:
            continue
        tname = f[5]
        if mask and masked(tname, int(f[7]), int(f[8])):
            continue
        qname = f[0]
        key = (qname.split(":", 1)[0], tname)
        rec = acc.get(key)
        if rec is None:
            acc[key] = [1, 1, qname]
        else:
            rec[0] += 1
            if rec[2] != qname:
                rec[1] += 1
                rec[2] = qname
        n_kept += 1

    out = sys.stdout.write
    n_out = 0
    for (bc, lr), (hits, reads, _last) in acc.items():
        if hits < args.min_hits:
            continue
        out(f"{bc}\t{lr}\t{hits}\t{reads}\n")
        n_out += 1

    print(f"[reduce] paf_in={n_in} kept={n_kept} pairs_out={n_out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
