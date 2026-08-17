#!/bin/bash
# Simulate HiFi reads for ALL genomes, then merge and verify.
#
# Resumable: any genome whose .hifi.fq.gz already exists is skipped, so an
# interrupted run continues where it stopped.
#
# Usage:
#   bash 09_simulate_hifi_all.sh            # full run
#   bash 09_simulate_hifi_all.sh --test 3   # first 3 genomes only
set -euo pipefail

ROOT=${ROOT:-$(pwd)}
JOBS=${JOBS:-8}
export DEPTH=${DEPTH:-15}
export PASSES=${PASSES:-10}
export MINRQ=${MINRQ:-0.99}
export MINPASS=${MINPASS:-3}
export LENMEAN=${LENMEAN:-15000}
export CCS_THREADS=${CCS_THREADS:-4}
# pbsim3 ships its HMM models inside the conda env; locate rather than hard-code
if [ -z "${ERRHMM:-}" ]; then
  ERRHMM=$(find "${CONDA_PREFIX:-/usr}" -name "ERRHMM-SEQUEL.model" 2>/dev/null | head -1)
fi
[ -n "$ERRHMM" ] && [ -f "$ERRHMM" ] || {
  echo "set ERRHMM to the path of ERRHMM-SEQUEL.model (ships with pbsim3)" >&2
  exit 1
}
export ERRHMM

SCRIPTS=$(cd "$(dirname "$0")" && pwd)
cd "$ROOT"
OUT=pacbio_hifi
mkdir -p "$OUT"

TEST=0
[ "${1:-}" = "--test" ] && TEST=${2:-3}

for b in pbsim ccs samtools; do
  command -v "$b" >/dev/null || { echo "$b not on PATH" >&2; exit 1; }
done
[ -f "$ERRHMM" ] || { echo "model not found: $ERRHMM" >&2; exit 1; }

ls genome/*.fna > "$OUT/genomes.txt"
if [ "$TEST" -gt 0 ]; then
  head -n "$TEST" "$OUT/genomes.txt" > "$OUT/genomes.test.txt"
  mv "$OUT/genomes.test.txt" "$OUT/genomes.txt"
fi
N=$(wc -l < "$OUT/genomes.txt")
echo "[sim] $N genomes  depth=$DEPTH passes=$PASSES min-rq=$MINRQ  jobs=$JOBS"
echo "[sim] estimated ${DEPTH}x: ~$((DEPTH * 37 * N / JOBS / 3600))h wall"

parallel -j "$JOBS" --bar --joblog "$OUT/joblog.txt" \
  "bash $SCRIPTS/09_simulate_hifi_one.sh {} $OUT" :::: "$OUT/genomes.txt"

# ---------------------------------------------------------------- merge
echo "[merge] concatenating gzip members"
mkdir -p pacbio_hifi_merged
cat "$OUT"/*.hifi.fq.gz > pacbio_hifi_merged/pbsim_long_read.fastq.gz

# ---------------------------------------------------------------- verify
echo "[verify] header contract"
zcat pacbio_hifi_merged/pbsim_long_read.fastq.gz | awk 'NR%4==1' | sed 's/^@//' \
  > "$OUT/all_headers.txt"
TOTAL=$(wc -l < "$OUT/all_headers.txt")
echo "  reads: $TOTAL"

# every header must be <genome>_<4-digit contig index>_<readnum>
BAD=$(grep -cvE '_[0-9]{4}_[0-9]+$' "$OUT/all_headers.txt" || true)
echo "  malformed headers: $BAD"

# genome field must resolve to a real genome, and contig index must be in range
awk -F'\t' '{n[$2]++} END{for(g in n) print g"\t"n[g]}' truth/contig_truth.tsv \
  | sort > "$OUT/ncontig.tsv"
sed -E 's/_[0-9]{4}_[0-9]+$//' "$OUT/all_headers.txt" | sort | uniq -c \
  | awk '{print $2"\t"$1}' | sort > "$OUT/perg.tsv"
echo "  genomes represented: $(wc -l < "$OUT/perg.tsv") (expected $N)"

sed -E 's/_[0-9]+$//' "$OUT/all_headers.txt" \
  | awk '{n=split($0,a,"_"); idx=a[n]+0; g=$0; sub(/_[0-9]{4}$/,"",g);
          if (idx>m[g]) m[g]=idx} END{for(g in m) print g"\t"m[g]}' \
  | sort > "$OUT/maxidx.tsv"
join -t $'\t' "$OUT/maxidx.tsv" "$OUT/ncontig.tsv" > "$OUT/idxcheck.tsv"
echo "  max contig index <= n_contigs: $(awk -F'\t' '$2<=$3' "$OUT/idxcheck.tsv" | wc -l) / $(wc -l < "$OUT/idxcheck.tsv")"
awk -F'\t' '$2>$3 {print "  OVERFLOW: "$0}' "$OUT/idxcheck.tsv" | head

# EVERY contig must have produced reads. A CCS failure on one contig silently
# removes it entirely -- headers stay well-formed and index ranges stay valid,
# so the checks above pass while the chromosome is missing. This is the check
# that catches it.
echo "[verify] per-contig completeness"
sed -E 's/_[0-9]+$//' "$OUT/all_headers.txt" | sort -u > "$OUT/seen_contigs.txt"
awk -F'\t' '{n[$2]++} END{for(g in n) for(i=1;i<=n[g];i++) printf "%s_%04d\n", g, i}' \
  truth/contig_truth.tsv | sort > "$OUT/expected_contigs.txt"
# restrict expectation to the genomes actually simulated
cut -f1 "$OUT/perg.tsv" | sed 's/$/_/' > "$OUT/simulated.txt"
grep -F -f "$OUT/simulated.txt" "$OUT/expected_contigs.txt" | sort \
  > "$OUT/expected_sim.txt"
MISSING=$(comm -23 "$OUT/expected_sim.txt" "$OUT/seen_contigs.txt" | wc -l)
echo "  contigs expected: $(wc -l < "$OUT/expected_sim.txt")"
echo "  contigs with reads: $(wc -l < "$OUT/seen_contigs.txt")"
echo "  MISSING CONTIGS: $MISSING"
if [ "$MISSING" -gt 0 ]; then
  echo "  ---- missing (first 20) ----"
  comm -23 "$OUT/expected_sim.txt" "$OUT/seen_contigs.txt" | head -20
  echo "  Do NOT use this dataset. Check pacbio_hifi/failed_*/ for CCS logs."
fi

# achieved depth per genome, as a sanity check on yield
echo "[verify] achieved depth"
zcat pacbio_hifi_merged/pbsim_long_read.fastq.gz \
  | awk 'NR%4==1{h=$0; sub(/^@/,"",h); sub(/_[0-9]{4}_[0-9]+$/,"",h); g=h}
         NR%4==2{b[g]+=length($0)} END{for(x in b) print x"\t"b[x]}' \
  | sort > "$OUT/bases.tsv"
awk -F'\t' '{s[$2]+=$4} END{for(g in s) print g"\t"s[g]}' truth/contig_truth.tsv \
  | sort > "$OUT/glen.tsv"
join -t $'\t' "$OUT/bases.tsv" "$OUT/glen.tsv" \
  | awk -F'\t' -v d="$DEPTH" '{c=$2/$3; printf "  %-50s %.1fx\n", $1, c;
      if (c < d*0.5) print "    LOW: expected ~" d "x"}'

echo
echo "Next:"
echo "  mv pacbio pacbio_clr_backup"
echo "  mv pacbio_hifi_merged pacbio"
echo "  python3 scripts/07_longread_contigs.py --truth truth"
echo "  rm -f work/longreads.fa work/longreads.mmi   # rebuild index on HiFi"
