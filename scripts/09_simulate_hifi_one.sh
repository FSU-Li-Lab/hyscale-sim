#!/bin/bash
# Simulate PacBio HiFi reads for ONE genome, with fully traceable headers.
#
# HEADER CONTRACT -- this is the point of the script.
# Every read is named:
#
#     <genome>_<contigindex>_<readnum>
#     e.g. GCA_000167875.2_ASM16787v2_genomic_0002_147
#
# where <contigindex> is the 1-based, zero-padded position of the source contig
# in the genome's .fna, in file order. That is exactly the format the existing
# truth chain depends on: 07_longread_contigs.py joins <contigindex> against
# truth/contig_truth.tsv to recover the contig accession and hence plasmid vs
# chromosome. Verified against 6,514,591 reads of the previous dataset with a
# 100% resolution rate, so do not change the separator or the padding.
#
# pbsim already writes one BAM per contig named <prefix>_<NNNN>.bam, numbered by
# reference order, so <contigindex> comes straight from the filename. CCS names
# reads S<zmw>/<n>/ccs; the <n> becomes <readnum>.
#
# Pipeline: pbsim (multi-pass CLR subreads) -> ccs (consensus -> HiFi) ->
#           samtools fastq -> rename -> concatenate -> gzip
#
# MAF files are deleted: 700 MB per genome would be ~340 GB across 473, and
# exact read coordinates are cheaper to recover afterwards by mapping the HiFi
# reads back to their own reference (see 10_longread_coords.sh).
#
# Usage: simulate_hifi_one.sh <genome.fna> <outdir>
set -euo pipefail

FNA=${1:?usage: simulate_hifi_one.sh <genome.fna> <outdir>}
OUTDIR=${2:?usage: simulate_hifi_one.sh <genome.fna> <outdir>}

DEPTH=${DEPTH:-15}
PASSES=${PASSES:-10}
MINRQ=${MINRQ:-0.99}
MINPASS=${MINPASS:-3}
LENMEAN=${LENMEAN:-15000}
CCS_THREADS=${CCS_THREADS:-4}
ERRHMM=${ERRHMM:?set ERRHMM to the ERRHMM-SEQUEL.model path}

GENOME=$(basename "$FNA" .fna)
# pbsim runs inside a temp dir, so a relative genome path would break there
FNA=$(readlink -f "$FNA")
[ -f "$FNA" ] || { echo "[FAIL] genome not found: $FNA" >&2; exit 1; }
OUTDIR=$(mkdir -p "$OUTDIR" && readlink -f "$OUTDIR")
DEST="$OUTDIR/${GENOME}.hifi.fq.gz"
[ -s "$DEST" ] && { echo "[skip] $GENOME"; exit 0; }

WORK=$(mktemp -d "${TMPDIR:-/tmp}/hifi_${GENOME}_XXXXXX")
trap 'rm -rf "$WORK"' EXIT

# CCS writes per-thread scratch into TMPDIR, naming it from the OUTPUT file's
# basename: thread.<n>_0.<output>.bam. Every genome here produces hifi_0001.bam,
# so concurrent jobs sharing one TMPDIR collide on that scratch name and
# silently corrupt each other -- surfacing as
#     "BAM reader ERROR: cannot read from corrupted file ... probably truncated"
# Pointing TMPDIR at this genome's own workdir isolates it. This also keeps CCS
# scratch out of any shared temp directory other jobs may be using.
export TMPDIR="$WORK"

cd "$WORK"
pbsim --strategy wgs --method errhmm --errhmm "$ERRHMM" \
      --genome "$FNA" --depth "$DEPTH" --length-mean "$LENMEAN" \
      --pass-num "$PASSES" --prefix sim > pbsim.log 2>&1 || {
    echo "[FAIL] pbsim $GENOME" >&2; cat pbsim.log >&2; exit 1; }

rm -f sim_*.maf.gz sim_*.maf sim_*.ref

shopt -s nullglob
n_ok=0
n_fail=0
for bam in sim_*.bam; do
    idx=${bam#sim_}; idx=${idx%.bam}          # 0001, 0002, ...

    if ! ccs "$bam" "hifi_${idx}.bam" \
        --min-passes "$MINPASS" --min-rq "$MINRQ" \
        -j "$CCS_THREADS" --report-file "ccs_${idx}.txt" \
        > "ccs_${idx}.log" 2>&1; then
        echo "[warn] ccs failed on $GENOME contig $idx" >&2
        echo "------- ccs stderr -------" >&2
        tail -20 "ccs_${idx}.log" >&2
        echo "--------------------------" >&2
        # keep the workdir so the failure can be reproduced
        KEEP="$OUTDIR/failed_${GENOME}_${idx}"
        mkdir -p "$KEEP"
        cp "ccs_${idx}.log" "$KEEP/" 2>/dev/null || true
        echo "$WORK" > "$KEEP/workdir.txt"
        trap - EXIT                            # do NOT delete the workdir
        n_fail=$((n_fail + 1))
        continue
    fi
    n_ok=$((n_ok + 1))

    # S<zmw>/<n>/ccs  ->  <genome>_<contigindex>_<n>
    samtools fastq "hifi_${idx}.bam" 2>/dev/null \
      | awk -v g="$GENOME" -v i="$idx" '
          NR%4==1 { n=split($1,a,"/"); print "@" g "_" i "_" a[2]; next }
          NR%4==3 { print "+"; next }
          { print }' \
      >> all.fq
done

[ -s all.fq ] || { echo "[FAIL] no reads for $GENOME" >&2; exit 1; }
if [ "$n_fail" -gt 0 ]; then
    echo "[FAIL] $GENOME: $n_fail of $((n_ok + n_fail)) contigs failed CCS." >&2
    echo "       Refusing to write a partial genome -- reads would be missing" >&2
    echo "       for entire contigs and the coverage truth would be wrong." >&2
    exit 1
fi

gzip -c all.fq > "${DEST}.tmp" && mv "${DEST}.tmp" "$DEST"

reads=$(( $(wc -l < all.fq) / 4 ))
bases=$(awk 'NR%4==2 {b+=length($0)} END {print b}' all.fq)
echo "[done] $GENOME contigs=$n_ok reads=$reads bases=$bases"
