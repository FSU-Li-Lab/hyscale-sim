#!/bin/bash
# Stage 1: build every truth table from metadata and read headers.
# Run once. Everything downstream joins against truth/ AFTER scoring.
set -euo pipefail

ROOT=${ROOT:-$(pwd)}
CLONE_THRESHOLD=${CLONE_THRESHOLD:-99.99}
SCRIPTS=$(cd "$(dirname "$0")" && pwd)

cd "$ROOT"
mkdir -p truth

echo "[1/5] barcode -> (tier, genome)"
awk '
  /:$/ {
    line=$0; sub(/:$/,"",line)
    n=split(line,a,"/"); genome=a[n]; tier=a[n-1]; next
  }
  /\.fq\.gz$/ {
    bc=$0; sub(/[12]\.fq\.gz$/,"",bc)
    key=tier"\t"bc"\t"genome
    if (!(key in seen)) { seen[key]=1; print key }
  }
' metadata/sc_bar_*.txt > truth/barcode_truth.tsv

# a barcode ID must never appear under two tiers, or the flat single_cell/
# directory has silently overwritten files and the reads are compromised
dup=$(cut -f2 truth/barcode_truth.tsv | sort | uniq -d | wc -l)
if [ "$dup" -ne 0 ]; then
  echo "  FATAL: $dup barcode IDs appear in more than one tier" >&2
  cut -f2 truth/barcode_truth.tsv | sort | uniq -d | head >&2
  exit 1
fi
echo "  $(wc -l < truth/barcode_truth.tsv) barcodes, no tier collisions"
cut -f1 truth/barcode_truth.tsv | sort | uniq -c | sed 's/^/  /'

echo "[2/5] long read -> genome   (decompresses ~80 GB, 20-30 min)"
zcat pacbio/pbsim_long_read.fastq.gz \
  | awk 'NR%4==1' \
  | sed 's/^@//' \
  | awk -v OFS='\t' '{ acc=$1; sub(/_[0-9]+_[0-9]+$/,"",acc); print $1, acc }' \
  > truth/longread_truth.tsv
echo "  $(wc -l < truth/longread_truth.tsv) long reads"
echo "  $(cut -f2 truth/longread_truth.tsv | sort -u | wc -l) source genomes"

echo "[3/5] contig -> (genome, chromosome|plasmid, length)"
for f in genome/*.fna; do
  b=$(basename "$f" .fna)
  awk -v g="$b" '
    /^>/ {
      if (acc != "") print acc"\t"g"\t"type"\t"len
      acc=substr($1,2)
      type = (tolower($0) ~ /plasmid|extrachromosomal|megaplasmid/) ? "plasmid" : "chromosome"
      len=0; next
    }
    { len += length($0) }
    END { if (acc != "") print acc"\t"g"\t"type"\t"len }
  ' "$f"
done > truth/contig_truth.tsv
echo "  $(wc -l < truth/contig_truth.tsv) contigs"
cut -f3 truth/contig_truth.tsv | sort | uniq -c | sed 's/^/  /'

echo "[4/5] skani matrices -> pair tables"
python3 "$SCRIPTS/parse_skani.py" \
  --genome-matrix  genome_ani.txt \
  --plasmid-matrix plasmid_ani.txt \
  --outdir truth

echo "[5/5] collapse indistinguishable genomes into strain groups"
python3 "$SCRIPTS/05_strain_groups.py" \
  --pairs truth/genome_pairs.tsv \
  --all-genomes truth/genome_list.txt \
  --outdir truth \
  --threshold "$CLONE_THRESHOLD"

echo
echo "truth/ contents:"
ls -la truth/
echo
echo "Next: review truth/pilot_genomes.txt, then"
echo "  python3 scripts/06_inject_negatives.py --root . --tier 1x"
echo "  bash scripts/02_link.sh 1x pilot"
