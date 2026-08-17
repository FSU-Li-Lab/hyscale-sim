#!/bin/bash

FASTA_DIR=/data10/xli/ref/Enterobacterales/fna_for_analysis
OUTDIR=/data10/xli/ref/Enterobacterales/fna_for_analysis/art_single_cell

# Number of simulated cells (barcodes) per genome
NCELL=30

# Per-cell coverages
COVERAGES=("0.01" "0.1" "0.2" "0.4" "0.8" "1" "10")

barcode_id=1

for cov in "${COVERAGES[@]}"
do

    echo "==============================="
    echo "Coverage = ${cov}X"
    echo "==============================="

    mkdir -p "${OUTDIR}/${cov}x"

    for fasta in "${FASTA_DIR}"/*.fna
    do

        genome=$(basename "$fasta" .fna)

        mkdir -p "${OUTDIR}/${cov}x/${genome}"

        echo "Genome: $genome"

        for ((i=1; i<=NCELL; i++))
        do

            barcode=$(printf "BC%06d" "$barcode_id")
            barcode_id=$((barcode_id+1))

            art_illumina \
                -ss HS25 \
                -p \
                -na \
                -l 150 \
                -f "$cov" \
                -m 350 \
                -s 50 \
                -i "$fasta" \
                -o "${OUTDIR}/${cov}x/${genome}/${barcode}"

        done

    done

done