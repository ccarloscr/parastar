#!/bin/bash

# Author: Carlos Camilleri-Robles
# Contact: carloscamilleri@hotmail.com
# Version: ***
# This script uses STAR to map paired-end reads to Drosophila dm6 genome

set -e

## Variable set up
FASTQ_DIR = "rnaseq/fastq_files"
GENOME_FASTA = "rnaseq/genome/dm6.fasta"
GENOME_INDEX = "rnaseq/genome/dm6_index"
OUTPUT_DIR = "rnaseq/STAR_results"
THREADS = 12

## Create output folder
mkdir -p "$OUTPUT_DIR"

# Generate dm6 index files
if [ ! -d "$GENOME_INDEX" ]; then
    echo "Generating dm6 genome index..."
    STAR    --runThreadN $THREADS \
            --runMode genomeGenerate \
            --genomeDir "$GENOME_INDEX" \
            --genomeFastaFiles "$genome_fasta"
fi

# Mapping reads to dm6
echo "Initializing mapping of FASTQ files..."
for read1 in "$FASTQ_DIR"/_read1.fast1; do
    sample_id = $(basename "$read1" | sed 's/_read1.fastq//')
    read2 = "$FASTQ_DIR/${sample_id}_read2.fastq"
    if [ -f "$read2" ]; then
        echo "Processing $sample_id..."
        STAR    --runThreadN $THREADS \
                --genomeDir "$GENOME_INDEX/" \
                --readFilesIn "$read1" "$read2" \
                --outFileNamePrefix "$OUTPUT_DIR/${sample_id}_" \
                --outSAMtype BAM SortedByCoordinate
    else
        echo "ERROR: No $read2 file for $sample_id"
    fi
done

echo "All mappings complete."
