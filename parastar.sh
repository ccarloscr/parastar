#!/bin/bash

# Author: Carlos Camilleri-Robles
# Contact: carloscamilleri@hotmail.com
# Version: Not tested yet
# This script uses GNU Parallel and STAR to map paired-end reads to a reference genome

#SBATCH --job-name=parastar       # Job name
#SBATCH --partition=irbio01       # Slurm queue
#SBATCH --nodes=1                 # Number of nodes
#SBATCH --ntasks=1                # Number of tasks
#SBATCH --cpus-per-task=24        # CPUs per task
#SBATCH --output=parastar_%j.out  # Output log file
#SBATCH --error=parastar_%j.err   # Error log file

set -e

## Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate star_env

## Variable set up
FASTQ_DIR="fastq_files"
GENOME_FASTA="Genomes/dm6/dm6.fa"
GENOME_INDEX="Genomes/dm6/dm6_index"
GTF_DIR="Genomes/dm6/dmel-all-r6.62.gtf"
OUTPUT_DIR="Results"
NUM_FASTQ=$(find "$FASTQ_DIR" -name "*_R1_001.fastq.gz" | wc -l)
THREADS=12

## Create output folder
mkdir -p "$OUTPUT_DIR"

## Edit the gtf file to match with fasta chromosome naming
if ! head -n 1 "$GTF_DIR" | awk '{print $1}' | grep -q '^chr'; then
    sed -i 's/^\([^ ]*\)/chr\1/' "$GTF_DIR"
fi


## Check the existence of input directories
if [ ! -d "$FASTQ_DIR" ]; then
    echo "ERROR: FASTQ directory $FASTQ_DIR not found."
    exit 1
fi

if [ ! -f "$GENOME_FASTA" ]; then
    echo "ERROR: Genome FASTA file $GENOME_FASTA not found."
    exit 1
fi

if [ ! -f "$GTF_DIR" ]; then
    echo "ERROR: GTF file $GTF_DIR not found."
    exit 1
fi


# Generate genome index files
if [ ! -d "$GENOME_INDEX" ]; then
    echo "Generating genome index..."
    STAR    --runThreadN "$THREADS" \
            --runMode genomeGenerate \
            --genomeDir "$GENOME_INDEX" \
            --genomeFastaFiles "$GENOME_FASTA" \
            --sjdbGTFfile "$GTF_DIR"
fi


## Resources optimization
if [ "$NUM_FASTQ" -eq 0 ]; then
    echo "ERROR: No FASTQ files found in $FASTQ_DIR."
    exit 1
elif [ "$NUM_FASTQ" -lt "$THREADS" ]; then
    threads_per_job=$((THREADS / NUM_FASTQ))
    JOBS="$NUM_FASTQ"
else
    threads_per_job=1
    JOBS="$THREADS"
fi


## Check that threads_per_job is correctly set
if [ -z "$threads_per_job" ] || [ "$threads_per_job" -le 0 ]; then
    echo "ERROR: Invalid value for threads_per_job."
    exit 1
fi


## Define the mapping function
echo "Initializing mapping of FASTQ files..."
mapping() {
    read1=$1
    threads_per_job=$2
    sample_id=$(basename "$read1" | sed 's/_R1_001.fastq.gz//')
    read2="$FASTQ_DIR/${sample_id}_R2_001.fastq.gz"
    
    if [ -f "$read2" ]; then
        echo "Processing $sample_id using $threads_per_job threads..."
        STAR    --runThreadN "$threads_per_job" \
                --genomeDir "$GENOME_INDEX/" \
                --readFilesIn "$read1" "$read2" \
                --readFilesCommand gunzip -c \
                --outFileNamePrefix "$OUTPUT_DIR/${sample_id:0:5}_" \
                --outSAMtype BAM SortedByCoordinate
    else
        echo "WARNING: No read2 file for $sample_id. Skipping sample..."
        return 1
    fi
}


## Exports to Parallel
export -f mapping
export FASTQ_DIR GENOME_INDEX OUTPUT_DIR threads_per_job


## Use Parallel for the mapping
find "$FASTQ_DIR" -name "*_R1_001.fastq.gz" | parallel -j "$JOBS" mapping {} "$threads_per_job"

