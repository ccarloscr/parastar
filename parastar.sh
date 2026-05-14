#!/usr/bin/env bash

# Author: Carlos Camilleri-Robles
# Contact: carloscamilleri@hotmail.com
# Version: 2026-05-12
# Description: Batch mapping of paired-end FASTQ files using STAR and GNU Parallel
# License: MIT License


# Set strict mode
set -Eeuo pipefail

########################################
# LOAD CONFIG
########################################

# Load configuration variables from parastar.conf
CONFIG_FILE="${1:-parastar.conf}"

# Check if config file exists
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

# Source the config file
source "$CONFIG_FILE" || error_exit "Failed to source config: $CONFIG_FILE"

########################################
# FUNCTIONS
########################################

# Logging function with timestamp
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# Error handling function
error_exit() {
    echo "ERROR: $*" >&2
    exit 1
}

# Check if a command exists
check_command() {
    command -v "$1" >/dev/null 2>&1 || \
        error_exit "Required command not found: $1"
}

# Validate that a file exists
validate_file() {
    [[ -f "$1" ]] || error_exit "Missing file: $1"
}

# Validate that a directory exists
validate_directory() {
    [[ -d "$1" ]] || error_exit "Missing directory: $1"
}

# Validate dependencies
validate_dependencies() {
    log "Checking dependencies"

    check_command STAR
    check_command parallel
    check_command find
    check_command gzip
    check_command pigz
}

# Prepare output, log, and temporary directories
prepare_directories() {
    mkdir -p "$OUTPUT_DIR"
    mkdir -p "$LOG_DIR"
    mkdir -p "$TMP_DIR"
}

# Validate input files and directories
validate_inputs() {
    log "Validating input files"

    validate_directory "$FASTQ_DIR"
    validate_file "$GENOME_FASTA"
    validate_file "$GTF_FILE"
}

# Validate that a gzip file is not corrupted
validate_gzip_file() {
    local file="$1"

    if [[ "$VALIDATE_GZIP" == true ]]; then
        if ! gzip -t "$file" 2>/dev/null; then
            log "WARNING: Corrupted gzip file, skipping: $file" >&2
            return 1
        fi
    fi
}

# Validate that the requested resources do not exceed SBATCH limits
validate_resources() {
    local total=$(( PARALLEL_JOBS * STAR_THREADS ))
    if (( total > SBATCH_CPUS_PER_TASK )); then
        log "WARNING: PARALLEL_JOBS ($PARALLEL_JOBS) x STAR_THREADS ($STAR_THREADS) = $total exceeds SBATCH_CPUS_PER_TASK ($SBATCH_CPUS_PER_TASK)"
    fi
}

# Build STAR genome index if it doesn't exist or if forced
build_genome_index() {

    # Check if genome index already exists and is valid, unless FORCE_REBUILD_INDEX is true
    if [[ "$CREATE_GENOME_INDEX" != true ]]; then
        log "Skipping genome index creation"
        return
    fi

    # Check for required STAR index files to determine if the index is already built
    local required_files=(Genome SA SAindex chrLength.txt)
    local valid_index=true

    # Check if all required index files exist
    for f in "${required_files[@]}"; do
        if [[ ! -f "$GENOME_INDEX_DIR/$f" ]]; then
            valid_index=false
        fi
    done

    # If FORCE_REBUILD_INDEX is true, we will rebuild the index regardless of existing files
    if [[ "$FORCE_REBUILD_INDEX" == true ]]; then
        valid_index=false
    fi

    # If the index is valid, we can skip building it
    if [[ "$valid_index" == true ]]; then
        log "Genome index already exists"
        return
    fi

    log "Building STAR genome index"

    mkdir -p "$GENOME_INDEX_DIR"

    # Run STAR genomeGenerate to build the index
    STAR \
        --runMode genomeGenerate \
        --runThreadN "$STAR_THREADS" \
        --genomeDir "$GENOME_INDEX_DIR" \
        --genomeFastaFiles "$GENOME_FASTA" \
        --sjdbGTFfile "$GTF_FILE" \
        --sjdbOverhang "$STAR_OVERHANG" \
        --limitGenomeGenerateRAM "$STAR_RAM_LIMIT"
}

# Mapping function to be run in parallel
map_sample() {

    local r1="$1"
    local r2
    local sample
    local sample_output
    local sample_log

    # Derive R2 filename from R1 using the specified patterns
    r2="${r1/$R1_PATTERN/$R2_PATTERN}"

    # Check if R2 file exists
    if [[ ! -f "$r2" ]]; then
        echo "Missing R2 pair for: $r1"
        return 1
    fi

    # Validate that both R1 and R2 files are valid gzip files if validation is enabled
    validate_gzip_file "$r1" || return 1
    validate_gzip_file "$r2" || return 1

    # Extract sample name from R1 filename by removing the R1 pattern
    sample=$(basename "$r1" "$R1_PATTERN")

    # Define output and log file paths for the sample
    sample_output="$OUTPUT_DIR/${sample}"
    sample_log="$LOG_DIR/${sample}.log"

    # Check if the sample has already been processed and skip if SKIP_COMPLETED is true
    if [[ "$SKIP_COMPLETED" == true ]] && \
       [[ -f "${sample_output}Aligned.sortedByCoord.out.bam" ]]; then

        log "Skipping completed sample: $sample"
        return
    fi

    log "Mapping sample: $sample"

    # Dry run option to print the STAR command without executing it
    if [[ "$DRY_RUN" == true ]]; then
        echo "DRY RUN: STAR mapping for $sample"
        return
    fi

    # Run STAR for the sample, redirecting stdout and stderr to the sample log file
    STAR \
        --runThreadN "$STAR_THREADS" \
        --genomeDir "$GENOME_INDEX_DIR" \
        --readFilesIn "$r1" "$r2" \
        --readFilesCommand pigz -dc \
        --outTmpDir "$TMP_DIR/${sample}_tmp" \
        --outFileNamePrefix "$sample_output" \
        ${STAR_EXTRA_ARGS_SERIALIZED} \
        > "$sample_log" 2>&1
}

# Export functions and variables for GNU Parallel
export -f map_sample
export -f log
export -f validate_gzip_file
export -f error_exit

########################################
# MAIN
########################################

main() {
    # Set up trap for cleanup on exit
    trap 'log "Pipeline interrupted or finished. Cleaning up..."; rm -rf "$TMP_DIR"' EXIT
    
    log "Starting Parastar"
    validate_resources
    prepare_directories
    validate_dependencies
    validate_inputs
    build_genome_index
    log "Launching parallel mapping"

    # Export necessary variables for parallel execution 
    export \
        OUTPUT_DIR \
        LOG_DIR \
        TMP_DIR \
        GENOME_INDEX_DIR \
        STAR_THREADS \
        R1_PATTERN \
        R2_PATTERN \
        VALIDATE_GZIP \
        SKIP_COMPLETED \
        DRY_RUN \
        STAR_EXTRA_ARGS_SERIALIZED="${STAR_EXTRA_ARGS[*]}"   # Serialize STAR_EXTRA_ARGS for export to parallel

        # Find R1 files and map in parallel
        mapfile -d '' r1_files < <(find "$FASTQ_DIR" -name "*${R1_PATTERN}" -print0)

        # Check if any R1 files were found
        if [[ ${#r1_files[@]} -eq 0 ]]; then
            error_exit "No FASTQ files matching pattern '*${R1_PATTERN}' found in: $FASTQ_DIR"
        fi

        # Calculate sample count and log progress
        log "Found ${#r1_files[@]} sample(s) to process"
        # Stream the array of R1 files to GNU Parallel using null delimiters to handle spaces in filenames
        printf '%s\0' "${r1_files[@]}" |
        parallel -0 -j "$PARALLEL_JOBS" \
            --joblog "$LOG_DIR/parallel_joblog.tsv" \
            map_sample

    log "Pipeline completed successfully"
}

main "$@"

# Check for any failed samples in the parallel job log
if grep -qP '\t[1-9]\d*\t' "$LOG_DIR/parallel_joblog.tsv" 2>/dev/null; then
    log "WARNING: Some samples failed. Check $LOG_DIR/parallel_joblog.tsv"
fi
