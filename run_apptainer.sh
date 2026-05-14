#!/usr/bin/env bash

set -euo pipefail

# This script runs the Apptainer container for Parastar with the specified configuration file.
# If no configuration file is provided as an argument, it defaults to "parastar.conf".
CONFIG_FILE="${1:-parastar.conf}"

# Source the configuration file to load any necessary environment variables or settings.
source "$CONFIG_FILE"

# Run the Apptainer container with the specified configuration file
# The --cleanenv option ensures that the container environment is clean
# The --bind option mounts the current working directory into the container
# The --bind options also mount the directories containing the FASTQ files, genome FASTA file, and GTF file into the container
# The --pwd option sets the working directory inside the container to the current directory
# The parastar.sif file is the Apptainer image that contains the necessary environment and dependencies
# The parastar.sh script is executed inside the container with the provided configuration file as an argument
apptainer exec \
    --cleanenv \
    --bind "$PWD":"$PWD" \
    --bind "$FASTQ_DIR":"$FASTQ_DIR" \
    --bind "$(dirname "$GENOME_FASTA")":"$(dirname "$GENOME_FASTA")" \
    --bind "$(dirname "$GTF_FILE")":"$(dirname "$GTF_FILE")" \
    --pwd "$PWD" \
    parastar.sif \
    bash parastar.sh "$CONFIG_FILE"