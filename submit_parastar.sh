#!/usr/bin/env bash

source parastar.conf

mkdir -p "$LOG_DIR"

sbatch \
  --job-name="$SBATCH_JOB_NAME" \
  --partition="$SBATCH_PARTITION" \
  --nodes="$SBATCH_NODES" \
  --cpus-per-task="$SBATCH_CPUS_PER_TASK" \
  --mem="$SBATCH_MEM" \
  --time="$SBATCH_TIME" \
  --output="${LOG_DIR}/${SBATCH_OUTPUT}" \
  --error="${LOG_DIR}/${SBATCH_ERROR}" \
  --wrap="./parastar.sh parastar.conf"