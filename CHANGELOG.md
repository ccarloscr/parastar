# Changelog


## Version [2.0.0] — 2026-05-12

Complete rewrite of the pipeline. The original single-script design has been replaced by a modular, configurable architecture.

### Added

- `parastar.conf`: dedicated user-editable configuration file — no script editing required
- `submit_parastar.sh`: SLURM submission wrapper that reads all job parameters from `parastar.conf`
- `parastar_env.yaml`: reproducible conda environment (`STAR 2.7.11b`, `parallel`, `pigz`, `seqkit`)
- Dry-run mode (`DRY_RUN=true`) to validate config without executing STAR
- Skip-completed logic (`SKIP_COMPLETED=true`) for resumable runs
- Per-sample temp directories (`--outTmpDir`) to prevent parallel collisions
- gzip integrity validation (`VALIDATE_GZIP=true`) with per-sample failure reporting
- GNU Parallel joblog (`parallel_joblog.tsv`) recording per-sample exit codes
- Genome index rebuild control: `CREATE_GENOME_INDEX` and `FORCE_REBUILD_INDEX` flags
- Runtime resource sanity check: warns if `PARALLEL_JOBS × STAR_THREADS` exceeds `SBATCH_CPUS_PER_TASK`
- Zero-sample guard: exits with error if no FASTQ files are found matching the pattern
- Configurable file naming patterns (`R1_PATTERN`, `R2_PATTERN`)
- Configurable STAR extra arguments (`STAR_EXTRA_ARGS` array in conf)
- `EXIT` trap for temp directory cleanup and conda deactivation on pipeline end or error
- Explicit error handling on `source` calls for config and conda init script
- `pigz` added to dependency validation check

### Changed

- Conda environment renamed from `star_env` to `parastar_env` for consistency
- SLURM log paths now derived from `LOG_DIR` in conf rather than hardcoded
- Genome index check now validates specific required files (`Genome`, `SA`, `SAindex`, `chrLength.txt`) rather than just directory existence
- Thread allocation logic redesigned: threads per job and number of parallel jobs are now set explicitly in conf rather than computed dynamically
- `gunzip -c` replaced with `pigz -dc` for faster parallel decompression
- Sample R2 discovery now uses pattern substitution on R1 path instead of constructing path from `FASTQ_DIR`

### Removed

- Inline SBATCH directives (replaced by `submit_parastar.sh`)
- Automatic GTF chr-prefix modification (users should pre-process their GTF if needed)
- Dynamic thread-per-job calculation based on sample count
- `--quantMode TranscriptomeSAM` default (replaced by `--quantMode GeneCounts` in `STAR_EXTRA_ARGS`)
- `--genomeSAindexNbases` hardcoded value (users set genome-specific parameters in conf)

---

## Version [1.0.0] — 2025-09-26

Initial release.

- Single Bash script with inline SBATCH directives for SLURM submission
- Paired-end FASTQ mapping with STAR via GNU Parallel
- Automatic genome index generation if index directory is absent
- Dynamic thread allocation: redistributes CPUs based on number of samples found
- Optional GTF chr-prefix modification for chromosome name compatibility
- Basic input validation for FASTQ directory, genome FASTA, and GTF file
