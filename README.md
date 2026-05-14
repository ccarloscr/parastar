# parastar

![Apptainer](https://img.shields.io/badge/Apptainer-Compatible-green?logo=apptainer)
![Slurm](https://img.shields.io/badge/Workload-Slurm-blue?logo=linux)
![Built with Apptainer](https://img.shields.io/badge/Built%20with-Apptainer-00599c?logo=apptainer&logoColor=white)

**Batch RNA-seq alignment using STAR and GNU Parallel**

Parastar is a Bash pipeline for parallelised paired-end FASTQ mapping with [STAR](https://github.com/alexdobin/STAR). It handles genome index generation, sample discovery, and parallel alignment in a single configurable run, designed for HPC environments with SLURM.

---

## Features

- Parallel alignment of multiple samples via GNU Parallel
- Fully compatible with Apptainer/Singularity
- Automatic genome index generation (with rebuild/skip logic)
- Configurable via a single `.conf` file, no script editing needed
- SLURM submission wrapper included
- Dry-run mode for testing without executing STAR
- Skip-completed logic for resumable runs
- Per-job temp directories to avoid parallel collisions
- Joblog for post-run failure inspection



## Repository structure

```
parastar/
├── parastar.sh          # Main pipeline logic
├── parastar.conf        # User-editable configuration file
├── run_apptainer.sh     # Apptainer container wrapper
├── submit_parastar.sh   # SLURM submission wrapper
├── README.md
├── CHANGELOG.md
└── LICENSE
```


## Requirements

- Linux HPC with SLURM
- Apptainer (formerly Singularity)
- Reference genome in FASTA format
- GTF annotation file



## Installation

**1. Clone the repository**

```bash
git clone https://github.com/ccarloscr/parastar.git
cd parastar
```

**2. Get the Apptainer image**

```bash
# Download from Zenodo
wget -O parastar.sif https://zenodo.org/records/20187976/files/parastar.sif
```



## Configuration

All user-facing parameters are in `parastar.conf`. Edit this file before running, **do not edit `parastar.sh` directly**.

### Key parameters

| Section | Parameter | Description |
|---|---|---|
| SLURM | `SBATCH_CPUS_PER_TASK` | Total CPUs requested |
| SLURM | `SBATCH_MEM` | Memory per job |
| SLURM | `SBATCH_TIME` | Wall time limit |
| Input | `FASTQ_DIR` | Directory containing paired FASTQ files |
| Input | `GENOME_FASTA` | Path to reference genome FASTA |
| Input | `GTF_FILE` | Path to GTF annotation file |
| Output | `OUTPUT_DIR` | Directory for STAR alignment outputs |
| Output | `GENOME_INDEX_DIR` | Directory for STAR genome index |
| Output | `LOG_DIR` | Directory for logs |
| STAR | `STAR_THREADS` | Threads per STAR job |
| STAR | `STAR_OVERHANG` | `sjdbOverhang` = read length − 1 (default: 99 for 100 bp reads) |
| STAR | `STAR_RAM_LIMIT` | RAM limit for genome generation (bytes) |
| Parallel | `PARALLEL_JOBS` | Number of samples to align simultaneously |
| Option | `SKIP_COMPLETED` | Skip samples with existing BAM output |
| Option | `VALIDATE_GZIP` | Check gzip integrity before mapping |
| Option | `CREATE_GENOME_INDEX` | Build index if missing |
| Option | `FORCE_REBUILD_INDEX` | Force index rebuild even if it exists |
| Option | `DRY_RUN` | Print actions without executing STAR |

### Resource planning

Total CPUs used = `PARALLEL_JOBS × STAR_THREADS`. This should not exceed `SBATCH_CPUS_PER_TASK`. Parastar will warn you at runtime if there is a mismatch.

**Example** (default config): 2 jobs × 12 threads = 24 CPUs

### STAR overhang

`STAR_OVERHANG` should be set to **read length − 1**. The default is `99` (for 100 bp reads). Adjust this to match your data.

### File naming patterns

Parastar expects paired FASTQ files matching:
```
{sample}_R1_001.fastq.gz
{sample}_R2_001.fastq.gz
```

These patterns are configurable via `R1_PATTERN` and `R2_PATTERN` in `parastar.conf`.



## Usage

### On a SLURM cluster (recommended)

```bash
bash submit_parastar.sh
```

This submits `parastar.sh` as a SLURM job using the parameters defined in `parastar.conf`.

### Locally via container

```bash
bash run_apptainer.sh parastar.conf
```

### Dry run (test config without mapping)

Set `DRY_RUN=true` in `parastar.conf`, then run as above.



## Output

For each sample, STAR outputs are written to `OUTPUT_DIR/` with the sample name as prefix, e.g.:

```
results/
├── sampleA/Aligned.sortedByCoord.out.bam
├── sampleA/ReadsPerGene.out.tab
├── sampleA/Log.final.out
├── sampleB/...
logs/
├── sampleA.log
├── sampleB.log
├── parallel_joblog.tsv
└── slurm-JOBID.out/err
```

To inspect failed samples after a run:

```bash
awk '$7 != 0' logs/parallel_joblog.tsv
```


## License

MIT License. See [LICENSE](LICENSE) for details.
