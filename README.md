# parastar

parastar is a Bash script that uses GNU Parallel and STAR to efficiently map paired-end fastq.gz files to a reference transcriptome. It is optimized for _Drosophila melanogaster_ (dm6, release r6.62), but can be adapted for other species by modifying the input files and parameters.


## Installation

First, install the required tools using mamba (or conda):
```bash
mamba create -n star_env -c bioconda -c conda-forge star parallel
conda activate star_env
```

Then, clone this repository:
```bash
git clone https://github.com/ccarloscr/parastar.git
```

## Reference genome setup
This script requires both the FASTA and GTF files of the reference genome.

#### Download dm6 FASTA from UCSC
```bash
mkdir -p parastar/Genomes/dm6
cd parastar/Genomes/dm6
wget http://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
gunzip dm6.fa.gz
```

#### Download dm6 GTF from FlyBase (version r6.62, as of 20/02/2025):
```bash
wget http://ftp.flybase.org/genomes/dmel/current/gtf/dmel-all-r6.62.gtf.gz
gunzip dmel-all-r6.62.gtf.gz
```

## Input files
Place your compressed paired-end FASTQ files inside the [fastq_files](fastq_files) directory:
```bash
mkdir -p parastar/fastq_files
```
Each sample must have two files:
- sampleID_R1_001.fastq.gz
- sampleID_R2_001.fastq.gz

## Configuration
- Read length: 50 bp → change line 30 in parastar.sh if needed.
- Reference genome index will be generated automatically if not found.
- Output directory: [Results](Results) (created automatically).


