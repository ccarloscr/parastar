# parastar

Bash script using GNU Parallel and STAR to map paired-end fastq files to reference genome using STAR.
Default parameters are defined for _Drosophila melanogaster_ dm6 genome.


## Installation

First, install GNU Parallel and STAR using mamba:
```bash
mamba create -n star_env -c bioconda -c conda-forge star
```

Then, clone this repository:
```bash
git clone https://github.com/ccarloscr/parastar.git
cd parastar
```

### Download the fasta and gtf files of the reference genome
This script uses a reference genome to map fastq files, for which it needs the fasta and the gtf files of the mapped genome. 

Use the code below to download the fasta file of the dm6 genome from UCSC:
```bash
mkdir -p Genomes/dm6
cd Genomes/dm6
wget http://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
gunzip dm6.fa.gz
```

Use the code below to download the latest version of the dm6 gtf file (20/02/2025):
```bash
wget http://ftp.flybase.org/genomes/dmel/current/gtf/dmel-all-r6.62.gtf.gz
gunzip dmel-all-r6.62.gtf.gz
```

## Configuration





