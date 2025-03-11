# parastar

Bash script that uses GNU Parallel and STAR to map paired-end fastq.gz files to a reference transcriptome using STAR.
The default parameters are defined for the _Drosophila melanogaster_ dm6-based r6.62 transcriptome and may be changed for other species.


## Installation

First, install GNU Parallel and STAR using mamba:
```bash
mamba create -n star_env -c bioconda -c conda-forge star parallel
conda activate star_env
```

Then, clone this repository:
```bash
git clone https://github.com/ccarloscr/parastar.git
```

### Download the fasta and gtf files of the reference genome
This script uses a reference genome to map fastq files, for which it needs the fasta and the gtf files of the mapped genome. 

Use the code below to download the fasta file of the dm6 genome from UCSC:
```bash
mkdir -p parastar/Genomes/dm6
cd parastar/Genomes/dm6
wget http://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
gunzip dm6.fa.gz
```

Use the code below to download the latest version of the dm6 gtf file (20/02/2025):
```bash
wget http://ftp.flybase.org/genomes/dmel/current/gtf/dmel-all-r6.62.gtf.gz
gunzip dmel-all-r6.62.gtf.gz
```

### Folder containing input files
Compressed fastq.gz input files should be located in a directory named fastq_files, within the parastar repository:
```bash
mkdir -p parastar/fastq_files
```

