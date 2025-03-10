# star_mapping

Simple script used to map paired-end fastq files to reference genome using STAR.
Default parameters are defined for _Drosophila melanogaster_ dm6 genome.


## Installation

First, install STAR using mamba:
```bash
mamba create -n star_env -c bioconda -c conda-forge star
```

Then, clone this repository:
```bash
git clone https://github.com/ccarloscr/star_mapping.git
cd star_mapping
```

## Configuration

# Download the fasta genome from UCSC (default dm6 genome)
cd ~/star_mapping/Genomes/dm6
wget http://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
gunzip dm6.fa.gz







