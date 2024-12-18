# hapBSA
Haplotype-based Bulked segregant analysis

## Introduction
hapBSA, a haplotype-based bulked segregant analysis tool for QTL mapping in half-sib populations. HapBSA uses maternal haplotype information to phase the short reads of the two sample pools and focuses on calculating the frequency bias of maternal haplotypes between two sample pools.
## Requirements
SAMtools v1.16.1 or higher  
BWA v0.7.17-r1188 or higher  
Hisat2 v2.2.1 or higher  
Whatshap v1.7 or higher  
Perl module Parallel::ForkManager v2.03  

## Step by step tutorial
## 1. download the scripts and test data
```
# set the working directory to 'hapBSA_Dir'
mkdir -p hapBSA_Dir/{00parent_bam,01pool_bam,02parent_SNP}
export HAPBSA_DIR=`realpath hapBSA_Dir`
# download test data
wget https://figshare.com/ndownloader/files/51266927 -O ->> test_data.tar.gz
tar zxf test_data.tar.gz
mv maternal_fq pool_fq test_genome $HAPBSA_DIR
```
Three directories, 'maternal_fq', 'pool_fq' and 'test_genome', will be generated after this step, and containing fastq files of meternal accession, fastq files of the two sample pools and a fasta file of test genome, respectively.

## 2. build docker image
```
# you can find a Dockerfile in https://github.com/zwycooky/hapBSA, and use it to build a docker container
docker build -t hapbsa-container .
```

## 3. Mapping short reads to reference genome
```
## build index for reference genome
docker run --user $(id -u):$(id -g) -v ${HAPBSA_DIR}/test_genome/:/input hapbsa-container bwa index /input/test.genome.fa
docker run --user $(id -u):$(id -g) -v ${HAPBSA_DIR}/test_genome/:/input hapbsa-container samtools faidx /input/test.genome.fa
docker run --user $(id -u):$(id -g) -v ${HAPBSA_DIR}/test_genome/:/input hapbsa-container samtools dict -o /input/test.genome.dict /input/test.genome.fa
docker run --user $(id -u):$(id -g) -v ${HAPBSA_DIR}/test_genome/:/input hapbsa-container hisat2-build /input/test.genome.fa /input/test.genome.fa

## mapping parent reads to reference genome
docker run --user $(id -u):$(id -g) \
	-v ${HAPBSA_DIR}/test_genome/:/genome \
	-v ${HAPBSA_DIR}/maternal_fq/:/input \
	-v ${HAPBSA_DIR}/00parent_bam/:/output hapbsa-container \
	02_mapit.pl \
	-1 /input/mat.pair1.fq.gz \
	-2 /input/mat.pair2.fq.gz \
	-f /genome/test.genome.fa \
	-b /genome/test.genome.fa \
	-p 10 \
	-o /output/parent
	
## mapping pool1 and pool2 reads to reference genome (RNA-seq)
docker run --user $(id -u):$(id -g) \
	-v ${HAPBSA_DIR}/test_genome/:/genome \
	-v ${HAPBSA_DIR}/pool_fq/:/input \
	-v ${HAPBSA_DIR}/01pool_bam/:/output hapbsa-container \
	hisat2.pl \
	/genome/test.genome.fa \
	/input/ \
	10 \
	/output/
```

## 4. phasing parent SNPs by whatshap
```

```