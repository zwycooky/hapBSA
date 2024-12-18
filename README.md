# hapBSA
Haplotype-based Bulked segregant analysis

## Introduction
hapBSA, a haplotype-based bulked segregant analysis tool for QTL mapping in half-sib populations. HapBSA uses maternal haplotype information to phase the short reads of the two sample pools and focuses on calculating the frequency bias of maternal haplotypes between two sample pools.
## Requirements
SAMtools v1.16.1 or higher  
BWA v0.7.17-r1188 or higher  
Hisat2 v2.2.1 or higher  
Whatshap v1.7 or higher  
Perl v5.30.0 or higher  
Perl module Parallel::ForkManager v2.03  

## Step by step tutorial
## 1. download the scripts and test data
```
# download test data
wget https://figshare.com/ndownloader/files/51266927 -O ->> test_data.tar.gz
tar zxf test_data.tar.gz
```

## 2. build docker image
```
# you can find a Dockerfile in https://github.com/zwycooky/hapBSA, and use it to build a docker container
docker build -t hapbsa-container .
```
Three directories, maternal_fq, pool_fq and test_genome, will be generated after this step, and containing fastq files of meternal accession, fastq files of the two sample pools and a fasta file of test genome, respectively.

## 3. Mapping parent reads to reference genome
```

```

