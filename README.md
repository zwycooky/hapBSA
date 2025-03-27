# hapBSA
Haplotype-based Bulked segregant analysis

## Table of contents

- [Introduction](#introduction)
- [Requirements](#requirements)
- [Usage of hapBSA](#usage-of-hapbsa)
- [Quick start](#quick-start)
- [Format of input file](#format-of-input-file)
- [Format of output file](#format-of-output-file)
- [Step by step tutorial (docker)](#step-by-step-tutorial)
- [Plot hapBSA result](#plot-hapbsa-result)

## Introduction
hapBSA, a haplotype-based bulked segregant analysis tool for QTL mapping in half-sib populations. HapBSA uses maternal haplotype information to phase the short reads of the two sample pools and focuses on calculating the frequency bias of maternal haplotypes between two sample pools.

## Requirements
SAMtools v1.16.1 or higher  
BWA v0.7.17-r1188 or higher  
Hisat2 v2.2.1 or higher  
Whatshap v1.7 or higher  
Perl module Parallel::ForkManager v2.03  
Perl module List::Util

## Usage of hapBSA
```
Usage:
  scripts/hapBSA_V4.pl -OPTIONS VALUES

options:
--input options
	-1 FILE			bam file of pool1
	-2 FILE			bam file of pool2
	-p FILE			hap file
	-r REFGENOME		reference genome for mapping

--sub progrem
	-e SCRIPT		path of separating_reads_by_haplotype.binarySearch.hapBSA.block.pl
	-m SCRIPT		path of snpMapper_sub.pl
	
--output options         
	-o PREFIX		prefix of output file

--criteria options
	-w INT			window length (bp) [default: 1000000]
	-s INT			window step (bp)   [default: 600000]
	-d INT			minimum depth for calculating SNP-index [default: 10]
	-D INT			minimum phased read numbers of a window [default: 50]
	-N INT			minimum read numbers of hap block [default: 12]
	-n INT			minimum numbers of phased snps for a read [default: 3]

--perfromanat options
	-a INT			cpus cores used for the analysis [default: 10]
```

## Quick start
```
perl scripts/hapBSA_V4.pl \
	-1 pool1.bam \
	-2 pool2.bam \
	-p haplotype.txt \
	-r reference.fa \
	-e scripts/separating_reads_by_haplotype.binarySearch.hapBSA.block.pl \
	-m scripts/snpMapper_sub.pl \
	-o output_prefix
```
Before starting, please ensure that the required tools are correctly installed and added to the $PATH to run this tool.

## Format of input file
Example of 'haplotype.txt' file:
```
chr1    36142   A       C       36142
chr1    36243   T       C       36142
chr1    36315   A       G       36142
chr1    36316   A       T       36142
chr1    36333   C       A       36142
chr1    36335   C       G       36142
chr1    78839   A       G       78839
chr1    78884   C       T       78839
chr1    126320  G       T       126320
chr1    126404  T       C       126320
chr1    126432  C       A       126320
chr1    139327  C       T       139327
chr1    139365  T       C       139327
chr1    255193  G       A       255193
chr1    255936  C       T       255936
chr1    256020  T       A       255936
```
The columns indicate:  
1. chromosome ID
2. SNP position
3. haplotype A
4. haplotype B
5. haplotype block
  
**Notice:** The haplotype block (column 5) is optional if your haplotype information is obtained using PollenSeq or other chromosome-level phasing methods.  
'haplotype.txt' can be transfered by 'get_phased_SNP_from_whatshap.pl' from phased VCF file generated from whatshap
```
perl scripts/get_phased_SNP_from_whatshap.pl whatshap.vcf haplotype.txt
```

## Format of output file
The output file of hapBSA will be named as ${output_prefix}.hapBSA.sliding_window.txt
```
chr1    4100000 0.65752997002997	0.1223  0.0909  0.4772  0.3689  0.2366  0.1418
chr1    4700000 0.65752997002997	0.0537  0.0712  0.4620  0.3531  0.2199  0.1288
chr1    5300000 NA			0.0283  0.0641  0.4985  0.3847  0.2555  0.1524
chr1    5900000 NA			0.0743  0.0591  0.4998  0.3828  0.2574  0.1509
chr1    6500000 NA			0.0695  0.0740  0.4817  0.3709  0.2389  0.1417
chr1    7100000 NA			0.0433  0.0826  0.4650  0.3566  0.2229  0.1312
chr1    7700000 NA			-0.0028 0.0643  0.4654  0.3572  0.2228  0.1312
```
The columns are:  
1. chromosome ID
2. sliding window position
3. hap-index
4. SNP-index
5. ED4
6. 0.01 threshold for hap-index/SNP-index
7. 0.05 threshold for hap-index/SNP-index
8. 0.01 threshold for ED4
9. 0.05 threshold for ED4
  
**Notice:** Hap-index < 0 indicates the frequency of hap A in pool1 < frequency of hap A in pool2. **This is only meaningful if the haplotypes are chromosome-level phasing.**

## Step by step tutorial
## 1. download the scripts and test data
```
# set the working directory to 'hapBSA_Dir'
mkdir -p hapBSA_Dir/{00parent_bam,01pool_bam,02parent_SNP,parent_tmpdir}
export HAPBSA_DIR=`realpath hapBSA_Dir`

# download test data
wget https://figshare.com/ndownloader/files/51266927 -O ->> hapBSA_test_data.gz
tar zxf hapBSA_test_data.gz
mv maternal_fq pool_fq test_genome $HAPBSA_DIR

```
Three directories, 'maternal_fq', 'pool_fq' and 'test_genome', will be generated after this step, and containing fastq files of meternal accession, fastq files of the two sample pools and a fasta file of test genome, respectively.

## 2. build docker image
```
# you can find a Dockerfile in https://github.com/zwycooky/hapBSA, and use it to build a docker container
docker build -t hapbsa-container .
```
If your network connection is unstable, some errors may occur. Please retry the command 'docker build -t hapbsa-container .' until the Docker container is successfully built.

## 3. Mapping short reads to reference genome
```
## rename the seq id of  test.genome.fa. 
sed -i 's/lg1:1-30000000/chr1/' ${HAPBSA_DIR}/test_genome/test.genome.fa

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
## SNP calling of parent using GATK
docker run --user $(id -u):$(id -g) \
	-v ${HAPBSA_DIR}/00parent_bam:/input/ \
	-v ${HAPBSA_DIR}/parent_tmpdir:/tmpdir/ \
	-v ${HAPBSA_DIR}/02parent_SNP:/output/ \
	-v ${HAPBSA_DIR}/test_genome/:/genome/ \
	hapbsa-container SNP_caller.pl /input/ /tmpdir/ /output/ /genome/test.genome.fa 10
  
## merge SNP calling results
docker run --user $(id -u):$(id -g) -v ${HAPBSA_DIR}/parent_tmpdir:/tmpdir/ \
	-v ${HAPBSA_DIR}/02parent_SNP:/output/ \
	hapbsa-container merge_split_vcf.pl /tmpdir/PASS.vcf.list /output/Parent_merged.vcf

## phasing SNPs using whatshap
docker run --user $(id -u):$(id -g) \
	-v ${HAPBSA_DIR}/02parent_SNP:/input/ \
	-v ${HAPBSA_DIR}/test_genome:/genome/ \
	-v ${HAPBSA_DIR}/00parent_bam:/bam/ \
	hapbsa-container whatshap phase --ignore-read-groups \
	-o /input/parent.whatshap.ngs.phased.vcf \
	--reference=/genome/test.genome.fa \
	/input/Parent_merged.vcf \
	/bam/parent.sorted.marked.duplicates.bam

## format phased SNPs from whatshap 
docker run --user $(id -u):$(id -g) \
	-v ${HAPBSA_DIR}/02parent_SNP:/input/ \
	hapbsa-container get_phased_SNP_from_whatshap.pl \
	/input/parent.whatshap.ngs.phased.vcf \
	/input/parent.whatshap.ngs.phased.txt
```

## 5. run hapBSA
```
docker run --user $(id -u):$(id -g) \
	-v ${HAPBSA_DIR}/01pool_bam:/pool_bam/ \
	-v ${HAPBSA_DIR}/02parent_SNP:/snp/ \
	-v ${HAPBSA_DIR}/test_genome:/genome/ \
	-v ${HAPBSA_DIR}:/output/ \
	hapbsa-container hapBSA_V4.pl \
	-1 /pool_bam/A.sorted.bam \
	-2 /pool_bam/B.sorted.bam \
	-p /snp/parent.whatshap.ngs.phased.txt \
	-r /genome/test.genome.fa \
	-e separating_reads_by_haplotype.binarySearch.hapBSA.block.pl \
	-m snpMapper_sub.pl \
	-o /output/test.hapBSA	
```

## plot hapBSA result
```
Rscript scripts\hapBSA_plot_cmdArg.R example_file\11_87864065_ngs_block.hapBSA.sliding_window.txt example_plot.pdf
```
'11_87864065_ngs_block.hapBSA.sliding_window.txt' is an example file of output of hapBSA.  
**Notice:** Before you run this Rscript you need to transfer the chromosome id into numbers (e.g. chr12 -> 12) and remove unanchored contig/scafford in the output of hapBSA.
