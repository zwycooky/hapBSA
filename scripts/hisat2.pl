#!/usr/bin/perl

use strict;

my ($genome,$fqdir,$cpu,$res_dir) = @ARGV[0,1,2,3];
my $Usage = "\n\t$0 <genome> <fastq directory> <cpu> <results dir>
\n";
die $Usage unless (@ARGV == 4);
#####----- checking index -----######
if (!-e "$genome.2.ht2") {
	&command_mkbowtie_index($genome,$genome,$cpu);
}
######################################################
my $time = time();
chomp(my @fqlist = `find $fqdir -name "*.gz"`);
my $fq_counter = 0;

if (! -e $res_dir) {
	mkdir $res_dir;
}

my $reads;
foreach ( sort @fqlist ) {
	chomp;
	s/\s//g;
	if (/\A\s*\z/) { next }
	$fq_counter ++;
	if ($fq_counter == 1) {
		$reads = $_;
	} elsif ($fq_counter == 2) {
		&command_hisat($genome,$reads,$_,$cpu);
		$fq_counter = 0;
	}
}

sub command_mkbowtie_index {
        my ($refgenome,$hisat2_index,$cpu) = @_;
        system "hisat2-build -p $cpu $refgenome $hisat2_index";
}

sub command_hisat {
	my ($hisat_index,$s_left,$s_right,$cpu)  = @_;
	my $id = (split /\./,(split /\//,$s_left)[-1])[0];
        !system "hisat2 -x $hisat_index -1 $s_left -2 $s_right -p $cpu --dta --rg-id 1 --rg \"PL:BGI\" --rg \"LB:1\" --rg \"SM:$id\" --rg \"PU:1\" -S $res_dir/$id.sam 2> $res_dir/$id.hisat2_aln.out" or die "Error in hisat2:$!";
        !system "samtools view -bhS -@ $cpu $res_dir/$id.sam > $res_dir/$id.bam" or die "Error in samtools view (sam to bam):$!";
        !system "samtools sort -@ $cpu -o $res_dir/$id.sorted.bam $res_dir/$id.bam" or die "Error in samtools sort:$!";
        !system "samtools index -@ $cpu $res_dir/$id.sorted.bam" or die "Error in samtools index:$!";
        
	unlink "$res_dir/$id.sam"; unlink "$res_dir/$id.bam";
}

