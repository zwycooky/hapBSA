#!/usr/bin/perl

use strict;
my ($input_vcf,$output) = @ARGV[0,1];
my $Usage = "\n\t$0 <whatshap vcf> <output>
\n";
die $Usage unless (@ARGV == 2);

open IN,'<',"$input_vcf" or die;
open OUT,'>',"$output" or die;
while (<IN>) {
	chomp;
	if (/\A#/) { next };
	my ($chr,$pos,$ref,$alt,$info) = (split)[0,1,3,4,9];
	my ($geno,$block) = (split /:/,$info)[0,-1];
	if ($geno eq '0|1') {
		print OUT "$chr\t$pos\t$ref\t$alt\t$block\n";
	}elsif ($geno eq '1|0') {
		print OUT "$chr\t$pos\t$alt\t$ref\t$block\n";
	}
}
close IN;
close OUT;

