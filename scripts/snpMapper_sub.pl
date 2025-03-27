#!/usr/bin/perl

use strict;
use List::Util qw(sum);
use Math::Random qw(random_binomial);

my ($refgenome,$bam1,$bam2,$depth) = @ARGV[0,1,2,3];
my $Usage = "\n\t$0 <ref genome> <bam1> <bam2> <Depth>
\n";
die $Usage unless (@ARGV == 4);
#############################################

open my $fh_pileup,"samtools mpileup -q 45 -Q 20 -ABf $refgenome $bam1 $bam2|";
my ($snp_count,$sim_snp_count,$index_sum,$index05_sum,$index01_sum,$index1_sum,%SNPindex_thre);
while (<$fh_pileup>) {
	chomp;
	
	next if /^\[|^</;

	my @Contents = (split /\s+/,$_);

	# Omit any pos where only one sample mapped
	next if @Contents < 9;

	my ( $Chr, $Pos, $refBase ) = @Contents[ 0, 1, 2 ];

	# For convenience, Sample1 for mutant pool, Sample2 for wild pool
	my ( $mutcov0, $mutbases, $wtcov0, $wtbases ) = @Contents[ 3, 4, 6, 7 ];

	# Omit any intronic position?
	# next if ($Mut_Bases =~ />|</ or $Wt_Bases =~ />|</);

	# Omit low-coverage position

	# Calculate reads coverage for the position
	my @mutCounts = base_counter( $mutbases, $refBase );
	my @wtCounts  = base_counter( $wtbases,  $refBase );

	my %mut = (
		'A' => $mutCounts[0],
		'C' => $mutCounts[1],
		'G' => $mutCounts[2],
		'T' => $mutCounts[3]
	);
	$mut{$refBase} += $mutCounts[4];
	my %wt = (
		'A' => $wtCounts[0],
		'C' => $wtCounts[1],
		'G' => $wtCounts[2],
		'T' => $wtCounts[3]
	);
	$wt{$refBase} += $wtCounts[4];

	my $mutcov1 = $mut{'A'} + $mut{'C'} + $mut{'G'} + $mut{'T'};
	my $wtcov1  = $wt{'A'} + $wt{'C'} + $wt{'G'} + $wt{'T'};
	if ( $mutcov1 < $depth || $wtcov1 < $depth ) {
		next;
	}
	my @alleles =  sort { $mut{$b} <=> $mut{$a} or $wt{$b} <=> $wt{$a} } keys %mut;
	my $mut_freq = sprintf "%.4f", $mut{ $alleles[0] } / $mutcov1;    # 0.7
	my $wt_freq  = sprintf "%.4f", $wt{ $alleles[0] } / $wtcov1;      # 0.3
	my $deltaSNP = sprintf "%.4f", ( $mut_freq - $wt_freq );
	
	if ($mut_freq > 0.9 && $wt_freq > 0.9) {
		next;
	}
	# delta SNP-index great than user-defined value will be considered as putative SNPs
	my $chr_num;
        if ($Chr =~ /(\d+)\z/) {
        	$chr_num = $1;
        }
	
	## random threshold ##
	$sim_snp_count ++;
	if ($sim_snp_count <= 3000) {
    		my @thre_index;
    		$#thre_index = 999; # 预分配数组大小为1000

    		foreach my $i (0..999) {
        		my $mthre = random_binomial(1 ,25, 0.5) / 25;
        		my $wthre = random_binomial(1 ,25, 0.5) / 25;

        		my $mut_rand = random_binomial(1 ,$mutcov1, 1 - $mthre);
        		my $wt_rand  = random_binomial(1 ,$wtcov1, 1 - $wthre);

        		my $rand_SNPindex = abs( 
            			$mut_rand / $mutcov1 - $wt_rand / $wtcov1 
        		);
        		$thre_index[$i] = $rand_SNPindex;
    		}

    		@thre_index = sort { $b <=> $a } @thre_index;

    		my @sig_tmp0501 = (9, 49);
    		foreach my $sig (@sig_tmp0501) {
        		$SNPindex_thre{$sig} += $thre_index[$sig];
    		}
	}
=pod
	if ($sim_snp_count <= 3000) {
		my (@thre_index);
		foreach (1..1000) {
			my ($mut_rand, $wt_rand) = (0, 0);
			my $mthre = sum(map { rand() >= 0.5 } (1..25)) / 25;
			my $wthre = sum(map { rand() >= 0.5 } (1..25)) / 25;
			$mut_rand += sum(map { rand() >= $mthre } (1..$mutcov1));
			$wt_rand += sum(map { rand() >= $wthre } (1..$wtcov1));
		
			my $rm1 = $mut_rand / $mutcov1;
			my $rm2 = 1 - $rm1;
			my $rw1 = $wt_rand / $wtcov1;
			my $rw2 = 1 - $rw1;
			my $rand_SNPindex = sprintf "%.4f", abs($rm1 - $rw1);
			push @thre_index,$rand_SNPindex;
		}

		@thre_index = sort {$b <=> $a} @thre_index;
	
		my @sig_tmp0501 = (9,49);
		foreach (@sig_tmp0501) {
			$SNPindex_thre{$_} += $thre_index[$_];
		}
	}
=cut

	# output #
	$snp_count ++;
	$index_sum += $deltaSNP;
}
close $fh_pileup;

my ($index_window,@res);
if ($snp_count == 0) {
	foreach (1..2) {
		push @res,'NA';
	}
	$index_window = 'NA';
}else{
	$index_window = sprintf "%.4f", abs($index_sum / $snp_count);
	my @sig_0501 = (9,49);
	if ($sim_snp_count > 3000) {
		$sim_snp_count = 3000;
	}
	foreach (@sig_0501) {
		my $tmp_window = sprintf "%.4f", ($SNPindex_thre{$_} / $sim_snp_count);
                push @res, $tmp_window;
        }
}
my $res = join("\t",@res);
print "$index_window\t$res";

# base_counter: calculate base counts for each base in (A,C,G,T) order
sub base_counter {
	my ( $sample_bases, $refbase ) = @_;

	# Convert all dot and comma symbol to ref base
	#$sample_bases =~ s/\.|,/$refbase/gi;
	
	my $baseE = ($sample_bases =~ tr/.|,/./);
	# Remove patterns that represents INDELs
	while ( $sample_bases =~ /(.*?)[+-](\d+)[ATCG.,]+/ig ) {
		$sample_bases =~ s/(.*?)[+-](\d+)[ATCGNatcgn]{$2}(.*)/$1$3/i;
	}

	my $baseA = ($sample_bases =~ tr/A|a/A/);
	my $baseC = ($sample_bases =~ tr/C|c/C/);
	my $baseG = ($sample_bases =~ tr/G|g/G/);
	my $baseT = ($sample_bases =~ tr/T|t/T/);
	return ( $baseA, $baseC, $baseG, $baseT ,$baseE);
}

