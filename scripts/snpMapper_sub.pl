#!/usr/bin/perl

use strict;
my ($refgenome,$bam1,$bam2,$depth) = @ARGV[0,1,2,3];
my $Usage = "\n\t$0 <ref genome> <bam1> <bam2> <Depth>
\n";
die $Usage unless (@ARGV == 4);
#############################################

open my $fh_pileup,"samtools mpileup -q 45 -Q 20 -ABf $refgenome $bam1 $bam2|";
my ($snp_count,$index_sum,$ED4_sum,$index05_sum,$index01_sum,$ED405_sum,$ED401_sum,$index1_sum,$ED41_sum,%ED4_thre,%SNPindex_thre);
while (<$fh_pileup>) {
	chomp;
	
	next if /^\[/;
       	next if /^</;

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
	
	my $mutAfrq = $mut{'A'} / $mutcov1;
	my $mutCfrq = $mut{'C'} / $mutcov1;
	my $mutGfrq = $mut{'G'} / $mutcov1;
	my $mutTfrq = $mut{'T'} / $mutcov1;
	
	my $wtAfrq = $wt{'A'} / $wtcov1;
        my $wtCfrq = $wt{'C'} / $wtcov1;
        my $wtGfrq = $wt{'G'} / $wtcov1;
        my $wtTfrq = $wt{'T'} / $wtcov1;

	my $ED4 = sprintf "%.4f", ( $mut_freq - $wt_freq ) * ( $mut_freq - $wt_freq );

	if ($mut_freq > 0.9 && $wt_freq > 0.9) {
		next;
	}
	# delta SNP-index great than user-defined value will be considered as putative SNPs
	my $chr_num;
        if ($Chr =~ /(\d+)\z/) {
        	$chr_num = $1;
        }
	
	## random threshold ##
	my (@thre_index,@thre_ED4);
	foreach (1..1000) {
		my ($mthre,$wthre);
		for (1..25) {
			my $rand = rand();
			if ($rand >= 0.5) {
				$mthre ++;
			}
		}
		for (1..25) {
                        my $rand = rand();
                        if ($rand >= 0.5) {
                                $wthre ++;
                        }
                }
		$mthre = $mthre / 25;
		$wthre = $wthre / 25;

		my ($mut_rand,$wt_rand);
		foreach (1..$mutcov1) {
			my $rand = rand();
			if ($rand >= $mthre) {
				$mut_rand ++;
			}
		}
		foreach (1..$wtcov1) {
                        my $rand = rand();
                        if ($rand >= $wthre) {
                                $wt_rand ++;
                        }
                }
		my $rm1 = $mut_rand / $mutcov1;
		my $rm2 = 1 - $rm1;
		my $rw1 = $wt_rand / $wtcov1;
		my $rw2 = 1 - $rw1;
		my $rand_SNPindex = sprintf "%.4f", abs($rm1 - $rw1);
		my $rand_ED4 = sprintf "%.4f", ($rm1 - $rw1)*($rm1 - $rw1);
		push @thre_index,$rand_SNPindex;
		push @thre_ED4,$rand_ED4;
	}
	@thre_index = sort {$b <=> $a} @thre_index;
	@thre_ED4 = sort {$b <=> $a} @thre_ED4;
	foreach (0..99) {
		$SNPindex_thre{$_} += $thre_index[$_];
		$ED4_thre{$_} += $thre_ED4[$_];
	}

	# output #
	#print "$chr_num\t$Pos\t$refBase\t$mut_freq\t$wt_freq\t$deltaSNP\t$ED4\t$thre05_index\t$thre01_index\t$thre05_ED4\t$thre01_ED4\t$mutcov1\t$mut{'A'}\t$mut{'C'}\t$mut{'G'}\t$mut{'T'}\t$wtcov1\t$wt{'A'}\t$wt{'C'}\t$wt{'G'}\t$wt{'T'}\n";
	$snp_count ++;
	$index_sum += $deltaSNP;
	$ED4_sum += $ED4;
}
close $fh_pileup;

my ($index_window,$ED4_window,@res);
if ($snp_count == 0) {
	foreach (1..4) {
		push @res,'NA';
	}
	$index_window = 'NA';
	$ED4_window = 'NA';
}else{
	$index_window = sprintf "%.4f", ($index_sum / $snp_count);
	$ED4_window = sprintf "%.4f", ($ED4_sum / $snp_count);
	my @sig_0501 = (9,49);
	foreach (@sig_0501) {
		my $tmp_window = sprintf "%.4f", ($SNPindex_thre{$_} / $snp_count);
                push @res, $tmp_window;
        }
	foreach (@sig_0501) {
                my $tmp_window = sprintf "%.4f", ($ED4_thre{$_} / $snp_count);
                push @res, $tmp_window;
        }
}
my $res = join("\t",@res);
print "$index_window\t$ED4_window\t$res";

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

sub sum{
	my $sum = 0;
	for (@_){
		$sum += $_; 
	}
	return $sum;
}
