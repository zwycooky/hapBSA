#!/usr/bin/perl

use strict;
use Getopt::Std;
use Cwd 'abs_path';
use Parallel::ForkManager;

our ($opt_1,$opt_2,$opt_p,$opt_r,$opt_a,$opt_o,$opt_w,$opt_s,$opt_d,$opt_e,$opt_m,$opt_N,$opt_n,$opt_D);
getopt ("1:2:p:r:a:o:w:s:d:e:m:N:n:D:");

my $Usage = "Usage:\n  $0 -OPTIONS VALUES

options:
--input options
	-1 FILE			bam file of pool1
	-2 FILE			bam file of pool2
	-p FILE			hap file
	-r REFGENOME	reference genome for mapping

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
	-a INT          cpus cores used for the analysis [default: 10]     
\n";

die $Usage unless ($opt_1 && $opt_2 && $opt_p && $opt_r && $opt_o && $opt_e && $opt_m);

my $bam1                = $opt_1;
my $bam2                = $opt_2;
my $hap_file            = $opt_p;
my $refgenome           = $opt_r;
my $outprefix           = $opt_o;
my $cpu                 = (defined $opt_a)?$opt_a:6;
my $window_len          = (defined $opt_w)?$opt_w:1000000;
my $step_len            = (defined $opt_s)?$opt_s:600000;
my $depth               = (defined $opt_d)?$opt_d:10;
my $block_read_nums	= (defined $opt_N)?$opt_N:12;
my $snps_for_reads	= (defined $opt_n)?$opt_n:3;
my $read_nums_window	= (defined $opt_D)?$opt_D:50;
my $sep_script          = $opt_e;
my $snpMapper           = $opt_m;

die "ERROR:$!" unless (-e $sep_script);
die "ERROR:$!" unless (-e $snpMapper);

$sep_script = abs_path($sep_script);
$snpMapper  = abs_path($snpMapper);

my $time1 = time();

## read hap file ##
my ($hap,%chr_len,$hap_block);
open HAP,'<',"$hap_file" or die;
while (<HAP>) {
	chomp;
	my ($chr,$pos,$hap1,$hap2,$block) = (split)[0,1,2,3,4];
	if (defined ($block)) {
		$hap_block = 1;
		push @{$hap->{$chr}},"$pos\t$hap1\t$hap2\t$block";
	}else{
		push @{$hap->{$chr}},"$pos\t$hap1\t$hap2";
	}
	$chr_len{$chr} = $pos;
}
close HAP;

## generate tmpdir ##
my $randnum = int(rand(100000))+1;
my $tmpdir = "tmpdir$randnum";
if (!-e $tmpdir) {
	mkdir $tmpdir;
}
$tmpdir = abs_path($tmpdir);

## generating sliding window ##
#my @command;
open OUT,'>',"$outprefix.hapBSA.tmp.txt" or die "ERROR with output:$!";
foreach (sort keys %chr_len) {
	my $chr = $_;
	my $chr_end = $chr_len{$chr};
	my $window = $window_len;
	my $step = $step_len;
	my $win_start = 0;
	my $win_end = $win_start + $window;
	
	my @sliding_window;
	while ($win_end < $chr_end) {
		push @sliding_window,"$chr\t$win_start\t$win_end";
		$win_start += $step;
                $win_end += $step;
	}
	push @sliding_window,"$chr\t$win_start\t$win_end";

	my $pm = Parallel::ForkManager->new($cpu);

	foreach (@sliding_window) {
		
		$pm->start and next;
		
		my ($tmp_chr,$s,$e) = (split /\t/,$_)[0,1,2];
		!system "samtools view -bh -q 45 -o $tmpdir/$tmp_chr.$s.p1.bam $bam1 $tmp_chr:$s-$e" or die "ERROR with samtools:$!";
		!system "samtools view -bh -q 45 -o $tmpdir/$tmp_chr.$s.p2.bam $bam2 $tmp_chr:$s-$e" or die "ERROR with samtools:$!";
		!system "samtools index $tmpdir/$tmp_chr.$s.p1.bam" or die "ERROR with samtools:$!";
		!system "samtools index $tmpdir/$tmp_chr.$s.p2.bam" or die "ERROR with samtools:$!";

		my @tmp_hap = @{$hap->{$tmp_chr}};
		open TMPHAP,'>',"$tmpdir/$tmp_chr.$s.hap.txt";
		foreach (@tmp_hap) {
			my $pos = (split)[0];
			if ($pos >= $s && $pos <= $e) {
				print TMPHAP "$chr\t$_\n";
			}elsif ($pos > $e) {
				last;
			}
		}
		close TMPHAP;
		
		## hapBSA ##
		my $prefix = "$tmp_chr.$s";
        	!system "perl $sep_script $tmpdir/$tmp_chr.$s.p1.bam $tmpdir/$tmp_chr.$s.hap.txt $tmpdir/$prefix.p1 $snps_for_reads" or die "ERROR with phasing reads:$!";
       		!system "perl $sep_script $tmpdir/$tmp_chr.$s.p2.bam $tmpdir/$tmp_chr.$s.hap.txt $tmpdir/$prefix.p2 $snps_for_reads" or die "ERROR with phasing reads:$!";
		## snp index ##
		chomp(my $window_indexED4 = `perl $snpMapper $refgenome $tmpdir/$tmp_chr.$s.p1.bam $tmpdir/$tmp_chr.$s.p2.bam $depth`);
		
		## hap block mode for reads based phasing ##
		if (defined $hap_block) {
			my ($hap_block_count);
			
			open PH,'<',"$tmpdir/$prefix.p1.hap1.txt" or die;
			while (<PH>) {
				chomp;
				my ($block_id) = (split)[-1];
				$hap_block_count->{$block_id}[0] ++;
			}
			close PH;
			
			open PH,'<',"$tmpdir/$prefix.p1.hap2.txt" or die;
                        while (<PH>) {
                                chomp;
                                my ($block_id) = (split)[-1];
                                $hap_block_count->{$block_id}[1] ++;
                        }
                        close PH;
			
			open PH,'<',"$tmpdir/$prefix.p2.hap1.txt" or die;
                        while (<PH>) {
                                chomp;
                                my ($block_id) = (split)[-1];
                                $hap_block_count->{$block_id}[2] ++;
                        }
                        close PH;

			open PH,'<',"$tmpdir/$prefix.p2.hap2.txt" or die;
                        while (<PH>) {
                                chomp;
                                my ($block_id) = (split)[-1];
                                $hap_block_count->{$block_id}[3] ++;
                        }
                        close PH;
			
			my ($hapIndex_hap,$blockNum,$hapA_reads_num,$hapB_reads_num);
			foreach (sort keys %{$hap_block_count}) {
				my ($p1h1,$p1h2,$p2h1,$p2h2);
				if (!defined $hap_block_count->{$_}[0]) {
					$p1h1 = 0;
				}else{
					$p1h1 = $hap_block_count->{$_}[0];
				}

				if (!defined $hap_block_count->{$_}[1]) {
                                        $p1h2 = 0;
                                }else{
                                        $p1h2 = $hap_block_count->{$_}[1];
                                }
				
				if (!defined $hap_block_count->{$_}[2]) {
                                        $p2h1 = 0;
                                }else{
                                        $p2h1 = $hap_block_count->{$_}[2];
                                }
				
				if (!defined $hap_block_count->{$_}[3]) {
                                        $p2h2 = 0;
                                }else{
                                        $p2h2 = $hap_block_count->{$_}[3];
                                }
				
				#print "$p1h1\t$p1h2\t$p2h1\t$p2h2\n";

				if ($p1h1+$p2h1 > $block_read_nums && $p1h2+$p2h2 > $block_read_nums && $p1h1 >= $p1h2) {
					$hapIndex_hap += ($p1h1 / ($p1h1+$p1h2)) - ($p2h1 / ($p2h1+$p2h2));
					$blockNum ++;
					$hapA_reads_num += $p1h1 + $p2h1;
					$hapB_reads_num += $p1h2 + $p2h2;
				}elsif ($p1h1+$p2h1 > $block_read_nums && $p1h2+$p2h2 > $block_read_nums && $p1h1 < $p1h2) {
					$hapIndex_hap += ($p1h2 / ($p1h1+$p1h2)) - ($p2h2 / ($p2h1+$p2h2));
					$blockNum ++;
					$hapA_reads_num += $p1h1 + $p2h1;
                                        $hapB_reads_num += $p1h2 + $p2h2;
				}
			}
			if ($blockNum > 0 && $hapA_reads_num > $read_nums_window && $hapB_reads_num > $read_nums_window) {
				$hapIndex_hap = $hapIndex_hap / $blockNum;
			}else{
				$hapIndex_hap = 'NA';
			}

			my $win_pos = ($s + $e) / 2;
			my $chr_num = $tmp_chr;
			#if ($tmp_chr =~ /(\d+)\z/) {
			#       $chr_num = $1;
			#}
			print OUT "$chr_num\t$win_pos\t$hapIndex_hap\t$window_indexED4\n";

		}else{
			chomp(my $win_p1_Areads_num = `wc -l $tmpdir/$prefix.p1.hap1.txt`);
			chomp(my $win_p1_Breads_num = `wc -l $tmpdir/$prefix.p1.hap2.txt`);
			chomp(my $win_p2_Areads_num = `wc -l $tmpdir/$prefix.p2.hap1.txt`);
                	chomp(my $win_p2_Breads_num = `wc -l $tmpdir/$prefix.p2.hap2.txt`);
			$win_p1_Areads_num = (split /\s+/,$win_p1_Areads_num)[0];
			$win_p1_Breads_num = (split /\s+/,$win_p1_Breads_num)[0];
			$win_p2_Areads_num = (split /\s+/,$win_p2_Areads_num)[0];
			$win_p2_Breads_num = (split /\s+/,$win_p2_Breads_num)[0];
			my $win_pos = ($s + $e) / 2;
			my $chr_num = $tmp_chr;
			#if ($tmp_chr =~ /(\d+)\z/) {
			#	$chr_num = $1;
			#}
		
			my ($win_p1_Areads_ratio,$win_p2_Areads_ratio,$hapIndex,$hapIndex_hap);
			if ($win_p1_Areads_num + $win_p2_Areads_num > $read_nums_window && $win_p1_Breads_num + $win_p2_Breads_num > $read_nums_window) {
				$win_p1_Areads_ratio = $win_p1_Areads_num / ($win_p1_Areads_num + $win_p1_Breads_num);
				$win_p2_Areads_ratio = $win_p2_Areads_num / ($win_p2_Areads_num + $win_p2_Breads_num);
				$hapIndex = abs($win_p1_Areads_ratio - $win_p2_Areads_ratio);
				$hapIndex_hap = $win_p1_Areads_ratio - $win_p2_Areads_ratio;
			}else{
				$hapIndex = 'NA';
				$hapIndex_hap = 'NA';
			}
			print OUT "$chr_num\t$win_pos\t$hapIndex_hap\t$window_indexED4\n";
		}

		unlink "$tmpdir/$tmp_chr.$s.p1.bam";
		unlink "$tmpdir/$tmp_chr.$s.p2.bam";
		unlink "$tmpdir/$tmp_chr.$s.p1.bam.bai";
                unlink "$tmpdir/$tmp_chr.$s.p2.bam.bai";
		unlink "$tmpdir/$tmp_chr.$s.hap.txt";
		unlink "$tmpdir/$prefix.p1.hap1.txt";
		unlink "$tmpdir/$prefix.p1.hap2.txt";
		unlink "$tmpdir/$prefix.p2.hap1.txt";
		unlink "$tmpdir/$prefix.p2.hap2.txt";
	
		$pm->finish;
	}
	$pm->wait_all_children();
}
`rm -rf $tmpdir`;
close OUT;

my @res;
open IN,'<',"$outprefix.hapBSA.tmp.txt" or die "Cannot open hapBSA.tmp.txt:$!";
while (<IN>) {
	chomp;
	push @res,$_;
}
close IN;

open OUT,'>',"$outprefix.hapBSA.sliding_window.txt" or die;
foreach (sort {(split /\t/,$a)[0] <=> (split /\t/,$b)[0] or (split /\t/,$a)[1] <=> (split /\t/,$b)[1]} @res ) {
	print OUT "$_\n";
}
close OUT;

unlink "$outprefix.hapBSA.tmp.txt";

my $time2 = time();
my $time  = ( $time2 - $time1 ) / 60;
print "$0 finished! Total time elapsed: $time min\n";





