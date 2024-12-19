#!/usr/bin/perl

use strict;
use Cwd 'abs_path';

my ($sam,$phasing_data,$outprefix,$snp_nums_for_reads) = @ARGV[0,1,2,3];
my $Usage = "\n\t$0 <bam> <phased hap file> <output prefix> <snp numbers for reads>
\n";
die $Usage unless (@ARGV == 4);

#my $time1 = time();
## read hap data ##
my ($hap,$block);
open HAP,'<',"$phasing_data" or die "Error: Cannot open hap file:$!";
while (<HAP>) {
	chomp;
	my ($contig,$pos,$snp1,$snp2,$block_id) = (split)[0,1,2,3,4];
	$block = $block_id;
	#push @{$hap->{$contig}},"$pos\t$snp1\t$snp2\t$block";
	push @{$hap->{$contig}}, { pos => $pos, snp1 => $snp1, snp2 => $snp2, block => $block_id };
}
close HAP;

## start phasing pacbio/NGS reads ##
my ($pre_reads_id,@tmp,$fir);
open my $sam_file,"samtools view -q 45 $sam|" or die "Error: Cannot open sam/bam file:$!";

#open NOPHASING,'>',"$outprefix.nophase.txt";
open PHASE1,'>',"$outprefix.hap1.txt";
open PHASE2,'>',"$outprefix.hap2.txt";
#open LOWACC,'>',"$outprefix.low.accuracy.reads.txt";

while (<$sam_file>) {
	chomp;

	my ($reads_id,$mapQ) = (split /\t/,$_)[0,1];

	if ($fir == 0) {
		$fir = 1;
		push @tmp, $_;
		$pre_reads_id = $reads_id;
	}elsif ($pre_reads_id ne $reads_id) {
		my $sam_line = $tmp[0];
		my $phasing_res = &phasing_reads(@tmp);
		my $ifphased = (split /\t/,$phasing_res)[0];	
		#my ($reads_seq,$reads_q) = (split /\t/,$tmp[0])[9,10];
		
		if ($ifphased > 0) {
			my ($phased_hap,$acc,$reads_chr,$reads_start,$reads_end,$hap_block) = (split /\t/,$phasing_res)[0,1,2,3,4,5];
			if ($acc >= 90 && $phased_hap == 1) {
				print PHASE1 "$reads_chr\t$reads_start\t$reads_end\t$phased_hap\t$hap_block\n";
				#print "$reads_chr\t$reads_start\t$reads_end\t$phased_hap\n";
			}elsif ($acc >= 90 && $phased_hap == 2) {
				print PHASE2 "$reads_chr\t$reads_start\t$reads_end\t$phased_hap\t$hap_block\n";
				#print "$reads_chr\t$reads_start\t$reads_end\t$phased_hap\n";
			}
		}
	
		@tmp = ();
		push @tmp, $_;
		$pre_reads_id = $reads_id;		
	}else{
		push @tmp, $_;
	}
}
close $sam_file;

my $sam_line = $tmp[0];
my $phasing_res = &phasing_reads(@tmp);
my $ifphased = (split /\t/,$phasing_res)[0];
#my ($reads_seq,$reads_q) = (split /\t/,$tmp[0])[9,10];

if ($ifphased > 0) {
	my ($phased_hap,$acc,$reads_chr,$reads_start,$reads_end,$hap_block) = (split /\t/,$phasing_res)[0,1,2,3,4,5];
       	if ($acc >= 90 && $phased_hap == 1) {
		print PHASE1 "$reads_chr\t$reads_start\t$reads_end\t$phased_hap\t$hap_block\n";
       	}elsif ($acc >= 90 && $phased_hap == 2) {
		print PHASE2 "$reads_chr\t$reads_start\t$reads_end\t$phased_hap\t$hap_block\n";
      	}
}

#close NOPHASING;
close PHASE1;
close PHASE2;
#close LOWACC;

#my $time2 = time();
#my $time = ($time2 - $time1);
#print "$time s\n";

## sub progrem ##

sub phasing_reads {
	
	my @tmp = @_;
	my $full_reads = 0;
	my ($reads_name,$not_phase);
	my $n_tmp = @tmp;
	my ($accuracy,$snp1_count,$snp2_count,$total_count,$reads_len,$reads_chr,$reads_start,$reads_end,$max_block,$max_count);
	
	foreach (@tmp) {
		my ($reads_id,$contig,$start,$MQ,$cigar,$seq) = (split /\t/,$_)[0,2,3,4,5,9];
		#print "$reads_id\n";
		#print length($seq) . "\n";
		# get full reads #
		if ($full_reads == 0) {
			$full_reads = $seq;
			$reads_name = $reads_id;
			$reads_len = length($full_reads);
		}
		$reads_start = $start;
		$reads_chr = $contig;

		# find overlap with phasing region #
		if (exists $hap->{$contig}) {
			# get start and end position in contig with an alignment #
			my @cigar = &get_s_e_pos($start,$cigar);
			my $end = shift @cigar;
			$reads_end = $end;

			#print "end: $end\n\n";
			
			my @tmp_hap = @{$hap->{$contig}};
			my $search_s = &binarySearch($start,\@tmp_hap);
			my $search_e = &binarySearch($end,\@tmp_hap);
			@tmp_hap = @tmp_hap[($search_s-1)..($search_e+1)];
			
			## keep the biggest hap block for phaing ##
			my ($hap_block,%block_count);
			if (defined $block) {
				foreach (@tmp_hap) {
					my ($pos,$snp1,$snp2,$tmp_block) = ($_->{pos}, $_->{snp1}, $_->{snp2}, $_->{block});
					if ($pos > $end) { last };
					push @{$hap_block->{$tmp_block}}, {pos => $pos, snp1=> $snp1, snp2 => $snp2};
					$block_count{$tmp_block} ++;
				}
				$max_block = ();
				$max_count = 0;
				foreach (sort keys %block_count) {
					my $hap_count = $block_count{$_};
					if ($hap_count > $max_count) {
						$max_count = $hap_count;
						$max_block = $_;
					}
				}
				if ($max_count >= $snp_nums_for_reads) {
					@tmp_hap = @{$hap_block->{$max_block}};
				}else{
					@tmp_hap = ();
				}
			}
			
			my $overlap = 0;
			foreach (@tmp_hap) {
				my ($pos,$snp1,$snp2) = ($_->{pos}, $_->{snp1}, $_->{snp2});
				if ($pos >= $start && $pos <= $end) {
					$overlap = 1;
					## get base in reads ##
					my $relative_pos = $pos - $start + 1;
					my $target_base_num = &get_base_in_target($relative_pos,@cigar);
					if ($target_base_num eq 'no') { next };
					my $target_base = substr($seq,$target_base_num-1,1);
					#print "$contig\t$pos\t$target_base\t$snp1\t$snp2\n";
					
					if ($target_base eq $snp1) {
						$snp1_count ++;
						$total_count ++;
					}elsif ($target_base eq $snp2){
						$snp2_count ++;
						$total_count ++;
					}else {
						$total_count ++;
					}
					
				}elsif ($pos > $end) {
					last;
				}
			}
			
			if ($overlap == 0) {
				$not_phase += 1;
			}
			
			#print "$contig\t$start\t$end\n"
		}else{
			$not_phase += 1;
		}
	}
	
	if ($total_count >= $snp_nums_for_reads) {
		if ($snp1_count > $snp2_count) {
			$accuracy = sprintf("%.2f",$snp1_count / $total_count * 100);
			return("1\t$accuracy\t$reads_chr\t$reads_start\t$reads_end\t$max_block");
		}elsif ($snp1_count < $snp2_count) {
			$accuracy = sprintf("%.2f",$snp2_count / $total_count * 100);
			return("2\t$accuracy\t$reads_chr\t$reads_start\t$reads_end\t$max_block");
		}else{
			return("0\t0\t$reads_chr\t$reads_start\t$reads_end");
		}
	}else{
		return("0\t0\t$reads_chr\t$reads_start\t$reads_end");
	}
}

sub get_base_in_target {
	my $pos = shift @_;
	my @cigar = @_;
	
	#print "$pos\n";
	
	## position on reads ##
	my $base_num = 0;
	## position on contig ##
	my $contig_num = 0;
	
	foreach (@cigar) {
		if (/S/) {
			my $num = (split /S/,$_)[0];
			$base_num += $num;
		}elsif (/H/) {
			next;
		}elsif (/I/) {
			my $num = (split /I/,$_)[0];
			
			$base_num += $num;
		}elsif (/D/) {
			my $num = (split /D/,$_)[0];
			if ($pos > $contig_num && $pos < $contig_num + $num) {
				## pos in the Deletion region: is an ERROR ##
				## warn "$pos $contig_num $num in the Deletion region: is an ERROR:$!";
				return("no")
			}
			$contig_num += $num;
		}elsif (/N/) {
			my $num = (split /N/,$_)[0];
                        if ($pos > $contig_num && $pos < $contig_num + $num) {
                                ## pos in the Deletion region: is an ERROR ##
                                ## warn "$pos $contig_num $num in the Deletion region: is an ERROR:$!";
                                return("no")
                        }
                        $contig_num += $num;
		}elsif (/M/) {
			my $num = (split /M/,$_)[0];
			## find the pos ##
			if ($pos > $contig_num && $pos < $contig_num + $num) {
				my $add_num = $pos - $contig_num;
				my $reads_pos = $add_num + $base_num;
				#print "$reads_pos\n";
				return($reads_pos);
			}else{
				$base_num += $num;
				$contig_num += $num;
			}
		}else{
			die "Cannot recongnize $_:$!";
		}
	}
}

sub get_s_e_pos {
	my ($start_pos,$cigar) = @_;
	## split cigar ##
	my @cigar = (split /(\d+\w)/,$cigar);
	@cigar = grep (!/\A\s*\z/,@cigar);
	
	my $total_len = 0;
	foreach (@cigar) {
		if (/M/) {
			my $num = (split /M/,$_)[0];
			$total_len += $num;
		}elsif (/D/) {
			my $num = (split /D/,$_)[0];
			$total_len += $num;
		}elsif (/N/) {
			my $num = (split /N/,$_)[0];
                        $total_len += $num;
		}
	}
	
	my $end_pos = $start_pos + $total_len - 1;
	
	unshift @cigar,$end_pos;
	return(@cigar);
}

sub binarySearch {
	my ($pos,$arr) = @_;
	my ($low, $high) = (0, $#$arr);

	while ($low <= $high) {
        	my $mid = int(($low + $high) / 2);
        	my $mid_pos = $arr->[$mid]->{pos};
        	if ($mid_pos == $pos) {
            		return $mid;
        	} elsif ($mid_pos < $pos) {
            		$low = $mid + 1;
        	} else {
            		$high = $mid - 1;
        	}
    	}
    	return $low;
}



