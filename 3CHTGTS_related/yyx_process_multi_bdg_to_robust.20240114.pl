#!/usr/bin/perl

use strict;
use warnings;

my $version = '
Version: 0.1.0 (2024-01-14)
Author: Adam Yongxin Ye @ BCH
';
my $usage = "Usage: perl $0 <output_prefix> <input.multi.bed>
	<rep1.signif.merged.bed> <rep2.signif.merged.bed> ...

Output:
    <output_prefix>.robust.bed
	Column 1:  	chr of merged peak region
	Column 2:  	start of merged peak region (0-based)
	Column 3:  	end of merged peak region (0-based, end-excluded)
	Column 4:  	details of overlapping peak regions in each repeat
	Column 5:  	number of supporting repeats
	Column 6:  	strand is always *
	Column 7:  	best p-value
	Column 8:  	best summit start (1-based)
	Column 9:  	best summit end (1-based, end-included)
	Column 10: 	best summit height
	Column 11: 	average summit coordinate (1-based)
      Note: robust (among repeats) is defined by being supported by >50% (no equal) repeats
            best (among repeats) is defined by lowest p-value
    <output_prefix>.robust.bdg
    <output_prefix>.robust_summit.bdg
    <output_prefix>.robust_avg_summit.bdg
	score = number of supporting repeats
".$version;

if(@ARGV < 3){
	die $usage
}
my $output_prefix = shift(@ARGV);
my $input_multi_bed_filename = shift(@ARGV);
my @in_filenames = @ARGV;

my $N = scalar(@in_filenames);
print STDERR "Parsing $N repeats in command-line argumetns\n";

my @S = ([]) x $N;
my ($i, @Lines);
for($i=0; $i<$N; $i++){
	print STDERR "Now read in peak bed file for repeat ".($i+1)." ".$in_filenames[$i]."\n";
	open(IN, $in_filenames[$i]) or die "Error: cannot open peak bed file for repeat ".($i+1)." ".$in_filenames[$i]."\n";
	@Lines = <IN>;
	$S[$i] = [ map { s/[\r\n]+$//; [ split/\t/ ]; } @Lines ];
	close(IN);
}


my $output_robust_bed_filename = $output_prefix.".robust.bed";
my $output_robust_bdg_filename = $output_prefix.".robust.bdg";
my $output_robust_summit_bdg_filename = $output_prefix.".robust_summit.bdg";
my $output_robust_avg_summit_bdg_filename = $output_prefix.".robust_avg_summit.bdg";

open(OUT_BED, ">".$output_robust_bed_filename) or die "Error: cannot open robust.bed $output_robust_bed_filename for output\n";
open(OUT_BDG, ">".$output_robust_bdg_filename) or die "Error: cannot open robust.bdg $output_robust_bdg_filename for output\n";
open(SMT_BDG, ">".$output_robust_summit_bdg_filename) or die "Error: cannot open robust_summit.bdg $output_robust_summit_bdg_filename for output\n";
open(AVG_BDG, ">".$output_robust_avg_summit_bdg_filename) or die "Error: cannot open robust_avg_summit.bdg $output_robust_avg_summit_bdg_filename for output\n";

my (@F, @G, @H, $FF);
my ($rep, $idx);
my ($bpp, $bss, $bse, $bsh);   # best p-value, summit start, summit end, summit height
my ($sss, $nnn, $avg);   # sum of summit coordinate for later computing summit average
my ($pp, $ss, $se, $sh);   # p-value, summit start, summit end, summit height
print STDERR "Now read in multi.bed $input_multi_bed_filename\n";
open(IN, $input_multi_bed_filename) or die "Error: cannot open multi.bed $input_multi_bed_filename for input\n";
while(<IN>){
	s/[\r\n]+$//;
	@F = split/\t/;
	($bpp, $bss, $bse, $bsh) = (undef, undef, undef, undef);   # best p-value, summit start, summit end, summit height
	($sss, $nnn) = (0, 0);   # sum of summit coordinate for later computing summit average
	if($F[4] > $N/2){   # peak region called in >50% repeats
		@G = split(/,/, $F[3]);
		($rep, $idx) = (undef, undef);
		foreach (@G){
			@H = split/:/;
			if(@H==2){
				$rep = $H[0];
				$rep =~ s/^rep//;
				$rep--;   # 1-based to 0-based index
				$idx = $H[1];
				$idx--;   # 1-based to 0-based index
			}else{
				$idx = $H[0];
				$idx--;   # 1-based to 0-based index
			}
			if(defined($rep) && defined($idx)){
				$FF = $S[$rep]->[$idx];
				($pp, $ss, $se, $sh) = @{$FF}[6..9];   # p-value, summit start, summit end, summit height
				$nnn++;
				$sss += $ss + $se;
				if(!defined($bpp) || $pp < $bpp){
					($bpp, $bss, $bse, $bsh) = ($pp, $ss, $se, $sh);
				}
			}
		}

		push(@F, $bpp, $bss, $bse, $bsh, $sss/$nnn/2);
		print OUT_BED join("\t", @F)."\n";
		print OUT_BDG join("\t", @F[0..2,4])."\n";
		print SMT_BDG join("\t", $F[0], $F[7]-1, $F[8], $F[4])."\n";
		$avg = int($F[10]+0.5);
		print AVG_BDG join("\t", $F[0], $avg-1, $avg, $F[4])."\n";
	}
}
close(IN);

close(AVG_BDG);
close(SMT_BDG);
close(OUT_BDG);
close(OUT_BED);

0;

