#!/usr/bin/perl

use strict;
use warnings;


my $version = '
Version: 0.1.1 (2021-01-02)
Author: Adam Yongxin Ye @ BCH
';
my $usage = "Usage: $0 <sequence.fa>
	[should_skip_all_0 (default:1)]
Output: .tsv (10 columns)
	1) chr
	2) position (0-based)
	3) forward base
	4) forward heptamer score
	5) forward 12RSS score
	6) forward 23RSS score
	7) reverse base
	8) reverse heptamer score
	9) reverse 12RSS score
	10) reverse 23RSS score
Note:
	This script will scan for
		heptamer = CACAGTG
		spacer of 12 or 23 bp
		nonamer = ACAAAAACC
	The 12/23 RSS score is simply calculated as:
		match +2, transition +0, transversion +0
	Score of first 3bp of heptamer (CAC) are multiplied by 10, and minus 58, so that
		matching CAC only has score =  2
		matching CACA has score = 4 = 2 + 2
		matching CACAGTG has score = 10 = 2 + 4*2
		matching CACAGTG and nonamer = 28 = 10 + 9*2
		CAC + random bases has expected score = 8.5 = 2+(2/4)*(4+9) 
".$version;
#		CAC + random bases has expected score = 11.75 = 2+(2/4+1/4)*(4+9) 

if(@ARGV < 1){
	die $usage;
}
my ($input_fasta_filename) = @ARGV;
my $should_skip_all_0 = 1;
if(@ARGV > 1){
	$should_skip_all_0 = $ARGV[1];
}


sub max{
	my ($ans) = @_;
	foreach (@_){
		if($_ > $ans){
			$ans = $_;
		}
	}
	return $ans;
}


my $heptamer = "CACAGTG";
my @H = split(//, $heptamer);
my $nonamer = "ACAAAAACC";
my @N = split(//, $nonamer);

#my %mismatch_score_mat = (   # if we +1 for transition, +2 for transversion
#	"A" => {"A" => 0, "C" => 2, "G" => 1, "T" => 2.0},
#	"C" => {"A" => 2, "C" => 0, "G" => 2.0, "T" => 1},
#	"G" => {"A" => 1, "C" => 2.0, "G" => 0, "T" => 2},
#	"T" => {"A" => 2.0, "C" => 1, "G" => 2, "T" => 0},
#);
#my $missing_score = 2;

#my %match_score_mat = (   # if we +2 for match, +1 for transition
#	"A" => {"A" => 2, "C" => 0, "G" => 1, "T" => 0.0},
#	"C" => {"A" => 0, "C" => 2, "G" => 0.0, "T" => 1},
#	"G" => {"A" => 1, "C" => 0.0, "G" => 2, "T" => 0},
#	"T" => {"A" => 0.0, "C" => 1, "G" => 0, "T" => 2},
#);
my %match_score_mat = (   # if we +2 for match
	"A" => {"A" => 2, "C" => 0, "G" => 0, "T" => 0.0},
	"C" => {"A" => 0, "C" => 2, "G" => 0.0, "T" => 0},
	"G" => {"A" => 0, "C" => 0.0, "G" => 2, "T" => 0},
	"T" => {"A" => 0.0, "C" => 0, "G" => 0, "T" => 2},
);

sub judge_RSS{
	my ($seq) = @_;
#	$seq =~ s/\s+//g;
	my @S = split(//, $seq);
	my $i;
	my ($Hscore, $N12score, $N23score) = (0, 0, 0);
	for($i=0; $i<3; $i++){
		if($i < @S){
			if(exists($match_score_mat{$S[$i]})){
				$Hscore += $match_score_mat{$S[$i]}->{$H[$i]} * 10;
			}
#		if($S[$i] ne $H[$i]){
#			$Hscore--;
#			last;
		}
	}
#	$Hscore++;
#	if($Hscore >= 0){
	for($i=3; $i<@H; $i++){
#		if($i >= @S){
#			$Hscore += $missing_score * (@H-$i);
#			last;
#		}
		if($i < @S){
			if(exists($match_score_mat{$S[$i]})){
				$Hscore += $match_score_mat{$S[$i]}->{$H[$i]};
			}
		}
	}
	
	for($i=0; $i<@N; $i++){
		if($i+@H+12 >= @S){
#			$N12score += $missing_score;
		}else{
			if(exists($match_score_mat{$S[$i+@H+12]})){
				$N12score += $match_score_mat{$S[$i+@H+12]}->{$N[$i]};
			}
		}
		if($i+@H+23 >= @S){
#			$N23score += $missing_score;
		}else{
			if(exists($match_score_mat{$S[$i+@H+23]})){
				$N23score += $match_score_mat{$S[$i+@H+23]}->{$N[$i]};
			}
		}
	}
	
#	return ($Hscore, $N12score, $N23score);
	return ($Hscore-60+2, max(0, $Hscore+$N12score-60+2), max(0, $Hscore+$N23score-60+2));
}


my %rC_hash = (
	"A" => "T",
	"C" => "G",
	"G" => "C",
	"T" => "A",
);

sub reverseComplement{
	my ($seq) = @_;
	my @S = split(//, $seq);
	return join("", map { 
			if(exists($rC_hash{$_})){
				$rC_hash{$_};
			}else{
				$_;
			}
		} reverse @S);
}


sub scan_RSS_for_one_chr{
	my ($chr, $input_seq, $should_skip_all_0) = @_;
	#my $input_seq = join("", <STDIN>);
	#$input_seq =~ s/\s+//g;
	#print "original:" . $input . "\n";
	#print "reverseC:" . reverseComplement($input) . "\n";
	my $rC_seq = reverseComplement($input_seq);
	my $len = length($input_seq);

	my $i;
	my (@F, @R);
	for($i=0; $i<$len; $i++){
		@F = judge_RSS(substr($input_seq, $i, 7+9+23));
		@R = judge_RSS(substr($rC_seq, $len-1-$i, 7+9+23));
		if($should_skip_all_0){
			if($F[0]<=0 && $F[1]<=0 && $F[2]<=0 && $R[0]<=0 && $R[1]<=0 && $R[2]<=0){
				next;
			}
		}
		print join("\t", $chr, $i, 
				substr($input_seq, $i, 1), @F,
				substr($rC_seq, $len-1-$i, 1), @R,
			)."\n";

	}
}


my ($chr, $seq, $tmp) = ("", "", "");
open(IN, $input_fasta_filename) or die "Error: cannot open input fasta $input_fasta_filename\n";
while(<IN>){
	if(/^>(.*?)\s*$/){
		$tmp = $1;
		if($seq ne ""){
			print STDERR "Now scanning RSS for $chr ...\n";
			scan_RSS_for_one_chr($chr, $seq, $should_skip_all_0);
		}
		($chr, $seq) = ($tmp, "");
		print STDERR "Now reading chr = $chr ...\n";
	}else{
		s/\s+//g;
		tr/acgtn/ACGTN/;
		$seq .= $_;
	}
}
if($seq ne ""){
	print STDERR "Now scanning RSS for $chr ...\n";
	scan_RSS_for_one_chr($chr, $seq, $should_skip_all_0);
}

close(IN);


0;



