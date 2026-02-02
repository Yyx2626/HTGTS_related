#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw(max);

my $version = '
Version: 0.1.0 (2026-01-12)
Authors: Adam Yongxin Ye @ BCH & Gemini & ChatGPT
';
my $usage = "Usage: perl $0 <input.tlx> <ref.fa>
	[max_allowed_variants (default:2)]
	[check_prey_rel_start (default:3)] [check_prey_rel_end (default:22)]
	[filtered_out.tlx (default:undef)]

Description:
	Filters tlx records by evaluating the number of variants within a specific prey subregion.
	To count the variants, the prey sequence is re-aligned to tlx-reported reference region using the Needleman-Wunsch (global) algorithm.

Note: check_prey_rel_start, check_prey_rel_end: 1-based, end-included

Output: STDOUT (tlx format; read lines with too many mutations in the prey subregion should be filtered out)
".$version;

if(@ARGV < 2){
	die $usage;
}
my $input_tlx_filename = shift(@ARGV);
my $ref_fasta_filename = shift(@ARGV);
my $max_allowed_variants = 2;
if(@ARGV > 0){
	$max_allowed_variants = shift(@ARGV);
}
my $check_prey_rel_start = 3;
if(@ARGV > 0){
	$check_prey_rel_start = shift(@ARGV);
}
my $check_prey_rel_end = 22;
if(@ARGV > 0){
	$check_prey_rel_end = shift(@ARGV);
}
my $filtered_out_filename = undef;
if(@ARGV > 0){
	$filtered_out_filename = shift(@ARGV);
}

if(defined($filtered_out_filename)){
	open(FILTER_OUT, ">".$filtered_out_filename) or die "Error: cannot open $filtered_out_filename for output\n";
}

my $ref_hash = read_fasta_to_hash($ref_fasta_filename);

### ref: yyx_show_tlx_alignment.20210701.pl
my %f2i = ();
$f2i{"Qname"} = 0;
$f2i{"Rname"} = 2;
$f2i{"Junction"} = 3;
$f2i{"Strand"} = 4;
$f2i{"Rstart"} = 5;
$f2i{"Rend"} = 6;
$f2i{"B_Rname"} = 7;
$f2i{"B_Rstart"} = 8;
$f2i{"B_Rend"} = 9;
$f2i{"B_Strand"} = 10;
$f2i{"B_Qstart"} = 11;
$f2i{"B_Qend"} = 12;
$f2i{"Qstart"} = 13;
$f2i{"Qend"} = 14;
$f2i{"Seq"} = 18;

my (@F, $i);
my ($readSeq, $refPreySeq, $readPreySeq);
my ($chr, $start, $end, $strand, $bedStart, $bedEnd);
my ($readStart, $readEnd);
open(IN, $input_tlx_filename) or die "Error: cannot open tlx file $input_tlx_filename for input\n";
print STDERR "Now parse tlx file $input_tlx_filename ...\n";
while(<IN>){
	s/[\r\n]+$//;
	@F = split/\t/;
	if($F[0] eq "Qname"){
		%f2i = ();
		for($i=0; $i<@F; $i++){
			$f2i{$F[$i]} = $i;
		}
		print join("\t", @F, qw/prey_realign_check_Cigar/)."\n";
		if(defined($filtered_out_filename)){
			print FILTER_OUT join("\t", @F, qw/prey_realign_check_Cigar/)."\n";
		}
		next;
	}
	
	$readSeq = $F[$f2i{"Seq"}];
	$readStart = $F[$f2i{"Qstart"}] - 1;
	$readEnd = $F[$f2i{"Qend"}];
	
	$chr = $F[$f2i{"Rname"}];
	$start = $F[$f2i{"Rstart"}] - 1;
	$end = $F[$f2i{"Rend"}];
	$strand = $F[$f2i{"Strand"}] > 0 ? "+" : "-";
	if(!exists($ref_hash->{$chr})){
		print STDERR "Warning: chr $chr not exist in ref.fa\n";
	}
	$refPreySeq = substr($ref_hash->{$chr}, $start, $end-$start);
	if($strand eq "-"){
		$refPreySeq = reverse_complement($refPreySeq);
	}
	
	$readPreySeq = substr($readSeq, $readStart, $readEnd-$readStart);
	$refPreySeq =~ tr/a-z/A-Z/;
	$readPreySeq =~ tr/a-z/A-Z/;
	my ($aln_query, $aln_subject) = needleman_wunsch($readPreySeq, $refPreySeq);
	my $cigar = get_cigar_query_range($aln_query, $aln_subject, $check_prey_rel_start, $check_prey_rel_end);
	my $cigar_count_hash = parse_cigar_to_total_count_hash($cigar);
	my $total = 0;
	my $match = 0;
	foreach (keys %$cigar_count_hash){
		$total += $cigar_count_hash->{$_};
		$match += $cigar_count_hash->{$_} if $_ eq "=" or $_ eq "M";
	}
	if($total - $match <= $max_allowed_variants){
		print join("\t", @F, $cigar)."\n";
	}else{
		if(defined($filtered_out_filename)){
			print FILTER_OUT join("\t", @F, $cigar)."\n";
		}
	}

}
close(IN);
if(defined($filtered_out_filename)){
	close(FILTER_OUT);
}

0;


sub read_fasta_to_hash{
	my ($fasta_file) = @_;
	my %seq_hash;

	open(FASTA, '<', $fasta_file) or die "Error: cannot open fasta file $fasta_file for input: $!\n";
	print STDERR "Now read in fasta file $fasta_file ...\n";
	my $current_id;
	while (<FASTA>) {
		next if /^\s*$/;   # Skip empty lines
		s/[\r\n]+$//;
		if (/^>(\S+)/) {
			# New FASTA header
			$current_id = $1;
			$seq_hash{$current_id} = "";
		} else {
			# Sequence line
			next unless defined $current_id;  # safety check
			$seq_hash{$current_id} .= $_;
		}
	}
	close(FASTA);
	return \%seq_hash;
}

sub reverse_complement{
	my @ans = map { tr/ACGTacgt/TGCAtgca/; join("", reverse split//); } @_;
	if(@_ > 1){
		return @ans;
	}else{
		return $ans[0];
	}
}

# --- Subroutine: Needleman-Wunsch Global Alignment ---
sub needleman_wunsch {
	my ($seq_q, $seq_s) = @_;
	my @q = split("", $seq_q);
	my @s = split("", $seq_s);
	
	my $match = 1;
	my $mismatch = -1;
	my $gap = -1;
	
	my @score;
	$score[0][0] = 0;
	
	# Initialization
	for (my $i = 1; $i <= @q; $i++) { $score[$i][0] = $i * $gap; }
	for (my $j = 1; $j <= @s; $j++) { $score[0][$j] = $j * $gap; }
	
	# Fill Matrix
	for (my $i = 1; $i <= @q; $i++) {
		for (my $j = 1; $j <= @s; $j++) {
			my $s_score = ($q[$i-1] eq $s[$j-1]) ? $match : $mismatch;
			my $diag = $score[$i-1][$j-1] + $s_score;
			my $up   = $score[$i-1][$j] + $gap;
			my $left = $score[$i][$j-1] + $gap;
			$score[$i][$j] = max($diag, $up, $left);
		}
	}
	
	# Traceback
	my ($res_q, $res_s) = ("", "");
	my ($i, $j) = (scalar @q, scalar @s);
	
	while ($i > 0 || $j > 0) {
		if ($i > 0 && $j > 0 && $score[$i][$j] == ($score[$i-1][$j-1] + ($q[$i-1] eq $s[$j-1] ? $match : $mismatch))) {
			$res_q = $q[$i-1] . $res_q;
			$res_s = $s[$j-1] . $res_s;
			$i--; $j--;
		} elsif ($i > 0 && $score[$i][$j] == $score[$i-1][$j] + $gap) {
			$res_q = $q[$i-1] . $res_q;
			$res_s = "-" . $res_s;
			$i--;
		} else {
			$res_q = "-" . $res_q;
			$res_s = $s[$j-1] . $res_s;
			$j--;
		}
	}
	return ($res_q, $res_s);
}

# --- Subroutine: Extract CIGAR based on Query Coordinates ---
sub get_cigar_query_range {
    my ($aln_q, $aln_s, $q_start, $q_end) = @_;   # start, end: 1-based, end-included, on query
    my @q_chars = split("", $aln_q);
    my @s_chars = split("", $aln_s);
    
    my $current_query_pos = 0;
    my @ops;
    
    for (my $i = 0; $i < @q_chars; $i++) {
        # Increment coordinate only if it's NOT a gap in the query
        if ($q_chars[$i] ne "-") {
            $current_query_pos++;
        }
        
        # Check if current alignment column falls within the query coordinate range
        if ($current_query_pos >= $q_start && $current_query_pos <= $q_end) {
            if ($q_chars[$i] eq "-") {
                push @ops, "D"; # Deletion: Gap in Query, base in Subject
            } elsif ($s_chars[$i] eq "-") {
                push @ops, "I"; # Insertion: Base in Query, gap in Subject
            } elsif ($q_chars[$i] eq "N" || $s_chars[$i] eq "N") {
                # If either is 'N', use 'M' to represent a generic alignment match
                push @ops, "M";
            } elsif ($q_chars[$i] eq $s_chars[$i]) {
                push @ops, "="; # Sequence Match
            } else {
                push @ops, "X"; # Sequence Mismatch
            }
        }
        last if $current_query_pos > $q_end;
    }
    
    # Summarize into CIGAR string
    my $cigar_str = "";
    if (@ops) {
        my $count = 1;
        for (my $i = 1; $i <= @ops; $i++) {
            if ($i < @ops && $ops[$i] eq $ops[$i-1]) {
                $count++;
            } else {
                $cigar_str .= $count . $ops[$i-1];
                $count = 1;
            }
        }
    }
    
    return $cigar_str;
}

# --- Subroutine: Parse any valid CIGAR operator ---
sub parse_cigar_to_total_count_hash {
	my ($cigar) = @_;
	my %counts;

	# Standard CIGAR operators: MIDNSHP=X
	# \d+ matches the length
	# ([MIDNSHP=X]) captures any one of the standard operators
	while ($cigar =~ /(\d+)([MIDNSHP=X])/g) {
		my $length = $1;
		my $type = $2;
		
		$counts{$type} += $length;
	}

	return \%counts;
}

