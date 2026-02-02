#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $version = '
Version: 0.1.5 (2023-02-16)
Author: Adam Yongxin Ye @ BCH
';
my $usage = "Usage: $0 <input.tlx> <ref.fa> <VDJ.bed> <V.fa> <J.fa> <J.aux>
	[max_match_distance (default: 10)] [MH_shorten_prey_or_bait (default: bait)]
Input:
	<J.aux>  the optional file in IgBLAST to annotate J frame, 4 columns:
		gene name  (e.g. JH1 , IGHJ1*01)
		first coding frame start position (0-based)  (e.g. 0 , 1 , 2)
		chain type  (e.g. JH , JK , JL ; I will ignore this column)
		CDR3 end (e.g. 18 , 13 , 6 , 7)
Output: STDOUT   append several columns as follows
	MH_len   ( >0: microhomology length, <0: mid length  for each read in <input.tlx>)
	pre  bait  mid  prey  post   (segmented sequences on each read in <input.tlx>)
	Bfeature  Pfeature   (annotation of overlapping features in <VDJ.bed>)
	Vpart  midO  Jpart   (extended V(D)J sequence)
	InFrame  Stop  Productive
".$version;

if(@ARGV < 6){
	die $usage;
}

my $input_filename = shift(@ARGV);
my $ref_fa_filename = shift(@ARGV);
my $VDJ_bed_filename = shift(@ARGV);
my $V_fa_filename = shift(@ARGV);
my $J_fa_filename = shift(@ARGV);
my $J_aux_filename = shift(@ARGV);

### 2023-02-15, add two more optional arguments as below
my $max_match_distance = 10;
if(@ARGV > 0){
	$max_match_distance = shift(@ARGV);
}

my $MH_shorten_prey_or_bait = "bait";
if(@ARGV > 0){
	$MH_shorten_prey_or_bait = shift(@ARGV);
}

print STDERR "[ARG] input_filename=$input_filename\n";
print STDERR "[ARG] ref_fa_filename=$ref_fa_filename\n";
print STDERR "[ARG] VDJ_bed_filename=$VDJ_bed_filename\n";
print STDERR "[ARG] V_fa_filename=$V_fa_filename\n";
print STDERR "[ARG] J_fa_filename=$J_fa_filename\n";
print STDERR "[ARG] J_aux_filename=$J_aux_filename\n";
print STDERR "[ARG] max_match_distance=$max_match_distance\n";
print STDERR "[ARG] MH_shorten_prey_or_bait=$MH_shorten_prey_or_bait\n";

if($MH_shorten_prey_or_bait ne "bait" && $MH_shorten_prey_or_bait ne "prey"){
	die "Error: Unrecognized MH_shorten_prey_or_bait=$MH_shorten_prey_or_bait, which should be either 'bait' or 'prey'\n";
}


sub read_fasta{
	my ($input_fa_filename) = @_;
	my $title = undef;
	my $seq = "";
	my $ans = {};
	open(IN, $input_fa_filename) or die "Error: cannot open fasta file $input_fa_filename for input\n";
	print STDERR "Now reading fasta file $input_fa_filename ...\n";
	while(<IN>){
		s/[\r\n]+$//;
		if(/^>(.*?)\s*$/){
			if(defined($title)){
				if(exists($ans->{$title})){
					print STDERR "Warning: duplicated sequence title $title in fasta file $input_fa_filename\n";
				}
				$ans->{$title} = $seq;
			}
			$title = $1;
			$seq = "";
		}else{
			$seq .= $_;
		}
	}
	if(defined($title)){
		if(exists($ans->{$title})){
			print STDERR "Warning: duplicated sequence title $title in fasta file $input_fa_filename\n";
		}
		$ans->{$title} = $seq;
	}
	close(IN);
	return $ans;
}

my $ref_hash = read_fasta($ref_fa_filename);
my $V_hash = read_fasta($V_fa_filename);
my $J_hash = read_fasta($J_fa_filename);
my ($tmp_hash, $value);
$tmp_hash = $V_hash;
$V_hash = {};
foreach (keys %$tmp_hash){
	$value = $tmp_hash->{$_};
	s/^lcl[|]//;
	$V_hash->{$_} = $value;
}
$tmp_hash = $J_hash;
$J_hash = {};
foreach (keys %$tmp_hash){
	$value = $tmp_hash->{$_};
	s/^lcl[|]//;
	$J_hash->{$_} = $value;
}
my $VJ_hash = {};
foreach (keys %$V_hash){
	$VJ_hash->{$_} = $V_hash->{$_};
}
foreach (keys %$J_hash){
	$VJ_hash->{$_} = $J_hash->{$_};
}


sub read_aux_file{
	my ($input_aux_filename) = @_;
	my $Jaux_hash = {};
	my @F;
	open(IN, $input_aux_filename) or die "Error: cannot open aux file $input_aux_filename for input\n";
	print STDERR "Now reading aux file $input_aux_filename ...\n";
	while(<IN>){
		next if(/^#/);
		s/[\r\n]+$//;
		@F = split/\s+/;
		$Jaux_hash->{$F[0]} = [@F[1..3]];
	}
	close(IN);
	return $Jaux_hash;
}

my $Jaux_hash = read_aux_file($J_aux_filename);
#foreach (sort keys %$Jaux_hash){
#	print STDERR "[DEBUG] Jaux $_\n";
#}
foreach (sort keys %$J_hash){
	if(!exists($Jaux_hash->{$_})){
		print STDERR "Warning: J $_ does not exist in the aux file\n";
	}
}



my @four_bases = qw/ A C G T /;
my $i;
my %complement_base_hash = ();
my ($b1, $b2);
for($i=0; $i<=2; $i++){
	$b1 = $four_bases[$i];
	$b2 = $four_bases[3-$i];
	if($i==2){
		$b1 = "N";
		$b2 = "N";
	}
	$complement_base_hash{$b1} = $b2;
	$complement_base_hash{$b2} = $b1;
	$b1 =~ tr/A-Z/a-z/;
	$b2 =~ tr/A-Z/a-z/;
	$complement_base_hash{$b1} = $b2;
	$complement_base_hash{$b2} = $b1;
}
sub hash_get_default{
	my ($hash, $query, $default, $error_message) = @_;
	if(exists($hash->{$query})){
		return $hash->{$query};
	}
	if(defined($error_message)){
		print STDERR $error_message;
	}
	return $default;
}
sub reverse_complement{
	my @S;
	my @ans = map{
		join("", map{
			hash_get_default(\%complement_base_hash, $_, $_);
		} reverse split// );
	} @_;
	if(@ans==1){
		return $ans[0];
	}
	return @ans;
}
#print join("\n", reverse_complement("ACCGT", "AAAA"))."\n";



### Editing distance
### ref: https://www.geeksforgeeks.org/edit-distance-dp-5/
sub min{
	my $ans = undef;
	foreach (@_){
		if(!defined($ans)){
			$ans = $_;
		}elsif($_ < $ans){
			$ans = $_;
		}
	}
	return $ans;
}
sub editDistDP{
	my ($str1, $str2) = @_;
	my @S1 = split(//, $str1);
	my @S2 = split(//, $str2);
	my $m = @S1;
	my $n = @S2;
	my @dp = ();
	my ($i, $j);
	for($i=0; $i<=$m; $i++){
		push(@dp, [(0) x ($n+1)]);
	}
	for($i=0; $i<=$m; $i++){
		for($j=0; $j<=$n; $j++){
			if($i==0){
				$dp[$i]->[$j] = $j;
			}elsif($j==0){
				$dp[$i]->[$j] = $i;
			}elsif($S1[$i-1] eq $S2[$j-1]){
				$dp[$i]->[$j] = $dp[$i-1]->[$j-1];
			}else{
				$dp[$i]->[$j] = 1 + min($dp[$i]->[$j-1], $dp[$i-1]->[$j], $dp[$i-1]->[$j-1]);
			}
		}
	}
	return $dp[$m][$n];
}



sub read_bed{
	my ($input_bed_filename) = @_;
	my %H = ();
	my @S = ();
	my $NR = 0;
	my (@F, $chr, $start, $end, $name, $strand);
	open(IN, $input_bed_filename) or die "Error: cannot open bed file $input_bed_filename for input\n";
	print STDERR "Now reading bed file $input_bed_filename ...\n";
	while(<IN>){
		$NR++;
		s/[\r\n]+$//;
		@F = split/\t/;
#		push(@S, [@F]);
#		if($NR != (@S - 1)){
#			die "Error: strange unalignment between NR=$NR and \@S in read_bed\n";
#		}
		if(@F < 4){
			print STDERR "Warning: Line $NR in bed file ($_) does not have >= 4 columns, so I skip it\n";
			next;
		}
		($chr, $start, $end, $name) = @F[0..3];
		$strand = "*";
		if(@F >= 6){
			$strand = $F[5];
			if($strand ne "*" && $strand ne "+" && $strand ne "-"){
				print STDERR "Warning: Unrecognized strand $strand at Line $NR in bed file, so I set strand to *\n";
				$strand = "*";
			}
		}
#		$name =~ s/ [(][^()]*[)]$//g;
		if(!exists($H{$name})){
			$H{$name} = [$chr, $start, $end, $strand];
			push(@S, $name);
		}else{
			print STDERR "Warning: duplicated records for $name in bed file\n";
		}
#		$NR++;
	}
	close(IN);
	return (\@S, \%H);
}
my ($VDJ_S, $VDJ_H) = read_bed($VDJ_bed_filename);


### guess bed format (start 0/1-based, end-excluded/included, and strand ?)
sub best_match_and_guess_bed_format{
	my ($chr, $start, $end, $name, $chr_seq, $VJ_hash) = @_;
	my @VJ_name_vec = ();
	my ($name2, $name3, $name4);
	my @VJ_names = keys %$VJ_hash;
	$name2 = $name;
	if(exists($VJ_hash->{$name2})){
#		print STDERR "[DEBUG] find exact $name\n";
		push(@VJ_name_vec, $name2);
	}else{
		foreach $name3 (@VJ_names){
			$name4 = $name3;
			$name4 =~ s/[ (].*$//;
			$name4 =~ tr/a-z/A-Z/;
			$name4 =~ s/[*].*$//;
			$name4 =~ s/P$//;
			if($name2 eq $name4){
#				print STDERR "[DEBUG] find $name3 for $name\n";
				push(@VJ_name_vec, $name3);
			}
		}
	}
	if(@VJ_name_vec == 0){
		$name2 =~ s/[ (].*$//;
		$name2 =~ tr/a-z/A-Z/;
		$name2 =~ s/P$//;
		if(exists($VJ_hash->{$name2})){
#			print STDERR "[DEBUG] find $name2 for $name\n";
			push(@VJ_name_vec, $name2);
		}else{
			foreach $name3 (@VJ_names){
				$name4 = $name3;
				$name4 =~ s/[ (].*$//;
				$name4 =~ tr/a-z/A-Z/;
				$name4 =~ s/[*].*$//;
				$name4 =~ s/P$//;
				if($name2 eq $name4){
#					print STDERR "[DEBUG] find $name3 for $name\n";
					push(@VJ_name_vec, $name3);
				}
			}
		}
	}
	if(@VJ_name_vec == 0){
		$name2 =~ s/^IG(.)(.)/$2$1/;
		if(exists($VJ_hash->{$name2})){
#			print STDERR "[DEBUG] find $name2 for $name\n";
			push(@VJ_name_vec, $name2);
		}else{
			foreach $name3 (@VJ_names){
				$name4 = $name3;
				$name4 =~ s/[ (].*$//;
				$name4 =~ tr/a-z/A-Z/;
				$name4 =~ s/[*].*$//;
				$name4 =~ s/P$//;
				if($name2 eq $name4){
#					print STDERR "[DEBUG] find $name3 for $name\n";
					push(@VJ_name_vec, $name3);
				}
			}
		}
	}
	if(@VJ_name_vec == 0){
		$name2 =~ s/^([VJ])([HKL])/IG$2$1/;
		if(exists($VJ_hash->{$name2})){
#			print STDERR "[DEBUG] find $name2 for $name\n";
			push(@VJ_name_vec, $name2);
		}else{
			foreach $name3 (@VJ_names){
				$name4 = $name3;
				$name4 =~ s/[ (].*$//;
				$name4 =~ tr/a-z/A-Z/;
				$name4 =~ s/[*].*$//;
				$name4 =~ s/P$//;
				if($name2 eq $name4){
#					print STDERR "[DEBUG] find $name3 for $name\n";
					push(@VJ_name_vec, $name3);
				}
			}
		}
	}
	if(@VJ_name_vec == 0){
		print STDERR "Warning: cannot find name $name in V.fa nor J.fa, so I will check for all\n";
		push(@VJ_name_vec, sort keys %$VJ_hash);
	}
	my $best_VJ_name = undef;
	my $best_key = undef;
	my $best_ED = undef;
	my ($VJ_name, $VJ_seq, %target_seqs, $target_len, $has_exact_match, $start_base, $end_included, $strand);
	my ($now_best_key, $now_best_ED, $ref_seq, $key, $ED);
	foreach $VJ_name (@VJ_name_vec){
		$VJ_seq = $VJ_hash->{$VJ_name};
		$VJ_seq =~ tr/a-z/A-Z/;
		%target_seqs = (
			"+" => $VJ_seq,
			"-" => reverse_complement($VJ_seq)
		);
		$target_len = length($VJ_seq);
		
		$has_exact_match = 0;
		$now_best_key = undef;
		$now_best_ED = undef;
		foreach $end_included (0..1){
			next if(($end-$start+$end_included)!=$target_len);
			foreach $start_base (0..1){
				$ref_seq = substr($chr_seq, $start-$start_base, $target_len);
				$ref_seq =~ tr/a-z/A-Z/;
				foreach $strand (qw/+ -/){
					$key = $start_base.$end_included.$strand;
					$ED = editDistDP($ref_seq, $target_seqs{$strand});
					if(!defined($now_best_ED) || $ED < $now_best_ED){
						$now_best_ED = $ED;
						$now_best_key = $key;
					}
					if($ED==0){
						$has_exact_match = 1;
					}
					last if($has_exact_match);
				}
				last if($has_exact_match);
			}
			last if($has_exact_match);
		}
		if(defined($now_best_ED) && (!defined($best_ED) || $now_best_ED < $best_ED)){
			$best_ED = $now_best_ED;
			$best_key = $now_best_key;
			$best_VJ_name = $VJ_name;
			last if($has_exact_match);
		}
	}
#	print STDERR "[DEBUG] for $name, best match is $best_VJ_name, with distance $best_ED, format $best_key\n";
	return ($best_key, $best_VJ_name, $best_ED);
}

my ($chr_seq, $best_key, $best_name, $best_ED);
my ($chr, $start, $end, $strand, $name, $ref_seq);
my %bed_coord_guess_hash = ();
my $NR = 0;
foreach $name (@$VDJ_S){
	$NR++;
	($chr, $start, $end, $strand) = @{$VDJ_H->{$name}};
	$chr_seq = hash_get_default($ref_hash, $chr, "-", "Warning: ($NR/".scalar(@$VDJ_S).") for $name, chr $chr does not exist in ref.fa\n");
	if($chr_seq ne "-"){
		if($name =~ /[VJ]/i){
			($best_key, $best_name, $best_ED) = best_match_and_guess_bed_format($chr, $start, $end, $name, $chr_seq, $VJ_hash);
			if(defined($best_ED)){
				if($best_ED <= 0){
					push(@{$VDJ_H->{$name}}, $best_key, $best_name);
					$bed_coord_guess_hash{$best_key}++;
					print STDERR "[DEBUG] ($NR/".scalar(@$VDJ_S).")  with format $best_key ,  $name  can be exactly matched to  $best_name\n";
				}elsif($best_ED <= $max_match_distance){   ### 2023-02-15, allow roughly matching
					push(@{$VDJ_H->{$name}}, $best_key, $best_name);
					$bed_coord_guess_hash{$best_key}++;
					print STDERR "[DEBUG] ($NR/".scalar(@$VDJ_S).")  with format $best_key ,  $name  can be roughly matched to  $best_name with distance $best_ED\n";
				}else{
					print STDERR "Warning: ($NR/".scalar(@$VDJ_S).") cannot match $name in V.fa nor J.fa, best match to $best_name with distance $best_ED\n";
				}
			}else{
				print STDERR "Warning: ($NR/".scalar(@$VDJ_S).") cannot match $name in V.fa nor J.fa, no best match\n";
			}
		}
	}
}
if(scalar(keys %bed_coord_guess_hash)==0){
	die "Error: no one in VDJ.bed can be matched in V.fa nor J.fa, so I cannot guess the coordincate format of VDJ.bed.\n";
}
my ($key);
my %guess_by_strand = ();   ### 2023-02-15, change to guess by strand
my ($now_start_base, $now_end_included, $now_strand);
#$best_key = undef;
foreach $key (sort { $bed_coord_guess_hash{$b} <=> $bed_coord_guess_hash{$a} } keys %bed_coord_guess_hash){
	print STDERR join("\t", "[DEBUG] bed_format", $key, $bed_coord_guess_hash{$key})."\n";
	($now_start_base, $now_end_included, $now_strand) = split(//, $key);
	if(!exists($guess_by_strand{$now_strand})){
		$guess_by_strand{$now_strand} = [$now_start_base, $now_end_included, $now_strand];
	}
	if(!exists($guess_by_strand{"*"})){
		$guess_by_strand{"*"} = [$now_start_base, $now_end_included, $now_strand];
	}
}
#my ($guess_start_base, $guess_end_included, $guess_strand) = split(//, $best_key);
#die;



sub get_ref_seq{
	my ($chr, $start, $end, $strand) = @_;
	my $ref = "-";
	$ref = hash_get_default($ref_hash, $chr, "-", "Warning: $chr does not exist in reference fasta\n");
	if($ref ne "-"){
		$ref = substr($ref, $start-1, $end-$start+1);
	}
	if($strand eq "-"){
		$ref = reverse_complement($ref);
	}
	$ref =~ tr/a-z/A-Z/;
	return $ref;
}

sub calculate_overlap_len{
	my ($s1, $e1, $s2, $e2) = @_;
	## assume sorted $s1 <= $e2, $s2 <= $e2
	if($s1 > $e2 or $s2 > $e1){
		if($s1==$e2+1 or $s2==$e1+1){
			return 0.1;
		}
		return 0;
	}else{
		if($s1 <= $s2 and $s2 <= $e1 and $e1 <= $e2){
			return ($e1-$s2+1);
		}elsif($s2 <= $s1 and $s1 <= $e2 and $e2 <= $e1){
			return ($e2-$s1+1);
		}elsif($s1 <= $s2 and $e2 <= $e1){
			return ($e2-$s2+1);
		}elsif($s2 <= $s1 and $e1 <= $e2){
			return ($e1-$s1+1);
		}else{
			return -1;
		}
	}
}

sub get_or_guess_format{
	my ($chr, $start, $end, $bed_strand, $format_key) = @_;
	my ($start_base, $end_included, $strand);
	my $is_format_guessed = 0;
	if(defined($format_key)){
		($start_base, $end_included, $strand) = split(//, $format_key);
	}else{
		($start_base, $end_included, $strand) = @{$guess_by_strand{"*"}};
		if(exists($guess_by_strand{$bed_strand})){
			($start_base, $end_included, $strand) = @{$guess_by_strand{$bed_strand}};
		}
		$is_format_guessed = 1;
	}
	return ($start_base, $end_included, $strand, $is_format_guessed);
}
sub adjust_start_end{
	my ($chr, $start, $end, $bed_strand, $format_key) = @_;
	my ($start_base, $end_included, $strand, $is_format_guessed) = get_or_guess_format($chr, $start, $end, $bed_strand, $format_key);
	## convert to 1-based, end-included
	$start += 1 - $start_base;
	$end += 1 - $start_base;
	$end -= 1 - $end_included;
	return ($chr, $start, $end, $strand, $is_format_guessed);
}

foreach $name (@$VDJ_S){
	my ($start_base, $end_included, $strand, $is_format_guessed) = get_or_guess_format(@{$VDJ_H->{$name}});
	if($is_format_guessed){
		print STDERR "Warning: I guessed format " . join("", $start_base, $end_included, $strand) . " for $name\n";
	}
}

sub get_overlapping_features{
	### judge overlapping features (low efficiency)
	my ($Qchr, $Qstart, $Qend, $Qstrand, $Qjunction) = @_;
	my ($chr, $start, $end, $name, $format_key, $start_base, $end_included, $strand, $is_format_guessed);
	my ($overlap_len, $Jdistance, $Jstart, $Jend);
	my @ans = ();
	foreach $name (@$VDJ_S){
#		($chr, $start, $end, $format_key) = @{$VDJ_H->{$name}};
		($chr, $start, $end, $strand, $is_format_guessed) = adjust_start_end(@{$VDJ_H->{$name}});
		next if($chr ne $Qchr);
		$overlap_len = calculate_overlap_len($Qstart, $Qend, $start, $end);
		if($overlap_len >= 1){
			$Jdistance = 0;
			$Jstart = $Qjunction - $start;
			$Jend = $Qjunction - $end;
			if($Jstart * $Jend > 0){
				$Jdistance = min(abs($Jstart), abs($Jend));
			}
			push(@ans, [$Jdistance, $name, $Qstrand eq $strand]);
		}
	}
	return sort { $a->[0] <=> $b->[0] } @ans;
}

sub judge_VDJname_category{
	my ($name) = @_;
	my $name2 = $VDJ_H->{$name}->[5];
	if(defined($name2)){
		if(exists($V_hash->{$name2})){
			return "V";
		}elsif(exists($J_hash->{$name2})){
			return "J";
		}
	}else{
		if($name =~ /^V/i || $name =~ /^IG.V/i){
			return "V";
		}elsif($name =~ /^J/i || $name =~ /^IG.J/i){
			return "J";
		}
	}
	return "other";
}

sub str_left{
	my ($str, $len) = @_;
	return substr($str, 0, $len);
}
sub str_right{
	my ($str, $len) = @_;
	return substr($str, length($str)-$len, $len);
}
sub str_trim_left{
	my ($str, $len) = @_;
	return substr($str, $len);
}
sub str_trim_right{
	my ($str, $len) = @_;
	return substr($str, 0, length($str)-$len);
}

sub cigar2align{
	my ($Qstart, $Qend, $Rstart, $Rend, $Strand, $Cigar) = @_;
	my ($Qidxes, $Ridxes, $ops) = ([], [], []);
	my ($nowQpos, $nowRpos, $now_len, $now_op);
	if($Strand eq "+"){
		$nowQpos = $Qstart;
		$nowRpos = $Rstart;
		### ref: https://samtools.github.io/hts-specs/SAMv1.pdf
		while($Cigar =~ /([0-9]+)([MIDNSHPX=])/g){
			$now_len = $1;
			$now_op = $2;
			push(@$ops, ($now_op) x $now_len);
#			if($now_op eq "M" || $now_op eq "I" || $now_op eq "S" || $now_op eq "=" || $now_op eq "X"){
			## The CIGAR op N for strand + in tlx file seems to be misleading and lead to discordance on Qend
			if($now_op eq "M" || $now_op eq "I" || $now_op eq "N" || $now_op eq "S" || $now_op eq "=" || $now_op eq "X"){
				push(@$Qidxes, ($nowQpos .. ($nowQpos+$now_len-1)));
				$nowQpos += $now_len;
			}else{
				push(@$Qidxes, ($nowQpos - 0.5) x $now_len);
			}
			if($now_op eq "M" || $now_op eq "D" || $now_op eq "N" || $now_op eq "=" || $now_op eq "X"){
#			if($now_op eq "M" || $now_op eq "D" || $now_op eq "=" || $now_op eq "X"){
				push(@$Ridxes, ($nowRpos .. ($nowRpos+$now_len-1)));
				$nowRpos += $now_len;
			}else{
				push(@$Ridxes, ($nowRpos - 0.5) x $now_len);
			}
		}
		if($nowQpos != $Qend + 1){
			print STDERR "Warning: + nowQpos($nowQpos) != Qend($Qend)+1 for CIGAR=$Cigar at Line $NR\n";
		}
		if($nowRpos != $Rend + 1){
			print STDERR "Warning: + nowRpos($nowRpos) != Rend($Rend)+1 for CIGAR=$Cigar at Line $NR\n";
		}
	}else{
		$nowQpos = $Qstart;
		$nowRpos = $Rend;
		### ref: https://samtools.github.io/hts-specs/SAMv1.pdf
		while($Cigar =~ /([0-9]+)([MIDNSHPX=])/g){
			$now_len = $1;
			$now_op = $2;
			push(@$ops, ($now_op) x $now_len);
#			if($now_op eq "M" || $now_op eq "I" || $now_op eq "S" || $now_op eq "=" || $now_op eq "X"){
			## The CIGAR op N for strand - in tlx file seems to be misleading and lead to discordance on Qend
			if($now_op eq "M" || $now_op eq "I" || $now_op eq "N" || $now_op eq "S" || $now_op eq "=" || $now_op eq "X"){
				push(@$Qidxes, ($nowQpos .. ($nowQpos+$now_len-1)));
				$nowQpos += $now_len;
			}else{
				push(@$Qidxes, ($nowQpos - 0.5) x $now_len);
			}
			## The CIGAR op N for strand - in tlx file seems to work
			if($now_op eq "M" || $now_op eq "D" || $now_op eq "N" || $now_op eq "=" || $now_op eq "X"){
				push(@$Ridxes, reverse (($nowRpos-$now_len+1) .. $nowRpos));
				$nowRpos -= $now_len;
			}else{
				push(@$Ridxes, ($nowRpos + 0.5) x $now_len);
			}
		}
		if($nowQpos != $Qend + 1){
			print STDERR "Warning: - nowQpos($nowQpos) != Qend($Qend)+1 for CIGAR=$Cigar at Line $NR\n";
		}
		if($nowRpos != $Rstart - 1){
			print STDERR "Warning: - nowRpos($nowRpos) != Rstart($Rstart)-1 for CIGAR=$Cigar at Line $NR\n";
		}
	}
	return [$Qidxes, $Ridxes, $ops];
}
sub align2cigar{
	my ($align) = @_;
	my @ops = @{$align->[2]};
	my $ans = "";
	my $lastOp = undef;
	my $lastLen = 0;
	my ($i, $nowOp);
	for($i=0; $i<@ops; $i++){
		$nowOp = $ops[$i];
		if(defined($lastOp) && ($nowOp ne $lastOp)){
			$ans += $lastLen . $lastOp;
			$lastLen = 0;
		}
		$lastOp = $nowOp;
		$lastLen++;
	}
	if(defined($lastOp)){
		$ans += $lastLen . $lastOp;
	}
	return $ans;
}
sub reverse_align{
	my ($align) = @_;
	return [ [reverse @{$align->[0]}], [reverse @{$align->[1]}], [reverse @{$align->[2]}] ];
}
sub which{
	my @bool_vec = @_;
	my @ans = ();
	my $i;
	for($i=0; $i<@bool_vec; $i++){
		if($bool_vec[$i]){
			push(@ans, $i);
		}
	}
	return @ans;
}
sub set_align_range{
	my ($align, $idx, $start, $end) = @_;
	my $alignLen = scalar(@{$align->[$idx]});
#	print STDERR "[DEBUG] set_align_range() alignLen=$alignLen\n";
	my @idxes = which( map { (!defined($start) || $_ >= $start) && (!defined($end) || $_ <= $end) } @{$align->[$idx]} );
#	print STDERR "[DEBUG] set_align_range() idxes=" . join(",", @idxes) . "\n";
	if(@idxes <= 0){
		print STDERR "Warning: empty alignment returned from set_align_range()\n";
		return [ [], [], [] ];
	}
	my $startIdx = $idxes[0];
	my $endIdx = $idxes[@idxes-1];
#	print STDERR "[DEBUG] $startIdx..$endIdx=" . join(",", $startIdx..$endIdx) . "\n";
	return [ [@{$align->[0]}[$startIdx..$endIdx]], [@{$align->[1]}[$startIdx..$endIdx]], [@{$align->[2]}[$startIdx..$endIdx]] ];
}
sub shrink_align_ends_with{
	my ($align, @ends_with) = @_;
	my $alignLen = scalar(@{$align->[2]});
	my %ends_hash = ();
	foreach (@ends_with){
		$ends_hash{$_} = 1;
	}
	my $startIdx = 0;
	my $endIdx = $alignLen - 1;
	while($startIdx <= $endIdx && exists($ends_hash{$align->[2]->[$startIdx]})){
		$startIdx++;
	}
	while($startIdx <= $endIdx && exists($ends_hash{$align->[2]->[$endIdx]})){
		$endIdx--;
	}
	if($startIdx > $endIdx){
		print STDERR "Warning: empty alignment returned from shrink_align_ends_with()\n";
		return [ [], [], [] ];
	}
	return [ [@{$align->[0]}[$startIdx..$endIdx]], [@{$align->[1]}[$startIdx..$endIdx]], [@{$align->[2]}[$startIdx..$endIdx]] ];
}
sub get_align_range{
	my ($align, $idx) = @_;
	my $alignLen = scalar(@{$align->[$idx]});
#	print STDERR "[DEBUG] get_align_range() alignLen=$alignLen\n";
	if($alignLen <= 0){
		return (undef, undef);
	}
	my @ans = ($align->[$idx]->[0], $align->[$idx]->[$alignLen-1]);
	if($ans[0] > $ans[1]){
		@ans = reverse @ans;
	}
	return @ans;
}
sub get_align_strand{
	my ($align) = @_;
	my $alignLen = scalar(@{$align->[1]});
#	print STDERR "[DEBUG] get_align_strand() alignLen=$alignLen\n";
#	print STDERR "[DEBUG] get_align_strand() R range [" . join(", ", get_align_range($align, 1)) . "]\n";
	if($alignLen <= 0){
		return undef;
	}
#	print STDERR "[DEBUG] align->[1] = " . $align->[1] . "\n";
	return ($align->[1]->[$alignLen-1] >= $align->[1]->[0]) ? "+" : "-";
}
sub get_align_Qseq{
	my ($align, $Qseq, $strand) = @_;
	my $alignLen = scalar(@{$align->[1]});
#	print STDERR "[DEBUG] get_align_Qseq() alignLen=$alignLen\n";
	if($alignLen <= 0){
		return "";
	}
	my ($Qstart, $Qend) = get_align_range($align, 0);
	my $ans = substr($Qseq, $Qstart-1, $Qend-$Qstart+1);
	$ans =~ tr/a-z/A-Z/;
	return ($strand eq "+") ? $ans : reverse_complement($ans);
}
sub print_align{
	my ($align, $out, $Qseq, $Rseq) = @_;
	if(!defined($out)){
		$out = *STDOUT;
	}
	my $alignLen = scalar(@{$align->[1]});
	print $out join("\t", "#Qpos", "Rpos", "Op");
	if(defined($Qseq)){
		print $out "\tQseq";
	}
	if(defined($Rseq)){
		print $out "\tRseq";
	}
	print $out "\n";
	my ($i, $pos);
	for($i=0; $i<$alignLen; $i++){
		print $out join("\t", $align->[0]->[$i], $align->[1]->[$i], $align->[2]->[$i]);
		if(defined($Qseq)){
			$pos = $align->[0]->[$i];
			if(abs(round($pos)-$pos) < 0.001){
				print $out "\t".substr($Qseq, $pos-1, 1);
			}else{
				print $out "\t-";
			}
		}
		if(defined($Rseq)){
			$pos = $align->[1]->[$i];
			if(abs(round($pos)-$pos) < 0.001){
				print $out "\t".substr($Rseq, $pos-1, 1);
			}else{
				print $out "\t-";
			}
		}
		print $out "\n";
	}
}

sub round{
	my ($number) = @_;
	return int($number + 0.5);
}
sub retrieve_integer_idxes{
	my @idxes = @_;
	my @ans = ();
	my $idx;
	foreach $idx (@idxes){
		if(abs($idx - round($idx)) < 0.001){
			push(@ans, round($idx));
		}
	}
	return @ans;
}

sub get_V_part_seq{
	my ($name, $seq, $align, $NR) = @_;
	my $should_debug = 0;
#	if($NR==431){  $should_debug = 1;  }
#	print STDERR Dumper($align);
	if($should_debug){
		print STDERR "[DEBUG] align Q range [" . join(", ", get_align_range($align, 0)) . "] for Line $NR\n";
		print STDERR "[DEBUG] align R range [" . join(", ", get_align_range($align, 1)) . "] for Line $NR\n";
	}
	my $strand = get_align_strand($align);
	my ($V_chr, $V_Rstart, $V_Rend, $V_strand) = adjust_start_end(@{$VDJ_H->{$name}});
	my $V_seq = get_ref_seq($V_chr, $V_Rstart, $V_Rend, "+");   # $V_seq is genome-oriented
	if($should_debug){
		print STDERR "[DEBUG] (V_chr, V_Rstart, V_Rend, V_strand) = ($V_chr, $V_Rstart, $V_Rend, $V_strand) for Line $NR\n";
		print STDERR "[DEBUG] Line $NR len(V_seq) = " . length($V_seq) . "\n";
		print STDERR "[DEBUG] V_seq = " . $V_seq . "\n";
		print STDERR "[DEBUG] V_seq (rC) = " . reverse_complement($V_seq) . "\n";
	}
#	$V_seq =~ tr/a-z/A-Z/;
	my ($Qstart, $Qend, $Rstart, $Rend, $Qseq);
	if($V_strand eq "+"){
		$align = set_align_range($align, 1, $V_Rstart, undef);
		if($should_debug){
			print STDERR "[DEBUG] + after set_align_range(), align Q range [" . join(", ", get_align_range($align, 0)) . "] for Line $NR\n";
			print STDERR "[DEBUG] + after set_align_range(), align R range [" . join(", ", get_align_range($align, 1)) . "] for Line $NR\n";
			print_align($align, *STDERR, $seq, get_ref_seq($V_chr, 1, 1e9, "+"));
		}
		$align = shrink_align_ends_with($align, qw/N D/);
		if($should_debug){
			print STDERR "[DEBUG] + after shrink_align_ends_with(), align Q range [" . join(", ", get_align_range($align, 0)) . "] for Line $NR\n";
			print STDERR "[DEBUG] + after shrink_align_ends_with(), align R range [" . join(", ", get_align_range($align, 1)) . "] for Line $NR\n";
			print_align($align, *STDERR, $seq, get_ref_seq($V_chr, 1, 1e9, "+"));
		}
		($Qstart, $Qend) = get_align_range($align, 0);
		($Rstart, $Rend) = get_align_range($align, 1);
		$Qseq = get_align_Qseq($align, $seq, $strand);   # $seq is read-oriented, $Qseq is genome-oriented
		substr($V_seq, $Rstart-$V_Rstart) = $Qseq;
		if($should_debug){
			print STDERR "[DEBUG] + len(Qseq) = " . length($Qseq) . "\n";
			print STDERR "[DEBUG] + len(V_seq) = " . length($V_seq) . "\n";
			print STDERR "[DEBUG] + V_seq = " . $V_seq . "\n";
			print STDERR "[DEBUG] + V_seq (rC) = " . reverse_complement($V_seq) . "\n";
		}
		return ($V_seq, $Qstart, $Qend);   # return V-D-J oriented
	}else{
		$align = set_align_range($align, 1, undef, $V_Rend);
		if($should_debug){
			print STDERR "[DEBUG] + after set_align_range(), align Q range [" . join(", ", get_align_range($align, 0)) . "] for Line $NR\n";
			print STDERR "[DEBUG] + after set_align_range(), align R range [" . join(", ", get_align_range($align, 1)) . "] for Line $NR\n";
			print_align($align, *STDERR, $seq, get_ref_seq($V_chr, 1, 1e9, "+"));
		}
		$align = shrink_align_ends_with($align, qw/N D/);
		if($should_debug){
			print STDERR "[DEBUG] + after shrink_align_ends_with(), align Q range [" . join(", ", get_align_range($align, 0)) . "] for Line $NR\n";
			print STDERR "[DEBUG] + after shrink_align_ends_with(), align R range [" . join(", ", get_align_range($align, 1)) . "] for Line $NR\n";
			print_align($align, *STDERR, $seq, get_ref_seq($V_chr, 1, 1e9, "+"));
		}
		($Qstart, $Qend) = get_align_range($align, 0);
		($Rstart, $Rend) = get_align_range($align, 1);
		$Qseq = get_align_Qseq($align, $seq, $strand);   # $seq is read-oriented, $Qseq is genome-oriented
		substr($V_seq, 0, $Rend-$V_Rstart+1) = $Qseq;
		if($should_debug){
			print STDERR "[DEBUG] - len(Qseq) = " . length($Qseq) . "\n";
			print STDERR "[DEBUG] - len(V_seq) = " . length($V_seq) . "\n";
			print STDERR "[DEBUG] - V_seq = " . $V_seq . "\n";
			print STDERR "[DEBUG] - V_seq (rC) = " . reverse_complement($V_seq) . "\n";
		}
		return (reverse_complement($V_seq), $Qstart, $Qend);   # return V-D-J oriented
	}
}
sub get_J_part_seq{
	my ($name, $seq, $align, $NR) = @_;
	my $strand = get_align_strand($align);
	my ($J_chr, $J_Rstart, $J_Rend, $J_strand) = adjust_start_end(@{$VDJ_H->{$name}});
	my $J_seq = get_ref_seq($J_chr, $J_Rstart, $J_Rend, "+");   # $J_seq is genome-oriented
#	$J_seq =~ tr/a-z/A-Z/;
	my ($Qstart, $Qend, $Rstart, $Rend, $Qseq);
	if($J_strand eq "+"){
		$align = shrink_align_ends_with(set_align_range($align, 1, undef, $J_Rend), qw/N D/);
		($Qstart, $Qend) = get_align_range($align, 0);
		($Rstart, $Rend) = get_align_range($align, 1);
		$Qseq = get_align_Qseq($align, $seq, $strand);   # $seq is read-oriented, $Qseq is genome-oriented
		substr($J_seq, 0, $Rend-$J_Rstart+1) = $Qseq;
		return ($J_seq, $Qstart, $Qend);   # return V-D-J oriented
	}else{
		$align = shrink_align_ends_with(set_align_range($align, 1, $J_Rstart, undef), qw/N D/);
		($Qstart, $Qend) = get_align_range($align, 0);
		($Rstart, $Rend) = get_align_range($align, 1);
		$Qseq = get_align_Qseq($align, $seq, $strand);   # $seq is read-oriented, $Qseq is genome-oriented
		substr($J_seq, $Rstart-$J_Rstart) = $Qseq;
		return (reverse_complement($J_seq), $Qstart, $Qend);   # return V-D-J oriented
	}
}

sub judge_in_frame{
	my ($VpartAddMid_len, $Jstart, $Jend, $Jname) = @_;
	my ($chr, $start, $end, $strand) = adjust_start_end(@{$VDJ_H->{$Jname}});
#	if(@{$VDJ_H->{$Jname}} < 5){
#		print STDERR "Warning: cannot judge in-frame for J=$Jname because it cannot match any J in fa file\n";
#		return "-";
#	}
	my $Jname2 = $VDJ_H->{$Jname}->[5];
	if(!defined($Jname2) || !exists($Jaux_hash->{$Jname2})){
		print STDERR "Warning: cannot judge in-frame for J=$Jname because it is not in aux file\n";
		return "-";
	}
	my $first_codon_start_pos = $Jaux_hash->{$Jname2}->[0];
#	print STDERR "[DEBUG] first_codon_start_pos=$first_codon_start_pos for J $Jname2\n";
#	print STDERR "[DEBUG] strand=$strand\n";
#	print STDERR "[DEBUG] VpartAddMid_len=$VpartAddMid_len\n";
#	print STDERR "[DEBUG] Jstart=$Jstart\n";
#	print STDERR "[DEBUG] start=$start\n";
	if($strand eq "+"){
		if(($VpartAddMid_len + $first_codon_start_pos - ($Jstart-$start)) % 3 == 0){
			return "T";
		}else{
			return "F";
		}
	}else{
		if(($VpartAddMid_len + $first_codon_start_pos - ($end-$Jend)) % 3 == 0){
			return "T";
		}else{
			return "F";
		}
	}
}

my %stop_codon_hash = (
	"TAA" => 1,
	"TAG" => 1,
	"TGA" => 1,
);
sub any_stop_codon{
	my ($seq) = @_;
#	print STDERR "[DEBUG] any_stop_codon() input seq=$seq\n";
	my $len = length($seq);
	my ($ii);
	for($ii=0; $ii<$len; $ii+=3){
#		if($seq eq "GACATCCAGATGACTCAGTCTCCAGCCTCCCTATCTGCATCTGTGGGAGAAACTGTCACCATCACATGTCGAGCAAGTGAGAATATTTACAGTTATTTAGCATGGTATCAGCAGAAACAGGGAAAATCTCCTCAGCTCCTGGTCTATAATGCAAAAACCTTAGCAGAAGGTGTGCCATCAAGGTTCAGTGGCAGTGGATCAGGCACACAGTTTTCTCTGAAGATCAACAGCCTGCAGCCTGAAGATTTTGGGAATTATTACTGTCAACATCATTATGGTACTCCTCCGACGTTCGGTTGAGGCACCAAGCTGGAAATCAAAC"){
#			print STDERR "[DEBUG] $ii ".substr($seq, $ii, 3)." ".exists($stop_codon_hash{substr($seq, $ii, 3)})."\n";
#		}
		if(exists($stop_codon_hash{substr($seq, $ii, 3)})){
			return "T";
		}
	}
	return "F";
}


my @requried_fields = qw/B_Qstart B_Qend Qstart Qend Seq Rname Rstart Rend Strand B_Rname B_Rstart B_Rend B_Strand/;

my @fields;
my %f2i;
my (@F);
my @output_fields;
my ($seq, $B_Qstart, $B_Qend, $P_Qstart, $P_Qend, $MH_len, $prey, $bait, $mid, $pre, $post);
my ($P_Rname, $P_Rstart, $P_Rend, $P_Strand, $B_Rname, $B_Rstart, $B_Rend, $B_Strand, $Rbait, $Rprey);
my ($P_Cigar, $P_align, $B_Cigar, $B_align);
my ($P_Junction, $B_Junction, @P_overlapping_features, @B_overlapping_features);
my ($P_Jdistance, $P_name, $P_isSameStrand, $B_Jdistance, $B_name, $B_isSameStrand);
my ($V_part_seq, $mid_O_seq, $J_part_seq, $VpartAddMid_len, $is_in_frame, $has_any_stop_codon, $is_productive);
my ($V_part_Qstart, $V_part_Qend, $J_part_Qstart, $J_part_Qend);
open(IN, $input_filename) or die "Error: cannot open tlx file $input_filename for input\n";
print STDERR "Now parsing tlx file $input_filename ...\n";
$_ = <IN>;
s/[\r\n]+$//;
@F = split/\t/;
@fields = @F;
for($i=0; $i<@F; $i++){
	$f2i{$F[$i]} = $i;
}
foreach (@requried_fields){
	if(!exists($f2i{$_})){
		die "Error: cannot find colname $_ in tlx file $input_filename\n";
	}
}
@output_fields = @fields;
push(@output_fields, qw/MH_len pre bait mid prey post Bfeature Pfeature Vpart midO Jpart InFrame Stop Productive/);
print join("\t", @output_fields)."\n";
$NR = 1;
while(<IN>){
	$NR++;
	s/[\r\n]+$//;
	@F = split/\t/;
#	if($F[0] eq "M01407:487:000000000-BWKYL:1:1101:17524:3140"){  print STDERR "[DEBUG] M01407:487:000000000-BWKYL:1:1101:17524:3140 read NR = $NR\n";  }
	$seq = $F[$f2i{"Seq"}];
	$B_Qstart = $F[$f2i{"B_Qstart"}];
	$B_Qend = $F[$f2i{"B_Qend"}];
	$P_Qstart = $F[$f2i{"Qstart"}];
	$P_Qend = $F[$f2i{"Qend"}];
	if($P_Qstart < $B_Qstart){
		print STDERR "Warning: prey start < bait start for Line $NR\n";
	}
	
	$P_Rname  = $F[$f2i{"Rname"}];
	$P_Rstart = $F[$f2i{"Rstart"}];
	$P_Rend   = $F[$f2i{"Rend"}];
	$P_Strand = $F[$f2i{"Strand"}];
	$P_Strand = ($P_Strand > 0)? "+" : "-";
	$P_Cigar = $F[$f2i{"Cigar"}];
	$B_Rname  = $F[$f2i{"B_Rname"}];
	$B_Rstart = $F[$f2i{"B_Rstart"}];
	$B_Rend   = $F[$f2i{"B_Rend"}];
	$B_Strand = $F[$f2i{"B_Strand"}];
	$B_Strand = ($B_Strand > 0)? "+" : "-";
	$B_Cigar = $F[$f2i{"B_Cigar"}];
	
	### 2023-02-16, deal with CIGAR, especially when there are any indels in alignments
	$P_align = cigar2align($P_Qstart, $P_Qend, $P_Rstart, $P_Rend, $P_Strand, $P_Cigar, $NR);
	$B_align = cigar2align($B_Qstart, $B_Qend, $B_Rstart, $B_Rend, $B_Strand, $B_Cigar, $NR);
	
	### 2023-02-15, deal with microhomology (MH, overlapping of bait and prey on read): shorten bait or prey
	$MH_len = $B_Qend - $P_Qstart + 1;
	if($MH_len > 0){
		if($MH_shorten_prey_or_bait eq "bait"){
			$B_Qend -= $MH_len;
#			$B_align = trim_align_ends($B_align, $B_Qstart, $B_Qend);
			$B_align = set_align_range($B_align, 0, $B_Qstart, $B_Qend);
			($B_Qstart, $B_Qend) = get_align_range($B_align, 0);
			($B_Rstart, $B_Rend) = get_align_range($B_align, 1);
		}elsif($MH_shorten_prey_or_bait eq "prey"){
			$P_Qstart += $MH_len;
#			$P_align = trim_align_ends($P_align, $P_Qstart, $P_Qend);
			$P_align = set_align_range($P_align, 0, $P_Qstart, $P_Qend);
			($P_Qstart, $P_Qend) = get_align_range($P_align, 0);
			($P_Rstart, $P_Rend) = get_align_range($P_align, 1);
		}else{
			print STDERR "Warning: Unrecognized MH_shorten_prey_or_bait=$MH_shorten_prey_or_bait, which should be either 'bait' or 'prey'\n";
		}
	}
	$P_Junction = ($P_Strand eq "+")? $P_Rstart : $P_Rend;
	$B_Junction = ($B_Strand eq "+")? $B_Rend : $B_Rstart;
	
	$bait = substr($seq, $B_Qstart-1, $B_Qend-($B_Qstart-1));
	$prey = substr($seq, $P_Qstart-1, $P_Qend-($P_Qstart-1));
	$pre = "-";
	if($B_Qstart-1 > 0){
		$pre = substr($seq, 0, $B_Qstart-1);
	}
	$mid = "-";
	if(($P_Qstart-1)-$B_Qend > 0){
		$mid = substr($seq, $B_Qend, ($P_Qstart-1)-$B_Qend);
	}
	$post = "-";
	$post = substr($seq, $P_Qend);
	push(@F, $MH_len, $pre, $bait, $mid, $prey, $post);
#	$bait_to_prey = $bait.$mid.$prey;
	$bait =~ tr/a-z/A-Z/;
	$mid =~ tr/a-z/A-Z/;
	$prey =~ tr/a-z/A-Z/;
	
	### judge overlapping features (low efficiency)
	@P_overlapping_features = get_overlapping_features($P_Rname, $P_Rstart, $P_Rend, $P_Strand, $P_Junction);
	@B_overlapping_features = get_overlapping_features($B_Rname, $B_Rstart, $B_Rend, $B_Strand, $B_Junction);
	
	$P_name = "-";
	if(@P_overlapping_features > 0){
		($P_Jdistance, $P_name, $P_isSameStrand) = @{$P_overlapping_features[0]};
	}
	$B_name = "-";
	if(@B_overlapping_features > 0){
		($B_Jdistance, $B_name, $B_isSameStrand) = @{$B_overlapping_features[0]};
	}
	push(@F, $B_name, $P_name);
	
	### extract V J information, then judge In-frame , Stop codon
	$V_part_seq = "-";
	$J_part_seq = "-";
	$mid_O_seq = "";
	$is_in_frame = "-";
	$has_any_stop_codon = "-";
	$is_productive = "-";
#	print "[DEBUG] $B_name category = ".judge_VDJname_category($B_name)."\n";
#	print "[DEBUG] $P_name category = ".judge_VDJname_category($P_name)."\n";
	if(judge_VDJname_category($B_name) eq "J" && judge_VDJname_category($P_name) eq "V"){
#		print "[DEBUG] bait on J=$B_name, prey on V=$P_name\n";
		if($P_isSameStrand || $B_isSameStrand){
			print STDERR "Warning: bait on J=$B_name, prey on V=$P_name, read same strand as V-D-J (".($P_isSameStrand?"T":"F").($B_isSameStrand?"T":"F").") ?! for Line $NR\n";
		}else{   # not $P_isSameStrand && not $B_isSameStrand
#			$V_part_seq = get_V_part_seq($P_name, $P_isSameStrand ? $prey : reverse_complement($prey), $P_Rstart, $P_Rend, $NR);
			($V_part_seq, $V_part_Qstart, $V_part_Qend) = get_V_part_seq($P_name, $seq, $P_align, $NR);
#			$J_part_seq = get_J_part_seq($B_name, $B_isSameStrand ? $bait : reverse_complement($bait), $B_Rstart, $B_Rend, $NR);
			($J_part_seq, $J_part_Qstart, $J_part_Qend) = get_J_part_seq($B_name, $seq, $B_align, $NR);
			$mid = substr($seq, $J_part_Qend+1-1, $V_part_Qstart-$J_part_Qend-1);
#			if(length($mid) <= 0){  $mid = "-";  }
			$mid_O_seq = reverse_complement($mid);
#			print STDERR "[DEBUG] mid=$mid\n";
			$VpartAddMid_len = length($V_part_seq);
#			if($mid ne "-"){   # deal with empty mid sequence "-"
				$VpartAddMid_len += length($mid);
#			}
			$is_in_frame = judge_in_frame($VpartAddMid_len, $B_Rstart, $B_Rend, $B_name);
			$has_any_stop_codon = any_stop_codon($V_part_seq.$mid_O_seq.$J_part_seq);
			$is_productive = ($is_in_frame eq "T" && $has_any_stop_codon eq "F") ? "T" : "F";
		}
	}elsif(judge_VDJname_category($P_name) eq "J" && judge_VDJname_category($B_name) eq "V"){
#		print "[DEBUG] bait on V=$B_name, prey on J=$P_name\n";
		if(!$P_isSameStrand || !$B_isSameStrand){
			print STDERR "Warning: bait on V=$B_name, prey on J=$P_name, read not same strand as V-D-J (".($P_isSameStrand?"T":"F").($B_isSameStrand?"T":"F").") ?! for Line $NR\n";
		}else{   # $P_isSameStrand && $B_isSameStrand
#			$V_part_seq = get_V_part_seq($B_name, $B_isSameStrand ? $bait : reverse_complement($bait), $B_Rstart, $B_Rend, $NR);
			($V_part_seq, $V_part_Qstart, $V_part_Qend) = get_V_part_seq($B_name, $seq, $B_align, $NR);
#			$J_part_seq = get_J_part_seq($P_name, $P_isSameStrand ? $prey : reverse_complement($prey), $P_Rstart, $P_Rend, $NR);
			($J_part_seq, $J_part_Qstart, $J_part_Qend) = get_J_part_seq($P_name, $seq, $P_align, $NR);
			$mid = substr($seq, $V_part_Qend+1-1, $J_part_Qstart-$V_part_Qend-1);
#			if(length($mid) <= 0){  $mid = "-";  }
			$mid_O_seq = $mid;
			$VpartAddMid_len = length($V_part_seq);
#			if($mid ne "-"){
				$VpartAddMid_len += length($mid);
#			}
			$is_in_frame = judge_in_frame($VpartAddMid_len, $P_Rstart, $P_Rend, $P_name);
			$has_any_stop_codon = any_stop_codon($V_part_seq.$mid_O_seq.$J_part_seq);
			$is_productive = ($is_in_frame eq "T" && $has_any_stop_codon eq "F") ? "T" : "F";
		}
	}else{
		print STDERR "Warning: bait $B_name category = ".judge_VDJname_category($B_name)." and prey $P_name category = ".judge_VDJname_category($P_name)." are not a V-J pair for Line $NR\n";
	}
	
	if(length($mid_O_seq) <= 0){
		$mid_O_seq = "-";
	}
	push(@F, $V_part_seq, $mid_O_seq, $J_part_seq, $is_in_frame, $has_any_stop_codon, $is_productive);
	print join("\t", @F)."\n";
}
close(IN);

0;


