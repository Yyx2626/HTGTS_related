#!/usr/bin/perl

use strict;
use warnings;

my $version = '
Version: 0.1.0 (2020-08-31)
Author: Adam Yongxin Ye @ BCH
';
my $usage = "Usage: $0 <input.tlx> <ref.fa> <output_prefix>
	[slop_bp_upstream (default:15)] [slop_bp_downstream (default:slop_bp_upstream)]
	[motif (default: CACTGTG:GTG)] [aligned_shift_to_motif_pos (default: 6:2)]
	[should_output_bw (default:1)] [yyx_convert_tlx_to_bw.pl] [yyx_bdg_extract_multiply.pl]

Output:
	<output_prefix>.bed
	<output_prefix>.slop_*_*.<motif>_shift.bed
	<output_prefix>.slop_*_*.<motif>_shift.tlx
	<output_prefix>.slop_*_*.<motif>_pass.tlx
	<output_prefix>.slop_*_*.<motif>_filtered.tlx
	<output_prefix>.slop_*_*.<motif>_aligned.tlx
".$version;

if(@ARGV < 3){
	die $usage;
}
my ($input_filename, $ref_filename, $output_prefix) = @ARGV;
my $slop_bp_left = 15;
if(@ARGV > 3){
	$slop_bp_left = $ARGV[3];
}
my $slop_bp_right = $slop_bp_left;
if(@ARGV > 4){
	$slop_bp_right = $ARGV[4];
}
my $input_motif = "CACTGTG:GTG";
if(@ARGV > 5){
	$input_motif = $ARGV[5];
}
my $aligned_shift_to_motif_pos = "6:2";
if(@ARGV > 6){
	$aligned_shift_to_motif_pos = $ARGV[6];
}
my $should_output_bw = 1;
if(@ARGV > 7){
	$should_output_bw = $ARGV[7];
}
my $yyx_convert_tlx_to_bw = "yyx_convert_tlx_to_bw.20200904.pl";
if(@ARGV > 8){
	$yyx_convert_tlx_to_bw = $ARGV[8];
}
my $yyx_bdg_extract_multiply = "yyx_bdg_extract_multiply.20200120.pl";
if(@ARGV > 9){
	$yyx_bdg_extract_multiply = $ARGV[9];
}

my $motif_txt = $input_motif;
$motif_txt =~ s/:/_/g;


my $start_time = time();
print STDERR "[PERL-START] ".scalar(localtime())."\n";


my $command;
my @filenames = ();
$filenames[1] = $output_prefix . ".bed";
$filenames[2] = $output_prefix . ".slop_".$slop_bp_left."_".$slop_bp_right.".bed";
$filenames[3] = $output_prefix . ".slop_".$slop_bp_left."_".$slop_bp_right.".fa.tsv";
$filenames[4] = $output_prefix . ".slop_".$slop_bp_left."_".$slop_bp_right.".fa.bed";
$filenames[6]  = $output_prefix . ".slop_".$slop_bp_left."_".$slop_bp_right.".".$motif_txt."_shift.bed";
$filenames[11] = $output_prefix . ".slop_".$slop_bp_left."_".$slop_bp_right.".".$motif_txt."_shift.tlx";
$filenames[12] = $output_prefix . ".slop_".$slop_bp_left."_".$slop_bp_right.".".$motif_txt."_pass.tlx";
$filenames[13] = $output_prefix . ".slop_".$slop_bp_left."_".$slop_bp_right.".".$motif_txt."_filtered.tlx";
$filenames[14] = $output_prefix . ".slop_".$slop_bp_left."_".$slop_bp_right.".".$motif_txt."_aligned.tlx";

print STDERR "# Step 1. convert tlx to bed\n";
print STDERR "##  input=".$input_filename."  output=".$filenames[1]."\n";
if(! exist_file_or_dir(@filenames[1,4,6,11..13], "So skip Step 1")){
	open(OUT, ">".$filenames[1]) or die "Error: cannot output to ".$filenames[1]."\n";
	open(IN, $input_filename) or die "Error: cannot input from ".$input_filename."\n";

	my $retrieve_element = "junction";
	my $headline = <IN>;
	$headline =~ s/[\r\n]+$//;
	my @F = split(/\t/, $headline);
	my %f2i = ();
	my $i;
	for($i=0; $i<@F; $i++){
		$f2i{$F[$i]} = $i;
	}
	my @required_fields = qw/Qname Rname Strand/;
	if($retrieve_element eq "junction"){
		push(@required_fields, "Junction");
	}elsif($retrieve_element eq "prey"){
		push(@required_fields, qw/Rstart Rend/);
	}else{
		die "Error: unknown retrieve $retrieve_element\n";
	}
	foreach (@required_fields){
		if(!exists($f2i{$_})){
			die "Error: cannot find column $_ in input tlx file\n";
		}
	}
	my $rowidx = 1;
	my ($chr, $start,$end, $name, $strand);
	while(<IN>){
		s/[\r\n]+$//;
		if(/^\s*$/){ next; }
		@F = split/\t/;
		$chr = $F[$f2i{"Rname"}];
		$end = 0;
		$start = 0;
		if($retrieve_element eq "junction"){
			$end = $F[$f2i{"Junction"}];
			$start = $end - 1;
		}else{
			$end = $F[$f2i{"Rend"}];
			$start = $F[$f2i{"Rstart"}] - 1;
		}
		$name = $F[$f2i{"Qname"}];
		$strand = "+";
		if($F[$f2i{"Strand"}] < 0){
			$strand = "-";
		}
		print OUT join("\t", $chr, $start, $end, $name, $rowidx, $strand)."\n";
		$rowidx++;
	}

	close(IN);
	close(OUT);
}

print STDERR "# Step 2. slop upstream ".$slop_bp_left." bp and downstream ".$slop_bp_right." bp, and retrieve local sequence\n";
print STDERR "##  input=".$filenames[1]."  output=".$filenames[4]."\n";
if(! exist_file_or_dir(@filenames[4,6,11..13], "So skip Step 2")){
	if(! -f $ref_filename.".fai" ){
		die "Error: cannot find fai file ".$ref_filename.".fai\n";
	}
	$command = "bedtools slop -s -l ".$slop_bp_left." -r ".$slop_bp_right." -i ".$filenames[1]." -g ".$ref_filename.".fai >".$filenames[2];
	check_file_then_exec_command([@filenames[2..4]], $command, 1, 1, 0);

	$command = "bedtools getfasta -s -tab -fi ".$ref_filename." -bed ".$filenames[2]." >".$filenames[3];
	check_file_then_exec_command([@filenames[3..4]], $command, 1, 1, 0);

	$command = "paste ".$filenames[2]." ".$filenames[3]." >".$filenames[4];
	check_file_then_exec_command($filenames[4], $command, 1, 1, 0);

	check_final_file_then_remove_intermediate_file($filenames[4], [@filenames[2..3]]);
}

sub find_all{
	my ($pattern, $subject) = @_;
	my @ans = ();
	my $start;
	while($subject =~ /($pattern)/g){
		$start = pos($subject) - length($1);
		push(@ans, $start);
		pos($subject) = $start + 1;
	}
	return @ans;
}
sub min{
	my $min=undef;
	my $idx=undef;
	my $i;
	for($i=0; $i<@_; $i++){
		if(!defined($min) || $_[$i] < $min){
			$min = $_[$i] ;
			$idx = $i;
		}
	}
	return ($min, $idx);
}

print STDERR "# Step 3. find nearest motif ".$input_motif.", and shift\n";
print STDERR "##  input=".$filenames[4]."  output=".$filenames[6]."\n";
my @motifs = split(/:/, $input_motif);
my @motif_shift = split(/:/, $aligned_shift_to_motif_pos);
if(@motif_shift != @motifs){
	die "Error: aligned_shift_to_motif_pos should have the same length of ':'-separated values as motif\n";
}
my $i;
if(! exist_file_or_dir(@filenames[6,11..13], "So skip Step 3")){
	open(OUT, ">".$filenames[6]) or die "Error: cannot output to ".$filenames[6]."\n";
	open(IN, $filenames[4]) or die "Error: cannot input from ".$filenames[4]."\n";

	my (@F, @P, $abs_shift, $idx, $shift, $pos);
	while(<IN>){
		s/[\r\n]+$//; @F=split/\t/;
		for($i=0; $i<@motifs; $i++){
			@P = map{ $_ - $slop_bp_left + $motif_shift[$i]; } find_all($motifs[$i], uc($F[7]));
			if(@P > 0){
				last;
			}
		}
		if(@P <= 0){
			$shift = "NA";
		}else{
#			print STDERR join(", ", @P)."\n";
			($abs_shift, $idx) = min(map{ abs($_) } @P);
			$shift = $P[$idx];
#			print STDERR join(", ", $abs_shift, $idx, $shift)."\n";
			if($F[5] eq "+"){
				$pos = $F[1] + $slop_bp_left + 1;
				$pos += $shift;
			}else{
				$pos = $F[2] - $slop_bp_left;
				$pos -= $shift;
			}
		}
		print OUT join("\t", $F[0], $pos-1, $pos, @F[3..7], $shift)."\n";
	}

	close(IN);
	close(OUT);
	
	check_final_file_then_remove_intermediate_file($filenames[6], $filenames[4]);
}


print STDERR "# Step 4. merge original tlx, and split\n";
print STDERR "##  input=".$input_filename."  input=".$filenames[6]."\n";
for($i=1; $i<=4; $i++){
	print STDERR "##  output=".$filenames[10+$i]."\n";
}
if(! exist_file_or_dir(@filenames[11..13], "So skip Step 4")){
	open(OUT, ">".$filenames[11]) or die "Error: cannot output to ".$filenames[11]."\n";
	open(OUT2, ">".$filenames[12]) or die "Error: cannot output to ".$filenames[12]."\n";
	open(OUT3, ">".$filenames[13]) or die "Error: cannot output to ".$filenames[13]."\n";
	open(OUT4, ">".$filenames[14]) or die "Error: cannot output to ".$filenames[14]."\n";

	open(IN, $input_filename) or die "Error: cannot input from ".$input_filename."\n";
	open(IN2, $filenames[6]) or die "Error: cannot input from ".$filenames[6]."\n";

	my $headline = <IN>;
	$headline =~ s/[\r\n]+$//;
	my @F = split(/\t/, $headline);
	my %f2i = ();
	for($i=0; $i<@F; $i++){
		$f2i{$F[$i]} = $i;
	}
	push(@F, qw/getfasta_title getfasta_seq/, $motif_txt."_shift");
	print OUT join("\t", @F)."\n";
	print OUT2 join("\t", @F)."\n";
	print OUT3 join("\t", @F)."\n";
	print OUT4 join("\t", @F)."\n";
	my $rowidx = 1;
	my ($line, @G);
	while(<IN>){
		s/[\r\n]+$//;
		if(/^\s*$/){ next; }
		@F = split/\t/;
		$line = <IN2>;
		$line =~ s/[\r\n]+$//;
		@G = split(/\t/, $line);
		if($F[$f2i{"Qname"}] ne $G[3]){
			print STDERR "Warning: $rowidx-th reads are not ID matched !?\t".$F[$f2i{"Qname"}]."\t".$G[3]."\n"
		}
		push(@F, @G[6..8]);
		print OUT join("\t", @F)."\n";
		if($G[8] ne "NA"){
			print OUT2 join("\t", @F)."\n";
			$F[$f2i{"Junction"}] = $G[2];
			print OUT4 join("\t", @F)."\n";
		}else{
			print OUT3 join("\t", @F)."\n";
		}
		$rowidx++;
	}

	close(IN2);
	close(IN);
	close(OUT4);
	close(OUT3);
	close(OUT2);
	close(OUT);
}

if($should_output_bw){
	print STDERR "# Step 5. generate bdg and bw\n";
	if(! -f $ref_filename.".fai" ){
		die "Error: cannot find fai file ".$ref_filename.".fai\n";
	}
	if(! -f $yyx_convert_tlx_to_bw ){
		die "Error: cannot find yyx_convert_tlx_to_bw script file ".$yyx_convert_tlx_to_bw."\n";
	}
	if(! -f $yyx_bdg_extract_multiply ){
		die "Error: cannot find yyx_bdg_extract_multiply script file ".$yyx_bdg_extract_multiply."\n";
	}
	$filenames[21] = $output_prefix . ".slop_".$slop_bp_left."_".$slop_bp_right.".".$motif_txt."_shift";
	$filenames[22] = $output_prefix . ".slop_".$slop_bp_left."_".$slop_bp_right.".".$motif_txt."_pass";
	$filenames[23] = $output_prefix . ".slop_".$slop_bp_left."_".$slop_bp_right.".".$motif_txt."_filtered";
	$filenames[24] = $output_prefix . ".slop_".$slop_bp_left."_".$slop_bp_right.".".$motif_txt."_aligned";
	for($i=1; $i<=4; $i++){
		print STDERR "##  input=".$filenames[10+$i]."  output=".$filenames[20+$i]."\n";
		$command = "perl ".$yyx_convert_tlx_to_bw." ".$filenames[10+$i]." ".$ref_filename.".fai ".$filenames[20+$i]." 1 0 ".$yyx_bdg_extract_multiply;
		check_file_then_exec_command($filenames[20+$i], $command, 1, 1, 0);
	}
}


check_elapsed_time($start_time);
print STDERR "[PERL-END] ".scalar(localtime())."\n";

0;





sub generateRandomString{
	my ($length) = @_;
	my $ans = "";
	my ($random, $random2);
	my ($i, $tmp);
	for($i=0; $i<$length; $i++){
		$random = rand();
		$random2 = $random * 2 - int($random * 2);
		if($random * 2 < 1){
			$tmp = int(ord('A') + $random2 * 26);
			if($tmp > ord('Z')){  $tmp = ord('Z'); }
			$tmp = chr($tmp);
			$ans .= $tmp;
		}else{
			$tmp = int(ord('a') + $random2 * 26);
			if($tmp > ord('z')){  $tmp = ord('z'); }
			$tmp = chr($tmp);
			$ans .= $tmp;
		}
	}
	return($ans);
}



sub check_final_file_then_remove_intermediate_file{
	my ($final_filename, $intermediate_filename) = @_;
	my $command;
	if(ref($intermediate_filename) eq ""){
		$intermediate_filename = [$intermediate_filename];
	}
	if(exist_file_or_dir($final_filename, "So remove intermediate files...")){
		$command = "rm ".join(" ", @{$intermediate_filename});
		print STDERR "[PERL-SYSTEM] ".$command."\n";
		system($command);
	}
}


sub check_file_then_exec_command{
	my ($filename, $command, $should_time, $error_stop, $not_run) = @_;
	my $start_time = time();

	print STDERR "[PERL-SYSTEM] ".$command."\n";
	if(exist_file_or_dir($filename, "Skip this above command...")){
		return;
	}

	if(!(defined($not_run) && $not_run!=0)){
		if(system("/bin/bash", "-c", $command)!=0){
			if(defined($error_stop) && $error_stop!=0){
				die "Error: when exec last system command, return value = $?\n";
			}
		}
	}

	if(defined($should_time) && $should_time){
		check_elapsed_time($start_time);
	}
}

sub check_elapsed_time{
	my ($start_time, $end_time, $elapsed_time);
	my ($hour, $min, $sec, $day);
	$start_time = shift(@_);
	$end_time = time();
	$elapsed_time = $end_time - $start_time;
	$day = int($elapsed_time / (3600*24));
	$hour = int($elapsed_time % (3600*24) / 3600);
	$min = int($elapsed_time % 3600 / 60);
	$sec = $elapsed_time % 60;
	$elapsed_time = "";
	if($day>0){ $elapsed_time .= $day."day "; }
	if($hour>0){ $elapsed_time .= $hour."h"; }
	if($min>0){ $elapsed_time .= $min."min"; }
	if($sec>0 || $elapsed_time eq ""){ $elapsed_time .= $sec."s"; }
	print STDERR "[PERL-SYSTEM-TIME] ".$elapsed_time."\n";
}

sub exist_file_or_dir{
	my ($filenames, $str) = @_;
	my $returnValue = 0;
	if(ref($filenames) eq ""){
		$filenames = [$filenames];
	}
	if(!defined($str)){
		$str = "";
	}

	foreach my $filename (@{$filenames}){
		if(defined($filename) && -e $filename){
			if(-d $filename){
				if(! check_is_empty_dir($filename)){
					print STDERR "[CHECK-EXIST] Dir ".$filename." has already existed, and not empty. ".$str."\n";
					$returnValue = 1;
					last;
				}
			}elsif(-s $filename >= 100){
				print STDERR "[CHECK-EXIST] File ".$filename." has already existed. ".$str."\n";
				$returnValue = 1;
				last;
			}
		}
	}

	return $returnValue;
}


sub check_is_empty_dir{
	my ($dirname) = @_;
	my $dirHandle;
	my @contents;
	if(-d $dirname){
		if(opendir($dirHandle, $dirname)){
			@contents = readdir($dirHandle);
			closedir($dirHandle);
			if(scalar(@contents)>2){		# empty dir has . and ..
				return 0;		# not empty dir
			}else{
				return 1;		# is empty dir
			}
		}else{
			print STDERR ("Warning: cannot open dir $dirname\n");
			return -1;		# Cannot open dir
		}
	}else{
		print STDERR ("Warning: Not a dir is attemped to be checked\n");
		return -2;		# Not a dir
	}
}

