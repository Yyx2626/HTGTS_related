#!/usr/bin/perl

use strict;
use warnings;


my $version = '
Version: 0.1.2 (2021-11-02)
Author: Adam Yongxin Ye @ BCH
';
my $usage = "Usage: $0 <input.tlx> <output_prefix> <onTarget.bed> <goodRSS.bed> <ref.fa>
	[CE|SE (default:CE)]
	[robin_scripts_dir (default:robin_scripts)]
Options:
	option CE|SE for tlx baiting from coding-end (CE) or signal-end (SE) of RSS
	for CE, I will call yyx_HTGTS_align with parameters:
		slop 15(upstream) 0(downstream)
		motif CACTGTG:GTG , aligned_shift_to_motif_pos 6:2
	for SE, I will call yyx_HTGTS_align with parameters:
		slop 0(upstream) 15(downstream)
		motif CACAGTG:CAC , aligned_shift_to_motif_pos 0:0
Output: (major output files)
	<output_prefix>.onTarget.tlx
	<output_prefix>.noOnTarget.tlx
	<output_prefix>.noOnTarget.slop_15_0.CACTGTG_GTG_aligned.tlx
	(or <output_prefix>.noOnTarget.slop_0_15.CACAGTG_CAC_aligned.tlx)
	<output_prefix>.noOnTarget_CACaligned_noGoodRSS.tlx
".$version;


if(@ARGV < 5){
	die $usage;
}

my ($input_filename, $output_prefix, $onTarget_bed, $goodRSS_bed, $ref_fasta) = @ARGV;
print STDERR "[CMD-ARG] input_filename =\t$input_filename\n";
print STDERR "[CMD-ARG] output_prefix =\t$output_prefix\n";
print STDERR "[CMD-ARG] onTarget_bed =\t$onTarget_bed\n";
print STDERR "[CMD-ARG] goodRSS_bed =\t$goodRSS_bed\n";
print STDERR "[CMD-ARG] ref_fasta =\t$ref_fasta\n";

my $baitCeSe = "CE";
if(@ARGV > 5){
	$baitCeSe = $ARGV[5];
}
print STDERR "[CMD-ARG] baitCeSe =\t$baitCeSe\n";

my $robin_scripts_dir = "robin_scripts";
if(@ARGV > 6){
	$robin_scripts_dir = $ARGV[6];
}
print STDERR "[CMD-ARG] robin_scripts_dir =\t$robin_scripts_dir\n";



my $start_time = time();
print STDERR "[PERL-START] ".scalar(localtime())."\n";


my ($program, $command, $options);
my @filenames = ();


### Step 0. only keep on-target
print STDERR "[Step 0] only keep on-target\n";
$filenames[1] = $output_prefix . ".onTarget.tlx";
if(!exist_file_or_dir($filenames[1], "Skip this tlxbedintersect pipeline...")){
	tlxbedintersect($input_filename, $onTarget_bed, $output_prefix, "onTarget", "-u", $robin_scripts_dir);
}

### Step 1. filter out known on-target
print STDERR "[Step 1] filter out known on-target\n";
$filenames[2] = $output_prefix . ".noOnTarget.tlx";
if(!exist_file_or_dir($filenames[2], "Skip this tlxbedintersect pipeline...")){
	tlxbedintersect($input_filename, $onTarget_bed, $output_prefix, "noOnTarget", "-v", $robin_scripts_dir);
}


### Step 2. align to CAC (actually GTG last position, for Jk-CE-bait (not SE-bait))
print STDERR "[Step 2] align to CAC (actually GTG)\n";
$filenames[12] = $output_prefix . ".noOnTarget";
$options = "15 0 CACTGTG:GTG 6:2";
$filenames[13] = $output_prefix . ".noOnTarget.slop_15_0.CACTGTG_GTG_aligned.tlx";
if($baitCeSe eq "SE"){
	$options = "0 15 CACAGTG:CAC 0:0";
	$filenames[13] = $output_prefix . ".noOnTarget.slop_0_15.CACAGTG_CAC_aligned.tlx";
}
$program = "yyx_HTGTS_align.20200831.pl";
if(! -e $program){  die "Error: cannot find program file $program\n";  }
$command = "perl $program ".$filenames[2]." $ref_fasta ".$filenames[12]." $options 0";
check_file_then_exec_command($filenames[13], $command, 1, 1, 0);


### Step 3. filter out overlapping good RSS sites
print STDERR "[Step 3] filter out good RSS sites\n";
$filenames[22] = $output_prefix . ".noOnTarget_CACaligned_noGoodRSS.tlx";
if(!exist_file_or_dir($filenames[22], "Skip this tlxbedintersect pipeline...")){
	tlxbedintersect($filenames[13], $goodRSS_bed, $output_prefix, "noOnTarget_CACaligned_noGoodRSS", "-v", $robin_scripts_dir);
}


check_elapsed_time($start_time);
print STDERR "[PERL-END] ".scalar(localtime())."\n";

0;




sub tlxbedintersect{
	my ($input_filename, $bed_filename, $output_prefix, $output_suffix, $option, $robin_scripts_dir) = @_;
	### ref: /var/www/html/cgi-bin/tlxbedintersect.cgi
	# 1) /home/zhoudu/software/Scripts/tlx2BED.pl $tlxfile $tlxbed
	# 2) /usr/bin/bedtools intersect -u/-v -a $tlxbed -b $bedfile  >$output
	# 2) /usr/bin/bedtools intersect -c -a $bedfile -b $tlxbed  >$output
	# 3) /home/zhoudu/software/Scripts/pullTLXFromBED.pl $tlxfile $output $tlxoutput

	my ($program, $command);
	my @filenames = ();

	$filenames[2] = $output_prefix . ".tlx2BED.bed";
	$program = "$robin_scripts_dir/tlx2BED.pl";
	if(! -e $program){  die "Error: cannot find program file $program\n";  }
	$command = "perl $program $input_filename ".$filenames[2];
	check_file_then_exec_command($filenames[2], $command, 1, 1, 0);

	$filenames[4] = $output_prefix . ".$output_suffix.bed";
	if($option eq "-c"){
		$command = "bedtools intersect $option -a $bed_filename -b ".$filenames[2]." >".$filenames[4];
		check_file_then_exec_command($filenames[4], $command, 1, 1, 0);
		check_final_file_then_remove_intermediate_file($filenames[4], $filenames[2]);
	}else{
		$command = "bedtools intersect $option -a ".$filenames[2]." -b $bed_filename >".$filenames[4];
		check_file_then_exec_command($filenames[4], $command, 1, 1, 0);

		$filenames[6] = $output_prefix . ".$output_suffix.tlx";
		$program = "$robin_scripts_dir/pullTLXFromBED.pl";
		if(! -e $program){  die "Error: cannot find program file $program\n";  }
		$command = "perl $program $input_filename ".$filenames[4]." ".$filenames[6];
		check_file_then_exec_command($filenames[6], $command, 1, 1, 0);
		check_final_file_then_remove_intermediate_file($filenames[6], [@filenames[2..4]]);
	}
}



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

