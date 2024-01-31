#!/usr/bin/perl

use strict;
use warnings;
use List::Util "shuffle";

my $version = '
Version: 0.1.0 (2021-11-18)
Author: Adam Yongxin Ye @ BCH
';
my $usage = "Usage: $0 <from_number> <to_number> <input.tlx>
	[output.tlx(default:<input>.from_<from>_to_<to>.tlx)]
	[random_seed(default:1234567)]
".$version;

if(@ARGV < 3){
	die $usage;
}
my ($from, $to, $input_filename) = @ARGV;
my $output_filename = $input_filename;
if(@ARGV > 3){
	$output_filename = $ARGV[3];
}else{
	$output_filename =~ s/[.]tlx$//;
	$output_filename .= ".from_" . $from . "_to_" . $to .".tlx";
}
my $rand_seed = 1234567;
if(@ARGV > 4){
	$rand_seed = $ARGV[4];
}

print STDERR "Set random seed = $rand_seed\n";
srand($rand_seed);
if($to > $from){
	die "Error: Cannot downsampling from $from to $to ($to > $from)\n";
}
print STDERR "Random downsampling from $from to $to ...\n";
my @R = shuffle(1..$from);
my %H = ();
my $i;
for($i=1; $i<=$to; $i++){
  $H{$R[$i]} = 1;
}
print STDERR "Read $input_filename, and output to $output_filename ...\n";
open(OUT, ">".$output_filename) or die "Error: cannot open $output_filename for output\n";
open(IN, $input_filename) or die "Error: cannot open $input_filename for input\n";
my $hl=<IN>;
print OUT $hl;
my $NR=0;
while(<IN>){
  $NR++;
  if(exists($H{$NR})){
    print OUT $_;
  }
}

if($NR > $from){
	print STDERR "Warning: input tlx has $NR data lines, more than from=$from; therefore, from number may be misspecified\n";
}

0;

