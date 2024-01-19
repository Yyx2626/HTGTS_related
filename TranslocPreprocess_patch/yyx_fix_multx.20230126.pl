#!/usr/bin/perl

use strict;
use warnings;

my $version = '
Version: 0.1.0 (2023-01-26)
Author: Adam Yongxin Ye @ BCH
';
my $usage = "Usage: $0 <multx_folder_path> [should_run (default:0)]
".$version;

if(@ARGV < 1){
	die $usage;
}

my $multx_folder_path = shift(@ARGV);
my $should_run = 0;
if(@ARGV > 0){
	$should_run = shift(@ARGV);
}

my $filename;
my $command;

opendir(DIR, $multx_folder_path) or die "Error: cannot opendir $multx_folder_path: $!\n";
while($filename = readdir(DIR)){
	if($filename =~ /.fq.gz$/){
		print STDERR "[Note] Now processing '".$filename."' ...\n";
		
		$command = "mv -f \"$multx_folder_path/$filename\" \"$multx_folder_path/$filename.old\"";
		print STDERR "[CMD] $command\n";
		if($should_run){  system("/bin/bash", "-c", $command);  }
		
		$command = "zcat \"$multx_folder_path/$filename.old\" | perl -pe 's/ +.*\$//g;' | gzip -c >\"$multx_folder_path/$filename\"";
		print STDERR "[CMD] $command\n";
		if($should_run){  system("/bin/bash", "-c", $command);  }
		
		$command = "rm -f \"$multx_folder_path/$filename.old\"";
		print STDERR "[CMD] $command\n";
		if($should_run){  system("/bin/bash", "-c", $command);  }
	}
}
closedir(DIR);

0;
