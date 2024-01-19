# TranslocPreprocess_patch

A patch to [original HTGTS pipeline](https://robinmeyers.github.io/transloc_pipeline/) to fix the issue of empty Trim column in `preprocess_stats.txt` due to the change of the output format of `fastq-multx`

Author: Adam Yongxin Ye @ Boston Children's Hospital / Harvard Medical School


## Issue and Patch

For some unknown reason, the output format of recent version of `fastq-multx` is changed, which will lead to empty Trim column in `preprocess_stats.txt`.

In order to fix this issue, I made this patch perl script `yyx_fix_multx.20230126.pl` to fix the output format, and modified `TranslocPreprocess.pl` to add the patch step:
```
read_in_meta_file;

unless (defined $indir) {
	create_barcode_file;
	mkdir "$outdir/multx";
	my $multx_cmd = join(" ", "/usr/local/bin/fastq-multx -m $bc_mismatch -d $bc_min_dist -x -b -B $bcfile $read1",
														$paired_end ? "$read2 -o $outdir/multx/%_R1.fq.gz $outdir/multx/%_R2.fq.gz" :
																					"-o $outdir/multx/%.fq.gz");

	@bc_output = Capture($multx_cmd) or croak "Error: fastq-multx failed";
	print "\n@bc_output\n";

	# Yyx added below on 2023-01-26 to fix the problem of fastq-multx output
	$multx_cmd = "yyx_fix_multx.20230126.pl $outdir/multx 1";
	system("/bin/bash", "-c", $multx_cmd);
} else {
	check_existance_of_files;
}
```

Download `yyx_fix_multx.20230126.pl`, put in any folder in your PATH, and make it executable:
```chmod +x yyx_fix_multx.20230126.pl```

And replace `bin/TranslocPreprocess.pl` after installing [original HTGTS pipeline](https://robinmeyers.github.io/transloc_pipeline/) to this modified one; then you should be able to run through `TranslocPreprocess.pl`, and get numbers in Trim column in `preprocess_stats.txt`.



## Usage prompts

### yyx\_fix\_multx.20230126.pl

Usage: ```yyx_fix_multx.20230126.pl <multx_folder_path> [should_run (default:0)]```

