# off_targets_RSS_filter

Some scripts related to scanning cryptic RSS sites (CAC) on the genome, and filtering out strong RSS sites for cryptic RSS off-target analysis

Author: Adam Yongxin Ye @ Boston Children's Hospital / Harvard Medical School



## Installation

### Suggested Operation System (OS)

Unlix-like, such as:
- Linux
- MacOS

Prerequisite perl packages:
- Text::CSV

If users have not installed the required perl packages before, they may execute commands like `perl -MCPAN -e'install Text::CSV'` to install them (which may need sudo privilege or contact server administrator).

Users can download and run these perl scripts by `perl xxx.pl` directly in this folder. Or users may move or copy the scripts into any of the folders in their PATH (environmental variable, which can be seen be `echo $PATH`), then the scripts can be directly run by `xxx.pl`.

`yyx_tlx_filter_RSS.20211102.pl` will call tlxbedintersect (`robin_scripts/tlx2BED.pl` and `robin_scripts/pullTLXFromBED.pl`) and `yyx_HTGTS_align.20200831.pl`, which will further call `yyx_convert_tlx_to_bw.20200904.pl`.



## Pipeline

### Step 0. Scan RSS sites starting with CAC

Run `yyx_scan_RSS.20210102.pl` on your genome sequence, which will report all sites starting with CAC.

Note: Because the inefficiency of this scanning script, it may take long time to scan for the whole genome. Therefore, if you are using a server with multiple CPUs, then you may separate the whole genome into each chromosome, run the scanning script on each chromosome in parallel, and finally concatenate the results.

The output is in STDOUT, in tsv (tab-separated values) format (see Usage prompts section for the detailed description of each column and score definition), so you can redirect it into a output file, like
```time perl yyx_scan_RSS.20210102.pl mm9.chr6.fa >mm9.chr6.scan_RSS.tsv```

Then, use the following commands in bash to extract good RSS (with score >= 20) and convert into bed format
```
time cat mm9.chr6.scan_RSS.tsv | perl -ne 's/[\r\n]+$//; @F=split/\t/; $s=0; $n="";
 if($F[4]>=20){ $s=$F[4]; $n.=":F12RSS".$F[4]; }
 if($F[5]>=20){ if($F[5]>$s){ $s=$F[5]; } $n.=":F23RSS".$F[5]; }
 if($F[8]>=20){ if($F[8]>$s){ $s=$F[8]; } $n.=":R12RSS".$F[8]; }
 if($F[9]>=20){ if($F[9]>$s){ $s=$F[9]; } $n.=":R23RSS".$F[9]; }
 if($s > 0){ print join("\t", $F[0], $F[1], $F[1]+1, $F[2].$n, $s)."\n"; }' >mm9.chr6.goodRSS_scoreGe20.bed
```

Note: RSS score >= 20 means: compared to the ideal RSS site (heptamer(CACAGTG) + 12/23-bp spacer + nonamer(ACAAAACC)), it requires CAC and addtional >= 9bp matches to the remaining ideal heptamer (AGTG) and nonamer in the context of a 12-or-23-bp spacer, i.e. at most 4bp mismatches to the ideal RSS site.


### Step 1. Classify V(D)J-HTGTS junctions as on-targets and off-targets, and filter out good RSS for off-targets

Run `yyx_tlx_filter_RSS.20211102.pl` (see Usage prompts section for details on command-line arguments) to annotate junctions as on-targets and off-targets, and filter out off-targets on good RSS sites

First, it will call tlxbedintersect (`robin_scripts/tlx2BED.pl` and `robin_scripts/pullTLXFromBED.pl`, you need to specify the path to `robin_scripts` folder in the command-line argument `[robin_scripts_dir]`) to extract onTarget junctions located in the input `<onTarget.bed>` and output to `<output_prefix>.onTarget.tlx`; and also extract the remaining noOnTarget junctions and output to `<output_prefix>.noOnTarget.tlx`

Then, it will call `yyx_HTGTS_align.20200831.pl` to align the noOnTarget junctions to nearby CAC sites, while `yyx_HTGTS_align.20200831.pl` will further call `yyx_convert_tlx_to_bw.20200904.pl` and `yyx_bdg_extract_multiply.20200120.pl`; the output file is `<output_prefix>.noOnTarget.slop_15_0.CACTGTG_GTG_aligned.tlx` (if baiting from CE (coding end)) or `<output_prefix>.noOnTarget.slop_0_15.CACAGTG_CAC_aligned.tlx` (if baiting from SE (signal end)).

Finally, it will filter out the aligned noOnTarget junctions if located on the input `<goodRSS.bed>`; the final output file is `<output_prefix>.noOnTarget_CACaligned_noGoodRSS.tlx`.



## Demo

Run Step 1 for demo:
```
(date
echo perl yyx_tlx_filter_RSS.20211102.pl demo/one_sample.tlx demo/one_sample demo/IgK_mm9_ontargetIntersect_40updownRSS.bed goodRSS_scoreGe20_examples/mm9.chr6.goodRSS_scoreGe20.bed ../demo_reference/mm9.fa CE robin_scripts
time perl yyx_tlx_filter_RSS.20211102.pl demo/one_sample.tlx demo/one_sample demo/IgK_mm9_ontargetIntersect_40updownRSS.bed goodRSS_scoreGe20_examples/mm9.chr6.goodRSS_scoreGe20.bed ../demo_reference/mm9.fa CE robin_scripts
date) 2>&1 | tee demo/one_sample.yyx_tlx_filter_RSS.log
```
It may take ten seconds or less than one minute to run through it.

Check the number of lines of the demo results by `wc -l demo/one_sample.*.tlx`:
```
    1695 demo/one_sample.noOnTarget.slop_15_0.CACTGTG_GTG_aligned.tlx
     361 demo/one_sample.noOnTarget.slop_15_0.CACTGTG_GTG_filtered.tlx
    1695 demo/one_sample.noOnTarget.slop_15_0.CACTGTG_GTG_pass.tlx
    2055 demo/one_sample.noOnTarget.slop_15_0.CACTGTG_GTG_shift.tlx
    2055 demo/one_sample.noOnTarget.tlx
     784 demo/one_sample.noOnTarget_CACaligned_noGoodRSS.tlx
  111994 demo/one_sample.onTarget.tlx
```

Note: Due to file size limit (about 200MB) of github, I removed the genome file `../demo_reference/mm9.fa`. Users may download `mm9.fa.gz` from [UCSC](https://hgdownload.cse.ucsc.edu/goldenpath/mm9/bigZips/), then parse it to only keep the canonical chromosomes by `zless mm9.fa.gz | perl -ne 'BEGIN{ $so=0; } if(/^>/){ $so=0; if(/^>chr[^_]+$/){ $so=1; }} if($so){ print; }' >mm9.fa`, and generate fai index file by `samtools faidx mm9.fa`.



## Usage prompts

### yyx\_scan\_RSS.20210102.pl

```
Usage: yyx_scan_RSS.20210102.pl <sequence.fa>
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

Version: 0.1.1 (2021-01-02)
Author: Adam Yongxin Ye @ BCH
```

### yyx\_tlx\_filter\_RSS.20211102.pl

```
Usage: yyx_tlx_filter_RSS.20211102.pl <input.tlx> <output_prefix> <onTarget.bed> <goodRSS.bed> <ref.fa>
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

Version: 0.1.2 (2021-11-02)
Author: Adam Yongxin Ye @ BCH
```


