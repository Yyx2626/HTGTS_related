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

`yyx_tlx_filter_RSS.20211102.pl` will call tlxbedintersect (`robin_scripts/tlx2BED.pl` and `robin_scripts/pullTLXFromBED.pl`) and `yyx_HTGTS_align.20200831.pl`.



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


### Step 1. Classify HTGTS-V(D)J-seq junctions as on-targets and off-targets, and filter out good RSS for off-targets

Run `yyx_tlx_filter_RSS.20211102.pl` (see Usage prompts section for details on command-line arguments) to annotate junctions as on-targets and off-targets, and filter out off-targets on good RSS sites

First, it will call tlxbedintersect (`robin_scripts/tlx2BED.pl` and `robin_scripts/pullTLXFromBED.pl`, you need to specify the path to `robin_scripts` folder in the command-line argument `[robin_scripts_dir]`) to extract onTarget junctions located in the input `<onTarget.bed>` and output to `<output_prefix>.onTarget.tlx`; and also extract the remaining noOnTarget junctions and output to `<output_prefix>.noOnTarget.tlx`

Then, it will call `yyx_HTGTS_align.20200831.pl` to align the noOnTarget junctions to nearby CAC sites. The output file is `<output_prefix>.noOnTarget.slop_15_0.CACTGTG_GTG_aligned.tlx` (if baiting from CE (coding end)) or `<output_prefix>.noOnTarget.slop_0_15.CACAGTG_CAC_aligned.tlx` (if baiting from SE (signal end)).

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

This script is intended to scan all cryptic 12/23 RSS sites starting with CAC on the genome sequence, and give each site a score according to the number of matches to the ideal 12/23 RSS (heptamer CACAGTG + 12/23 spacer + nonamer ACAAAAACC).

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

Therefore, if good RSS is defined by score >= 20, it means matching CAC + additional >= 9 bp of the remaining ideal heptamer AGTG and/or nonamer ACAAAACC in the context of a 12-or-23-bp spacer.


### yyx\_scan\_heptamer.20220727.pl

This script is intended to scan all RSS heptamer starting with CAC on the genome sequence, and give each site a score according to the number of matches to the ideal heptamer CACAGTG.

```
Usage: yyx_scan_heptamer.20220727.pl <sequence.fa>
	[should_skip_all_0 (default:1)]
Output: .tsv (8 columns)
	1) chr
	2) position (0-based)
	3) forward base
	4) forward heptamer
	5) forward heptamer score
	6) reverse base
	7) reverse heptamer
	8) reverse heptamer score
Note:
	This script will scan for
		heptamer = CACAGTG
	The 12/23 RSS score is simply calculated as:
		match +2, transition +0, transversion +0
	Score of first 3bp of heptamer (CAC) are multiplied by 10, and minus 58, so that
		matching CAC only has score =  2
		matching CACA has score = 4 = 2 + 2
		matching CACAGTG has score = 10 = 2 + 4*2
		CAC + random bases has expected score = 4 = 2+(2/4)*4

Version: 0.1.0 (2022-07-27)
Author: Adam Yongxin Ye @ BCH
```


### yyx\_tlx\_filter\_RSS.20211102.pl

This script is the top-level wrapper to extract on-targets and off-targets from input HTGTS tlx file, align off-targets to nearby RSS heptamer (CACAGTG) or CAC, and filter out off-targets on good RSS sites.

It will call tlxbedintersect (`robin_scripts/tlx2BED.pl` and `robin_scripts/pullTLXFromBED.pl`) and `yyx_HTGTS_align.20200831.pl`.

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

Note: theoretically, only coding end (CE) may be resected, while signal end (SE) will not. Therefore, this pipeline only allows one-direction searching for CAC from the Junction site on genome: searching for GTG in upstream 15bp ~ downstream 0bp for CE, or searching for CAC in upstream 0bp ~ downstream 15bp for SE.


### yyx\_HTGTS\_align.20200831.pl

This script is intended to align the HTGTS-V(D)J-seq junctions to nearby RSS heptamer (CACAGTG) or CAC.

If the option `[should_output_bw]` is set to true (e.g. 1), it will call `yyx_convert_tlx_to_bw.20200904.pl` and `yyx_bdg_extract_multiply.20200120.pl` to output bw files.

```
Usage: yyx_HTGTS_align.20200831.pl <input.tlx> <ref.fa> <output_prefix>
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

Version: 0.1.0 (2020-08-31)
Author: Adam Yongxin Ye @ BCH
```

The output files: `*_shift.tlx` just appends 3 columns for annotation of the shift information to nearby heptamer or CAC. Then, `*_shift.tlx` is separated into `*_pass.tlx` (nearby CAC is found) and `*_filtered.tlx` (no nearby CAC). Finally, `*_aligned.tlx` is derived from `*_pass.tlx` with Junction coordinate moved to the start position of CAC.


### yyx\_convert\_tlx\_to\_bw.20200904.pl

This script is intended to convert a tlx file to a bw file: the depth of pileup of junctions (by default).

```
Usage: yyx_convert_tlx_to_bw.20200904.pl <input.tlx> <chromSize> <output_prefix>
	[should_strand (default:1)] [normalize_to (default:0)]
	[yyx_bdg_extract_multiply.pl (required for normalization)]
	[junction|prey (default:junction)]

Author: Adam Yongxin Ye @ BCH
Version: 0.1.1 (2020-09-04)
```

### yyx\_convert\_tlx\_to\_bw.20230919.pl

This script is intended to convert a tlx file to a bw file: the depth of pileup of junctions/preys/bait/baitEnd.

```
Usage: yyx_convert_tlx_to_bw.20230919.pl <input.tlx> <chromSize> <output_prefix>
	[should_strand (default:1)] [normalize_to (default:0)]
	[yyx_bdg_extract_multiply.pl (required for normalization)]
	[junction|prey|bait|baitEnd (default:junction)]

Author: Adam Yongxin Ye @ BCH
Version: 0.1.2 (2023-09-19)
```

### yyx\_bdg\_extract\_multiply.20200120.pl

This script is intended to scaling a bdg (bedgraph) file; for example, `<multiply_factor>` = -1 to make the values negative.

```
Usage: yyx_bdg_extract_multiply.20200120.pl <input.bdg|bw> <output.bdg> <multiply_factor> [chr] [start] [end]

Version: 0.1.0 (2020-01-20)
Author: Adam Yongxin Ye @ BCH
```


