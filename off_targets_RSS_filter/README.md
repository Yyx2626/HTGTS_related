# off_targets_RSS_filter

Some scripts related to scanning cryptic RSS sites (CAC) on the genome, and filtering out strong RSS sites for cryptic RSS off-target analysis

Author: Adam Yongxin Ye @ Boston Children's Hospital / Harvard Medical School


## Pipeline

### Step 1. Scan RSS sites starting with CAC

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


### Step 2. Classify V(D)J-HTGTS junctions as on-targets and off-targets, and filter out good RSS for off-targets

Run `yyx_tlx_filter_RSS.20211102.pl` (see Usage prompts section for details on command-line arguments) to annotate junctions as on-targets and off-targets, and filter out off-targets on good RSS sites

First, it will call tlxbedintersect (`robin_scripts/tlx2BED.pl` and `robin_scripts/pullTLXFromBED.pl`, you need to specify the path to `robin_scripts` folder in the command-line argument `[robin_scripts_dir]`) to extract onTarget junctions located in the input `<onTarget.bed>` and output to `<output_prefix>.onTarget.tlx`; and also extract the remaining noOnTarget junctions and output to `<output_prefix>.noOnTarget.tlx`

Then, it will call `yyx_HTGTS_align.20200831.pl` to align the noOnTarget junctions to nearby CAC sites, while `yyx_HTGTS_align.20200831.pl` will further call `yyx_convert_tlx_to_bw.20200904.pl` and `yyx_bdg_extract_multiply.20200120.pl`; the output file is `<output_prefix>.noOnTarget.slop_15_0.CACTGTG_GTG_aligned.tlx` (if baiting from CE (coding end)) or `<output_prefix>.noOnTarget.slop_0_15.CACAGTG_CAC_aligned.tlx` (if baiting from SE (signal end)).

Finally, it will filter out the aligned noOnTarget junctions if located on the input `<goodRSS.bed>`; the final output file is `<output_prefix>.noOnTarget_CACaligned_noGoodRSS.tlx`.



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


