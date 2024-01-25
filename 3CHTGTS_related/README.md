# 3CHTGTS_related

Some scripts related to 3C-HTGTS normalization and peak calling

Author: Adam Yongxin Ye @ Boston Children's Hospital / Harvard Medical School



## Installation

### Suggested Operation System (OS)

Unlix-like, such as:
- Linux
- MacOS

Prerequisite tools in PATH:
- [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)
- [bedGraphToBigWig](https://hgdownload.soe.ucsc.edu/admin/exe/)

Prerequisite R packages:
- tidyverse
- [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
- [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)

R packages like 'tidyverse' can be installed by `install.packages("tidyverse")` in R CMD

Bioconductor packages like 'GenomicRanges' and 'Biostrings' can be installed according to their webpage, like:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")
BiocManager::install("GenomicRanges")
```



## Pipeline 1. 3C-HTGTS normalization

Prerequisite tools in PATH:
- bedtools
- bedGraphToBigWig

First, run `yyx_normalize_3CHTGTS_tlx.20230501.pl` to remove bait peaks which are variable due to self-ligation level (see Usage prompts section for details in input and output). It will extract and count the junctions located in the specified `[signal_coordinate]` and exclude the junctions located in the specified `[rm_artifact_coordinate]`, and output to `<output_prefix>.junction_count.txt` and `<output_prefix>.signal_rm_artifact.tlx`. Then, it will do scaling normalization on the bigwig signal file, which will scale the signal junction number (excluding artifacts) to the specified `[normalized_to]` junctions.

Note: In order to reduce the impact of the level of self-ligation (circularization), you can visualize the original results in IGV to find out the `[rm_artifact_coordinate]` for the high peaks upstream of the bait-site, and then use this script to filter out these high peaks.

Then, downsample the signal\_rm\_artifact.tlx files by `normalizeTLX_specific.py` (see Usage prompts section for details in input and output).

Note: I added the code to fix the random seed to 1234567 in `normalizeTLX_specific.py` just for reproducibility of the demo dataset.  Users may change the random seed or remove it as they like.  Due to randomness in downsampling, the final number of downstream peak calling may vary a bit.

## Pipeline 2. 3C-HTGTS peak calling

Prerequisite R packages:
- tidyverse
- GenomicRanges
- Biostrings

### Step 1. Collapse (align) 3C-HTGTS junctions to the nearby enzyme cutting sites (CATG by NlaIII)

Run `yyx_3CHTGTS_tlx_to_CATG_dist_filtered_bdg.20240102.r` to collapse 3C-HTGTS junctions to the nearby enzyme cutting sites, and output the junction count on each cutting site.
```
for smpl in AAA BBB; do
 echo smpl=$smpl
 (date
 echo Rscript yyx_3CHTGTS_tlx_to_CATG_dist_filtered_bdg.20240102.r $smpl*.tlx $smpl.CATG_dist10_filtered mm9_AJ851868ins.fa CATG 10 mm9AJ_CATGsites.bed
 time Rscript yyx_3CHTGTS_tlx_to_CATG_dist_filtered_bdg.20240102.r $smpl*.tlx $smpl.CATG_dist10_filtered mm9_AJ851868ins.fa CATG 10 mm9AJ_CATGsites.bed
 date) |& cat >$smpl.CATG_dist10_filtered.log &
done
time wait
date
```
Note: in some bash, '|&' (the shorthand of '2>&1 |') will prompt error '-bash: syntax error near unexpected token &'; it can be simply fixed by replacing '|&' by '2>&1 |'.

The output `$smpl.CATG_dist10_filtered.bdg` will be the input of the next step.


### Step 2. Call 3C-HTGTS peaks using local Poisson tests

Run `yyx_3CHTGTS_callpeak.20240104.r` to call peak regions for each sample.
```
for w in 100; do
for smpl in AAA BBB; do
 echo smpl=$smpl
 (date
 echo Rscript yyx_3CHTGTS_callpeak.20240104.r $smpl.CATG_dist10_filtered.bdg $smpl.CATG_dist10_w${w}_signif $w
 time Rscript yyx_3CHTGTS_callpeak.20240104.r $smpl.CATG_dist10_filtered.bdg $smpl.CATG_dist10_w${w}_signif $w
 date) |& cat >$smpl.CATG_dist10_w${w}_signif.log &
done
time wait
date
done
```

It will calculate the median with a moving window of 100 cutting sites (actually 101 sites: one center, 50 left, and 50 right sites). Then, it will do a Poisson test for each site, with the median as a conservative over-estimation of the lambda parameter of Poisson distribution. Peak summits will be called at the sites with Bonferroni-adjusted p-values < 0.05, and the range of peak region will be determined by progressively extending the two sides to the sites that have local minimum raw p-value and also the raw p-values >= 0.05. Nearby overlapping peak regions were merged as one peak region, and only the “best” (defined by lowest p-value) summit was kept after merging.

The output `$smpl.CATG_dist10_w${w}_signif.merged.bed` will be the input of the next step.


### Step 3. Merge repeats to get robust peak regions

First, run `yyx_multi_bed_overlap.20240104.r` to find out merged peak regions between repeats. Then, run `yyx_process_multi_bdg_to_robust.20240114.pl` to get robust peak regions supported by >50% repeats.
```
mkdir merge_repeats

time echo "Group1 {AAA,BBB}
Group2 {CCC,DDD}" | while read grp smpls; do
 echo $grp $smpls
 in=`eval ls $smpls*.CATG_dist10_w100_signif.merged.bed`
 out=merge_repeats/$grp
 N=`echo $in | awk '{print NF}'`
 echo $N $in
 (date
 echo Rscript yyx_multi_bed_overlap.20240104.r $out.multi $in
 time Rscript yyx_multi_bed_overlap.20240104.r $out.multi $in
 echo perl yyx_process_multi_bdg_to_robust.20240114.pl $out $out.multi.bed $in
 time perl yyx_process_multi_bdg_to_robust.20240114.pl $out $out.multi.bed $in
 date) |& cat >$out.yyx_multi_bed_overlap.log
done
```

The output `$out.multi.bed` by `yyx_multi_bed_overlap.20240104.r` is the input of `yyx_process_multi_bdg_to_robust.20240114.pl`; and the final output are `$out.robust.bed`, `$out.robust.bdg`, `$out.robust_summit.bdg`, and `$out.robust_avg_summit.bdg` (see Usage prompts section for defailts of the output).


### Step 4. Annotate underlying features (optional)

Run `yyx\_anno\_underlying\_features.20240105.r` on annotation bed or bdg/bw files, each time one annotation file, to append annotation columns for each peak region.
```
mkdir check_underlying_features
mkdir check_underlying_features/anno_input

## prepare the annotation bed or bdg/bw files in check_underlying_features/anno_input/


mkdir check_underlying_features/robust_peaks_slop1kb

mkdir check_underlying_features/robust_peaks_slop1kb/query_input
rm -f check_underlying_features/robust_peaks_slop1kb/query_input/*

for grp in Group1 Group2; do
 echo $grp
 bedtools slop -i merge_repeats/$grp.robust.bdg -g mm9_AJ851868ins.chrominfo.txt -b 1000  >check_underlying_features/robust_peaks_slop1kb/query_input/$grp.robust_peaks_slop1kb.bdg
done


rm -f tmp{Input,Output}*

time for query_type in robust_peaks_slop1kb; do
 echo $query_type
 for grp in Group1 Group2; do
  echo $grp
  out=check_underlying_features/$query_type/$grp.anno_CBE_E2A_GROseq.bed
  log=check_underlying_features/$query_type/$grp.anno_CBE_E2A_GROseq.log
  cat check_underlying_features/$query_type/query_input/$grp.$query_type.bed >tmpInput.$grp.$query_type.bed
  ls -lh tmpInput.$grp.$query_type.bed
#  head tmpInput.$grp.$query_type.bed
  (date
  ## CBE motifs
  echo Rscript yyx_anno_underlying_features.20240105.r tmpOutput.$grp.$query_type 1 tmpInput.$grp.$query_type.bed check_underlying_features/anno_input/mm9AJ.fimo_JASPAR_CTCF.atLeast2repeats_parental.absScoreGt13.bed
  time Rscript yyx_anno_underlying_features.20240105.r tmpOutput.$grp.$query_type 1 tmpInput.$grp.$query_type.bed check_underlying_features/anno_input/mm9AJ.fimo_JASPAR_CTCF.atLeast2repeats_parental.absScoreGt13.bed
  cat tmpOutput.$grp.$query_type.bed >tmpInput.$grp.$query_type.bed
  ## E2A ChIPseq
  echo Rscript yyx_anno_underlying_features.20240105.r tmpOutput.$grp.$query_type 1 tmpInput.$grp.$query_type.bed check_underlying_features/anno_input/RAG1KO_E2A_ChIP_GSM546523_treat_pileup.bw
  time Rscript yyx_anno_underlying_features.20240105.r tmpOutput.$grp.$query_type 1 tmpInput.$grp.$query_type.bed check_underlying_features/anno_input/RAG1KO_E2A_ChIP_GSM546523_treat_pileup.bw
  cat tmpOutput.$grp.$query_type.bed >tmpInput.$grp.$query_type.bed
  ## GROseq
  for GROseqID in GROseq{1,2,3}; do
   echo Rscript yyx_anno_underlying_features.20240105.r tmpOutput.$grp.$query_type 1 tmpInput.$grp.$query_type.bed check_underlying_features/anno_input/$GROseqID.sum.bw
   time Rscript yyx_anno_underlying_features.20240105.r tmpOutput.$grp.$query_type 1 tmpInput.$grp.$query_type.bed check_underlying_features/anno_input/$GROseqID.sum.bw
   cat tmpOutput.$grp.$query_type.bed >tmpInput.$grp.$query_type.bed
  done
  cat tmpInput.$grp.$query_type.bed >$out
  rm -f tmpInput.$grp.$query_type.bed tmpOutput.$grp.$query_type.bed
  date) |& cat >$log &
 done
done
time wait
date

time for query_type in robust_peaks_slop1kb; do
 echo $query_type
 for grp in Group1 Group2; do
  echo $grp
  paste <(cat merge_repeats/$grp.robust.bed | perl -ne 's/[\r\n]+$//; $NR++; print $_."\tpeak$NR\n"; ') <(cat check_underlying_features/$query_type/$grp.anno_CBE_E2A_GROseq.bed | perl -pe 's/[.]\t0\t[*]\t//g;') >check_underlying_features/$query_type/$grp.peak_anno_CBE_E2A_GROseq.bed
 done
done
```

You can load the final annotation output `$grp.peak_anno_CBE_E2A_GROseq.bed` into Excel, and set thresholds to determine whether each peak has underlying CBE, E2A or transcription.



## Demo

### 3C-HTGTS normalization

(date
echo python3 normalizeTLX_specific.py 200000 demo/BMpreB_Cer_rep1_result.tlx demo/BMpreB_Cer_rep2_result.tlx
time python3 normalizeTLX_specific.py 200000 demo/BMpreB_Cer_rep1_result.tlx demo/BMpreB_Cer_rep2_result.tlx
date) 2>&1 | tee demo/normalizeTLX_specific.log


Run `yyx_normalize_3CHTGTS_tlx.20230501.pl` to remove bait peaks:
```
time for smpl in BMpreB_Cer_rep{1,2}; do
echo smpl=$smpl
(date
echo perl yyx_normalize_3CHTGTS_tlx.20230501.pl demo/${smpl}_result.200000.tlx ../demo_reference/mm9.fa.fai demo/${smpl} chr6:64515000-73877000 chr6:70659550-70659700
time perl yyx_normalize_3CHTGTS_tlx.20230501.pl demo/${smpl}_result.200000.tlx ../demo_reference/mm9.fa.fai demo/${smpl} chr6:64515000-73877000 chr6:70659550-70659700
date) 2>&1 | cat >demo/${smpl}.yyx_normalize_3CHTGTS_tlx.log &
done
time wait
date
rm -f demo/BMpreB_Cer_rep?.*.{both,pos,neg}.{bw,bdg}
```
It may take about two minutes to run for one sample.

Check the number of lines of the demo results by `wc -l demo/*.tlx`:
```
   44481 demo/BMpreB_Cer_rep1.signal_rm_artifact.tlx
  200001 demo/BMpreB_Cer_rep1_result.200000.tlx
   56253 demo/BMpreB_Cer_rep2.signal_rm_artifact.tlx
  200001 demo/BMpreB_Cer_rep2_result.200000.tlx
```

Downsample the signal\_rm\_artifact.tlx files by `normalizeTLX_specific.py`:
```
(date
echo python3 normalizeTLX_specific.py 40000 demo/BMpreB_Cer_rep1.signal_rm_artifact.tlx demo/BMpreB_Cer_rep2.signal_rm_artifact.tlx
time python3 normalizeTLX_specific.py 40000 demo/BMpreB_Cer_rep1.signal_rm_artifact.tlx demo/BMpreB_Cer_rep2.signal_rm_artifact.tlx
date) 2>&1 | tee demo/normalizeTLX_specific.log
```
It may take about one second to downsample.

Check the number of lines of the demo results by `wc -l demo/*signal_rm_artifact*.tlx`:
```
   40001 demo/BMpreB_Cer_rep1.signal_rm_artifact.40000.tlx
   44481 demo/BMpreB_Cer_rep1.signal_rm_artifact.tlx
   40001 demo/BMpreB_Cer_rep2.signal_rm_artifact.40000.tlx
   56253 demo/BMpreB_Cer_rep2.signal_rm_artifact.tlx
```


### 3C-HTGTS peak calling

First, run `yyx_3CHTGTS_tlx_to_CATG_dist_filtered_bdg.20240102.r` to align 3C-HTGTS junctions to enzyme cutting sites:
```
time for smpl in BMpreB_Cer_rep{1,2}; do
 echo smpl=$smpl
 (date
 echo Rscript yyx_3CHTGTS_tlx_to_CATG_dist_filtered_bdg.20240102.r demo/$smpl.signal_rm_artifact.40000.tlx demo/$smpl.CATG_dist10_filtered ../demo_reference/mm9.fa CATG 10 ../demo_reference/mm9_chr6_CATGsites.bed
 time Rscript yyx_3CHTGTS_tlx_to_CATG_dist_filtered_bdg.20240102.r demo/$smpl.signal_rm_artifact.40000.tlx demo/$smpl.CATG_dist10_filtered ../demo_reference/mm9.fa CATG 10 ../demo_reference/mm9_chr6_CATGsites.bed
 date) 2>&1 | cat >demo/$smpl.CATG_dist10_filtered.log &
done
time wait
date
```
It may take about half one minute to scan CATG (enzyme cutting sites) on the genome (if `mm9_chr6_CATGsites.bed` does not exist), and then about two minutes for each sample.

Note: Due to file size limit (about 200MB) of github, I removed the genome file `../demo_reference/mm9.fa` and limited `mm9_CATGsites.bed` to chr6 as `../demo_reference/mm9_chr6_CATGsites.bed`. Users may download `mm9.fa.gz` from [UCSC](https://hgdownload.cse.ucsc.edu/goldenpath/mm9/bigZips/), then parse it to only keep the canonical chromosomes by `zless mm9.fa.gz | perl -ne 'BEGIN{ $so=0; } if(/^>/){ $so=0; if(/^>chr[^_]+$/){ $so=1; }} if($so){ print; }' >mm9.fa`, and generate fai index file by `samtools faidx mm9.fa`.

Check the number of lines of the demo results by `wc -l demo/*.bdg`:
```
   10306 demo/BMpreB_Cer_rep1.CATG_dist10_filtered.bdg
   10011 demo/BMpreB_Cer_rep2.CATG_dist10_filtered.bdg
```

Second, run `yyx_3CHTGTS_callpeak.20240104.r` to call peaks for each sample:
```
for w in 100; do
for smpl in BMpreB_Cer_rep{1,2}; do
 echo smpl=$smpl
 (date
 echo Rscript yyx_3CHTGTS_callpeak.20240104.r demo/$smpl.CATG_dist10_filtered.bdg demo/$smpl.CATG_dist10_w${w}_signif $w
 time Rscript yyx_3CHTGTS_callpeak.20240104.r demo/$smpl.CATG_dist10_filtered.bdg demo/$smpl.CATG_dist10_w${w}_signif $w
 date) 2>&1 | cat >demo/$smpl.CATG_dist10_w${w}_signif.log &
done
time wait
date
done
```
It may take about half one minute to call peak regions for each sample.

Check the number of lines of the demo results by `wc -l demo/*signif*.bed`:
```
     450 demo/BMpreB_Cer_rep1.CATG_dist10_w100_signif.bed
     131 demo/BMpreB_Cer_rep1.CATG_dist10_w100_signif.merged.bed
     424 demo/BMpreB_Cer_rep2.CATG_dist10_w100_signif.bed
     140 demo/BMpreB_Cer_rep2.CATG_dist10_w100_signif.merged.bed
```

Third, run `yyx_multi_bed_overlap.20240104.r` and `yyx_process_multi_bdg_to_robust.20240114.pl` to compare called peaks and get robust peaks among repeats
```
time echo "BMpreB_Cer BMpreB_Cer_rep{1,2}" | while read grp smpls; do
 echo $grp $smpls
 in=`eval ls demo/$smpls*.CATG_dist10_w100_signif.merged.bed`
 out=demo/$grp
 N=`echo $in | awk '{print NF}'`
 echo $N $in
 (date
 echo Rscript yyx_multi_bed_overlap.20240104.r $out.multi $in
 time Rscript yyx_multi_bed_overlap.20240104.r $out.multi $in
 echo perl yyx_process_multi_bdg_to_robust.20240114.pl $out $out.multi.bed $in
 time perl yyx_process_multi_bdg_to_robust.20240114.pl $out $out.multi.bed $in
 date) 2>&1 | cat >$out.yyx_multi_bed_overlap.log
done
```
It may take about 5 seconds to merge peaks and get robust peaks among repeats.

Check the number of lines of the demo results by `wc -l demo/*robust*`:
```
     100 demo/BMpreB_Cer.robust.bdg
     100 demo/BMpreB_Cer.robust.bed
     100 demo/BMpreB_Cer.robust_avg_summit.bdg
     100 demo/BMpreB_Cer.robust_summit.bdg
```



## Usage prompts

### yyx\_normalize\_3CHTGTS\_tlx.20230501.pl

```
Usage: yyx_normalize_3CHTGTS_tlx.20230501.pl <input.tlx> <chromSize> <output_prefix>
	[signal_coordinate (default: all; example (Igk + nearby): chr6:64515000-73877000; or chr6)]
	[rm_artifact_coordinate (default: none; example: chr6:70675300-70675400 or chr6:1234-5678,chr6:2333-6666)]
	[normalized_to (default: 1000000)] [remove_intermediate_files (default: 1)]

Required tools:
	bedtools
	bedGraphToBigWig

Intermediate output:
	<output_prefix>.  both/pos/neg  .bed , .bdg , .bw (if bedGraphToBigWig exists)
	<output_prefix>.rm_artifact.  both/pos/neg  .bed , .bdg , .bw (if rm_artifact_coordinate != none)

Final output:
	<output_prefix>.junction_count.txt
	<output_prefix>.norm_from_*_to_*.  both/pos/neg  .bdg or .bw (if bedGraphToBigWig exists)
	<output_prefix>.rm_artifact.norm_from_*_to_*.  both/pos/neg  .bdg or .bw (if rm_artifact_coordinate != none)
	<output_prefix>.signal_rm_artifact.tlx

Author: Adam Yongxin Ye @ BCH
Version: 0.1.3 (2023-05-01)
```


### normalizeTLX\_specific.py

```
Usage: python3 normalizeTLX_specific.py <normalize_to> <input1.tlx> [input2.tlx] ...
```
Output:
- <input1>.<normalize_to>.tlx
- [input2].<normalize_to>.tlx



### yyx\_3CHTGTS\_tlx\_to\_CATG\_dist\_filtered\_bdg.20240102.r

```
Usage: Rscript yyx_3CHTGTS_tlx_to_CATG_dist_filtered_bdg.20240102.r <input.tlx> <output_preifx> <ref.fa>
  [motif (default: CATG)] [max_dist (default: 10)] [motif_sites.bed]
```
Output:
- `<output_prefix>`.bed
    Column 1-3: 	genomic coordinate of the cutting site
    Column 4:   	name = '.'
    Column 5:   	junction count (strand separated)
    Column 6:   	junction strand
- `<output_prefix>`.bdg
    Column 1-3: 	genomic coordinate of the cutting site
    Column 4:   	junction count (strand combined)


### yyx\_3CHTGTS\_callpeak.20240104.r
```
Usage: Rscript yyx_3CHTGTS_callpeak.20240104.r <input.CATG_dist10_filtered.bdg> <output_preifx>
  [local_idx_window_size (default: 100)]
  [adj_pVal_method (default: bonferroni)]
  [summit_adj_pVal_threshold (default: 0.05)]
  [foot_pVal_threshold (default: 0.05)]
  [max_output_log10_pVal (default: 40)]
```
Output:
- `<output_prefix>`.bed
- `<output_prefix>`.read\_count.bdg
- `<output_prefix>`.log10\_pVal.bdg
- `<output_prefix>`.merged.bed
- `<output_prefix>`.merged.read\_count.bdg
- `<output_prefix>`.merged.log10\_pVal.bdg
  Format for bed file:
    Column 1-3: 	genomic coordinate of the called peak region
    Column 4:   	name = '.'
    Column 5:   	junction count in the peak region
    Column 6:   	strand = '*'
    Column 7:   	local Poisson test p-value
    Column 8-9: 	genomic coordinate of the summit cutting site
    Column 10:  	junction count at the summit (summit height)


### yyx\_multi\_bed\_overlap.20240104.r
```
Usage: Rscript yyx_multi_bed_overlap.20240104.r <output_prefix> <rep1.bed|bdg> <rep2.bed|bdg> ...
```
Output:
- `<output_prefix>`.bed
    Column 1-3: 	genomic coordinate of the merged region
    Column 4:   	details of the supporting regions in each repeat for the merged region;
                	the number after colon is the line index in each repeat (1-based)
    Column 5:   	the number of repeats supporting the merged region
    Column 6:   	strand = '*'
- `<output_prefix>`.bdg
    Column 1-3: 	genomic coordinate of the merged region
    Column 4:   	the number of repeats supporting the merged region


### yyx\_process\_multi\_bdg\_to\_robust.20240114.pl
```
Usage: perl yyx_process_multi_bdg_to_robust.20240114.pl <output_prefix> <input.multi.bed>
	<rep1.signif.merged.bed> <rep2.signif.merged.bed> ...

Output:
    <output_prefix>.robust.bed
	Column 1:  	chr of merged peak region
	Column 2:  	start of merged peak region (0-based)
	Column 3:  	end of merged peak region (0-based, end-excluded)
	Column 4:  	details of overlapping peak regions in each repeat
	Column 5:  	number of supporting repeats
	Column 6:  	strand is always *
	Column 7:  	best p-value
	Column 8:  	best summit start (1-based)
	Column 9:  	best summit end (1-based, end-included)
	Column 10: 	best summit height
	Column 11: 	average summit coordinate (1-based)
      Note: robust (among repeats) is defined by being supported by >50% (no equal) repeats
            best (among repeats) is defined by lowest p-value
    <output_prefix>.robust.bdg
    <output_prefix>.robust_summit.bdg
    <output_prefix>.robust_avg_summit.bdg
	score = number of supporting repeats

Version: 0.1.0 (2024-01-14)
Author: Adam Yongxin Ye @ BCH
```


### yyx\_anno\_underlying\_features.20240105.r
```
Usage: Rscript yyx_anno_underlying_features.20240105.r  <output_prefix> <should_abs_score> <query.bed|bdg> <anno.bed|bdg|bw>
```
Output:
- `<output_prefix>`.bed

Note:
- For `<anno.bed>`, it will append four columns: (anno) count, density, highest signal, average signal
- For `<anno.bdg|bw>`, it will append two columns: (anno) highest signal, average signal

