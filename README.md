# HTGTS_related

Some scripts related to data analysis of HTGTS (High-Throughput Genome-Wide Translocation Sequencing)

Author: Adam Yongxin Ye @ Boston Children's Hospital / Harvard Medical School

License: MIT license

Citation: Zhang Y, Li X, Ba Z, Lou J, Gaertner KE, Zhu T, Lin X, Ye AY, Alt FW, Hu H. Molecular basis for differential Igk versus Igh V(D)J joining mechanisms. Nature. 2024 Jun;630(8015):189-197. doi: 10.1038/s41586-024-07477-y. Epub 2024 May 29. PMID: 38811728; PMCID: PMC11153149.


## [TranslocPreprocess_patch](https://github.com/Yyx2626/HTGTS_related/tree/master/TranslocPreprocess_patch)

A patch to [original HTGTS pipeline](https://robinmeyers.github.io/transloc_pipeline/) to fix the issue of empty Trim column in `preprocess_stats.txt` due to the change of the output format of `fastq-multx`


## [off_targets_RSS_filter](https://github.com/Yyx2626/HTGTS_related/tree/master/off_targets_RSS_filter)

Some scripts related to scanning cryptic RSS sites (CAC) on the genome, and filtering out strong RSS sites for cryptic RSS off-target analysis


## [3CHTGTS_related](https://github.com/Yyx2626/HTGTS_related/tree/master/3CHTGTS_related)

Some scripts related to 3C-HTGTS normalization and peak calling

