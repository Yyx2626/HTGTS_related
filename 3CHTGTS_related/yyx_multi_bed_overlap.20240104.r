#### Usage: Rscript this.r <output_prefix> <rep1.bed|bdg> <rep2.bed|bdg> ...

#### 2023-01-02 morning, Fred asked to develope a new peak calling program that can identify peaks out of local background noise level, not thresholding on signal height.

args = commandArgs(TRUE)
argIdx = 3
if(length(args) < argIdx){
	stop("not enough command-line parameters
Usage: Rscript this.r <output_prefix> <rep1.bed|bdg> <rep2.bed|bdg> ...")
}

output_prefix = args[1]
input_filenames = args[2:length(args)]



options(stringsAsFactors=FALSE)

library(tidyverse)
library(GenomicRanges)
#library(Biostrings)



`%.%` = function(x,y) paste0(x,y)

yyx_unlist_a_list_of_GRanges = function(GRanges_list){
	if(!is.list(GRanges_list)){
		return(GRanges_list)
	}
	ans_granges = GRanges_list[[1]]
	if(length(GRanges_list) > 1){
		for(k in 2:length(GRanges_list)){
			ans_granges = c(ans_granges, GRanges_list[[k]])
		}
	}
	return(ans_granges)
}

yyx_sort_GRanges = function(GRanges){
	return(sort(sortSeqlevels(GRanges), ignore.strand=TRUE))
}


# Function to export GenomicRanges object to a BED file with all metadata
exportToBEDWithMetadata <- function(granges, file_path, col_names=FALSE) {
  # Ensure the object is a GenomicRanges object
  if (!inherits(granges, "GenomicRanges")) {
    stop("Input must be a GenomicRanges object.")
  }

  # Convert GenomicRanges to a data frame
  df <- as.data.frame(granges)

  # Adjust the columns to fit the BED format: chromosome, start, end, name, score, strand
  bed_df <- data.frame(chromosome = df$seqnames,
                       start = df$start - 1, # BED format is 0-based
                       end = df$end,
                       name = ".",
                       score = 0,
                       strand = df$strand)

  # Bind the metadata columns
  metadata_df <- as.data.frame(mcols(granges))
  if("name" %in% colnames(metadata_df)){
    bed_df$name = metadata_df$name
    metadata_df$name = NULL
  }
  if("score" %in% colnames(metadata_df)){
    bed_df$score = metadata_df$score
    metadata_df$score = NULL
  }
  if(ncol(metadata_df) > 0){
    bed_df <- cbind(bed_df, metadata_df)
  }

  # Write to a BED file
  write_tsv(bed_df, file_path, col_names = col_names)
}

exportToBedgraph <- function(granges, file_path, col_names=FALSE) {
  # Ensure the object is a GenomicRanges object
  if (!inherits(granges, "GenomicRanges")) {
    stop("Input must be a GenomicRanges object.")
  }

  # Convert GenomicRanges to a data frame
  df <- as.data.frame(granges)

  # Adjust the columns to fit the bdg format: chromosome, start, end, score
  bdg_df <- data.frame(chromosome = df$seqnames,
                       start = df$start - 1, # BED format is 0-based
                       end = df$end,
                       score = 0)

  # Bind the metadata columns
  metadata_df <- as.data.frame(mcols(granges))
  if("score" %in% colnames(metadata_df)){
    bdg_df$score = metadata_df$score
    metadata_df$score = NULL
  }

  # Write to a BED file
  write_tsv(bdg_df, file_path, col_names = col_names)
}


yyx_read_bed_or_bdg_into_GRanges = function(filename){
	skip = 0
	if(endsWith(filename, "narrowPeak")){
		skip = 1   # skip headline of narrowPeak file
	}
	input = read_tsv(filename, col_names=FALSE, skip=skip)
#	input = input[!grepl("^track", input[[1]]),]   # remove headline of narrowPeak file
	ans = GRanges(
		seqnames = Rle(input[[1]]),
		ranges = IRanges(start=input[[2]]+1, end=input[[3]]),
		strand = "*"
	)
	if(ncol(input) >= 4){
		if(ncol(input) >= 6){
			strand(ans) = sub("[.]", "*", input[[6]])
		}
		if(ncol(input) >= 5){
			if(is.numeric(input[[5]])){
				ans$score = input[[5]]
				ans$name = input[[4]]
			}
		}else{
			if(is.numeric(input[[4]])){
			ans$score = input[[4]]
			}else{
				ans$name = input[[4]]
			}
		}
	}
	ans
}


yyx_GRanges_apply_per_chr <- function(granges, func) {
  if (!inherits(granges, "GenomicRanges")) {
    stop("Input must be a GenomicRanges object.")
  }
  lapply(split(granges, seqnames(granges)), func)
}

yyx_GRanges_apply_per_strand <- function(granges, func) {
  if (!inherits(granges, "GenomicRanges")) {
    stop("Input must be a GenomicRanges object.")
  }
  lapply(split(granges, strand(granges)), func)
}

# Function to reduce overlapping ranges but not adjacent ranges
#reduceOverlappingOnly <- function(granges) {
#  if (!inherits(granges, "GenomicRanges")) {
#    stop("Input must be a GenomicRanges object.")
#  }
#  
#  return(yyx_unlist_a_list_of_GRanges(yyx_GRanges_apply_per_chr(granges, function(now_chr_granges){
#	return(yyx_unlist_a_list_of_GRanges(yyx_GRanges_apply_per_strand(now_chr_granges, function(x){
#		if (length(x) <= 1) return(x)
#
#		x = sort(x)
#
#		start_positions <- start(x)
#		end_positions <- end(x)
#		keep <- rep(TRUE, length(x))
#
#		for (i in rev(seq_along(x)[-1])) {
#			if (start_positions[i] <= end_positions[i - 1]) {
#				end_positions[i - 1] <- max(end_positions[i - 1], end_positions[i])
#				keep[i] <- FALSE
#			}
#		}
#
#		GRanges(seqnames(x)[keep],
#			IRanges(start_positions[keep], end_positions[keep]),
#			strand = strand(x)[keep])
#	})))
#  })))
#}
reduceOverlappingOnly <- function(granges) {
	reduce(granges, min.gapwidth=0)
}

yyx_multi_GRanges_overlap = function(granges_list){
	reduced = reduceOverlappingOnly(yyx_unlist_a_list_of_GRanges(granges_list))
	reduced_anno = character(length(reduced))
	reduced_scores = numeric(length(reduced))
	for(k in 1:length(granges_list)){
		overlaps = findOverlaps(reduced, granges_list[[k]])
		for (i in seq_along(reduced)) {
			overlapping_original <- subjectHits(overlaps[queryHits(overlaps) == i])
			if(length(overlapping_original) > 0){
				reduced_anno[i] = reduced_anno[i] %.% ",rep" %.% k %.% ":" %.% paste(collapse=",", overlapping_original)
				reduced_scores[i] = reduced_scores[i] + 1
			}
		}
	}
	reduced$name = sub("^,", "", reduced_anno)
	reduced$score = reduced_scores
	reduced
}


input_granges_list = list()

for(k in 1:length(input_filenames)){
	cat("Now read bed or bdg " %.% input_filenames[k] %.% " ...\n")
	print(system.time({
	input_granges_list[[k]] = yyx_read_bed_or_bdg_into_GRanges(input_filenames[k])
	}))
}

cat("Now merge overlapping ...\n")
print(system.time({
merged_granges = yyx_multi_GRanges_overlap(input_granges_list)
}))

exportToBEDWithMetadata(merged_granges, output_prefix %.% ".bed")

exportToBedgraph(merged_granges, output_prefix %.% ".bdg")

