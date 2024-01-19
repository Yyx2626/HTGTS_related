#### Usage: Rscript this.r <input.CATG_dist10_filtered.bdg> <output_preifx>
####   [local_idx_window_size (default: 100)]
####   [adj_pVal_method (default: bonferroni)]
####   [summit_adj_pVal_threshold (default: 0.05)]
####   [foot_pVal_threshold (default: 0.05)]
####   [max_output_log10_pVal (default: 40)]


#### 2024-01-04 night, I also need to output summit height

#### 2024-01-04 afternoon, I output best summit coordinate and avg summit coordinate

#### 2023-01-02 morning, Fred asked to develope a new peak calling program that can identify peaks out of local background noise level, not thresholding on signal height.

args = commandArgs(TRUE)
argIdx = 2
if(length(args) < argIdx){
	stop("not enough command-line parameters
Usage: Rscript this.r <input.CATG_dist10_filtered.bdg> <output_preifx>
  [local_idx_window_size (default: 100)]
  [adj_pVal_method (default: bonferroni)]
  [summit_adj_pVal_threshold (default: 0.05)]
  [foot_pVal_threshold (default: 0.05)]
  [max_output_log10_pVal (default: 40)]")
}

input_bdg_filename = args[1]
output_prefix = args[2]

argIdx = argIdx + 1
local_idx_window_size = 100
if(length(args) >= argIdx){
	local_idx_window_size = as.numeric(args[argIdx])
}

argIdx = argIdx + 1
adj_pVal_method = "bonferroni"
if(length(args) >= argIdx){
	adj_pVal_method = args[argIdx]
}

argIdx = argIdx + 1
summit_adj_pVal_threshold = 0.05
if(length(args) >= argIdx){
	summit_adj_pVal_threshold = as.numeric(args[argIdx])
}

argIdx = argIdx + 1
foot_pVal_threshold = 0.05
if(length(args) >= argIdx){
	foot_pVal_threshold = as.numeric(args[argIdx])
}

argIdx = argIdx + 1
max_output_log10_pVal = 40
if(length(args) >= argIdx){
	max_output_log10_pVal = as.numeric(args[argIdx])
}



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

yyx_read_bdg_into_GRanges = function(filename){
	input = read_tsv(filename, col_names=FALSE)
	ans = GRanges(
		seqnames = Rle(input[[1]]),
		ranges = IRanges(start=input[[2]]+1, end=input[[3]]),
		strand = "*"
	)
	if(ncol(input) >= 4){
		ans$score = input[[4]]
	}
	ans
}


# Assuming your data is stored in a vector named 'data'
local_poisson_test = function(data, window_size=100){
	# Function: Calculate local median
	local_median <- function(data, index, window_size) {
	  start <- max(1, index - window_size %/% 2)
	  end <- min(length(data), index + window_size %/% 2)
	  median(data[start:end])
	}

	# Calculate local background noise level for each point
	local_medians <- sapply(1:length(data), function(i) local_median(data, i, window_size))

	# Poisson test
	pVals <- sapply(1:length(data), function(i) {
	  ppois(data[i], lambda=local_medians[i], lower.tail = FALSE)
	})
	
	pVals
}


define_peak_summit_and_feet = function(adj_pVal_vec, summit_adj_pVal_threshold, pVal_vec, foot_pVal_threshold=0.05){
	stopifnot(length(adj_pVal_vec)==length(pVal_vec))
	L = length(pVal_vec)
	peak_DF = data.frame(summit_idx = which(adj_pVal_vec < summit_adj_pVal_threshold))
	if(nrow(peak_DF) == 0){
		return(NULL)
	}
	peak_DF$left_foot_idx = NA
	peak_DF$right_foot_idx = NA
	for(i in 1:nrow(peak_DF)){
		idx = peak_DF$summit_idx[i]
		last_pVal = pVal_vec[idx]
		while(idx >= 1){
			if(pVal_vec[idx] < last_pVal && last_pVal >= foot_pVal_threshold){
				break
			}
			last_pVal = pVal_vec[idx]
			idx = idx - 1
		}
		if(idx < 1){ idx = 1 }
		peak_DF$left_foot_idx[i] = idx

		idx = peak_DF$summit_idx[i]
		last_pVal = pVal_vec[idx]
		while(idx <= L){
			if(pVal_vec[idx] < last_pVal && last_pVal >= foot_pVal_threshold){
				break
			}
			last_pVal = pVal_vec[idx]
			idx = idx + 1
		}
		if(idx > L){ idx = L }
		peak_DF$right_foot_idx[i] = idx
	}
	peak_DF
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

yyx_merge_overlapping_peaks = function(peak_granges, bdg_granges){
	return(yyx_unlist_a_list_of_GRanges(yyx_GRanges_apply_per_chr(peak_granges, function(x){
		x = sort(x)
#		start_positions <- start(x)
#		end_positions <- end(x)
#		keep <- rep(TRUE, length(x))
#
#		for (i in rev(seq_along(x)[-1])) {
#			if (start_positions[i] <= end_positions[i - 1]) {
#				end_positions[i - 1] <- max(end_positions[i - 1], end_positions[i])
#				x$score[i - 1] <- max(x$score[i-1], x$score[i])
#				x$pVal[i - 1] <- min(x$pVal[i-1], x$pVal[i])
#				keep[i] <- FALSE
#			}
#		}
#
#		reduced <- GRanges(seqnames(x)[keep],
#			IRanges(start_positions[keep], end_positions[keep]),
#			strand = strand(x)[keep])
		reduced = reduce(x, min.gapwidth=0)
		reduced_scores = numeric(length(reduced))
		overlaps = findOverlaps(reduced, bdg_granges)
		for (i in seq_along(reduced)) {
			overlapping_bdg <- subjectHits(overlaps[queryHits(overlaps) == i])
			reduced_scores[i] <- sum(bdg_granges$score[overlapping_bdg])
		}
		reduced$score = reduced_scores
#		reduced$pVal = x$pVal[keep]
		reduced_pVals = rep(NA, length(reduced))
		## 2023-01-04 added summit_start summit_end summit_score(height)
		reduced_summit_start = rep(NA, length(reduced))
		reduced_summit_end = rep(NA, length(reduced))
		reduced_summit_score = rep(NA, length(reduced))
		overlaps = findOverlaps(reduced, x)
		for (i in seq_along(reduced)) {
			overlapping_x <- subjectHits(overlaps[queryHits(overlaps) == i])
			reduced_pVals[i] = min(x$pVal[overlapping_x])
			tmp_idx = which.min(x$pVal[overlapping_x])
			reduced_summit_start[i] = x$summit_start[overlapping_x][tmp_idx]
			reduced_summit_end[i] = x$summit_end[overlapping_x][tmp_idx]
			reduced_summit_score[i] = x$summit_score[overlapping_x][tmp_idx]
		}
		reduced$pVal = reduced_pVals
		reduced$summit_start = reduced_summit_start
		reduced$summit_end = reduced_summit_end
		reduced$summit_score = reduced_summit_score

		return(reduced)
	})))
}


cat("Now read CATG_filtered bdg " %.% input_bdg_filename %.% " ...\n")
print(system.time({
bdg_granges = yyx_read_bdg_into_GRanges(input_bdg_filename)
}))
cat("In total, " %.% length(bdg_granges) %.% " cutting sites have signals\n")

cat("Now run local Poisson test with window_size = " %.% local_idx_window_size %.% " ...\n")
print(system.time({
bdg_granges = yyx_unlist_a_list_of_GRanges(yyx_GRanges_apply_per_chr(bdg_granges, function(x){
	x = sort(x)
	x$pVal = local_poisson_test(x$score, window_size=local_idx_window_size)
	x
}))
}))

cat("Bonferroni correction ...\n")
print(system.time({
bdg_granges$adj_pVal = p.adjust(bdg_granges$pVal, method=adj_pVal_method)
}))

cat("Now define peak summit and feet ...\n")
print(system.time({
peak_granges = yyx_unlist_a_list_of_GRanges(yyx_GRanges_apply_per_chr(bdg_granges, function(x){
#	x = sort(x)
	peak_DF = define_peak_summit_and_feet(x$adj_pVal, summit_adj_pVal_threshold, x$pVal, foot_pVal_threshold)
	if(is.null(peak_DF)){ return(NULL) }
	ans = x[peak_DF$summit_idx]
	left_granges = x[peak_DF$left_foot_idx]
	right_granges = x[peak_DF$right_foot_idx]
	start(ans) = pmin(start(ans), start(left_granges))
	end(ans) = pmax(end(ans), end(right_granges))
	for(i in 1:nrow(peak_DF)){
		ans$score[i] = sum(x$score[peak_DF$left_foot_idx[i]:peak_DF$right_foot_idx[i]])
	}
	## 2023-01-04 added summit_start summit_end summit_score(height); summit defined by adjusted_pVal
	ans$summit_start = start(x[peak_DF$summit_idx])
	ans$summit_end = end(x[peak_DF$summit_idx])
	ans$summit_score = x$score[peak_DF$summit_idx]
	ans
}))
}))

cat("After Bonferroni correction, " %.% length(peak_granges) %.% " cutting sites are significant\n")

exportToBEDWithMetadata(peak_granges, output_prefix %.% ".bed")

exportToBedgraph(peak_granges, output_prefix %.% ".read_count.bdg")

tmp_granges = peak_granges
tmp_granges$score = -log10(tmp_granges$pVal)
tmp_granges$score[tmp_granges$score > max_output_log10_pVal] = max_output_log10_pVal
exportToBedgraph(tmp_granges, output_prefix %.% ".log10_pVal.bdg")


cat("Now merge overlapping peaks ...\n")
merged_peak_granges = yyx_merge_overlapping_peaks(peak_granges, bdg_granges)

exportToBEDWithMetadata(merged_peak_granges, output_prefix %.% ".merged.bed")

exportToBedgraph(merged_peak_granges, output_prefix %.% ".merged.read_count.bdg")

tmp_granges = merged_peak_granges
tmp_granges$score = -log10(tmp_granges$pVal)
tmp_granges$score[tmp_granges$score > max_output_log10_pVal] = max_output_log10_pVal
exportToBedgraph(tmp_granges, output_prefix %.% ".merged.log10_pVal.bdg")


