#### Usage: Rscript this.r <output_prefix> <should_abs_score> <query.bed|bdg> <anno.bed|bdg|bw>

#### 2024-01-05 morning, this.r sumScores() is too slow for adding pos.bw and neg.bw, so I remove this function
####  Please use external yyx_overlapping_bed_to_nonoverlapping_bdg.20240105.pl sum instead

####   I plan to implement  yyx_anno_underlying_features.r  to quantify one anno feature for each peak region:
####      (anno) count, density, highest signal, average signal

args = commandArgs(TRUE)
argIdx = 3
if(length(args) < argIdx){
	stop("not enough command-line parameters
Usage: Rscript this.r <output_prefix> <should_abs_score> <query.bed|bdg> <anno.bed|bdg|bw>
")
}


`%.%` = function(x,y) paste0(x,y)

# Define the function
isTrueOrFalse <- function(inputString) {
    # Define lists of strings that mean TRUE and FALSE
    trueStrings <- c("true", "yes", "1", "t", "y")
    falseStrings <- c("false", "no", "0", "f", "n")

    # Convert the input string to lower case for case-insensitive comparison
    inputString <- tolower(inputString)

    if (inputString %in% trueStrings) {   # Check if the string means TRUE
        return(TRUE)
    }else if (inputString %in% falseStrings) {
        return(FALSE)
    }else {
		stop("Cannot determine True or False on " %.% inputString)
        return(NA)   # or you could choose to return a message like "Invalid input"
    }
}


output_prefix = args[1]
should_abs_score = isTrueOrFalse(args[2])
query_bed_filename = args[3]
anno_filename = args[4]

#argIdx = argIdx + 1
#query_slop_bp = 0
#if(length(args) >= argIdx){
#	query_slop_bp = as.numeric(args[argIdx])
#}



options(stringsAsFactors=FALSE)

library(tidyverse)
library(GenomicRanges)



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

yyx_read_query_bed_into_GRanges = function(filename, skip=0){
	input = read_tsv(filename, col_names=FALSE, skip=skip)
	ans = GRanges(
		seqnames = Rle(input[[1]]),
		ranges = IRanges(start=input[[2]]+1, end=input[[3]]),
		strand = "*"
	)
	if(ncol(input) >= 4){
		for(colidx in 4:ncol(input)){
			mcols(ans)[["V" %.% colidx]] = input[[colidx]]
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
reduceOverlappingOnly <- function(granges) {
	reduce(granges, min.gapwidth=0)
}


library(GenomicRanges)

# Define the function to sum the scores
sumScores <- function(grA, grB) {
	# Check if both GRanges have a 'score' column
	if (!"score" %in% colnames(mcols(grA)) || !"score" %in% colnames(mcols(grB))) {
		stop("Both GRanges must have a 'score' column.")
	}

	# Find all overlapping ranges
	overlaps <- findOverlaps(grA, grB)
	
	# Combine the GRanges, keeping all ranges
	combinedGR <- c(grA, grB)

	# Sum the signal values for overlapping ranges
	if (length(overlaps) > 0) {
		# Create a new GRanges object to store the summed scores
		summedGR <- GRanges()

		# Loop through overlaps
		for (i in seq_along(queryHits(overlaps))) {
			# Find the range of the overlap
			overlapGRanges = intersect(grA[queryHits(overlaps)[i]], grB[subjectHits(overlaps)[i]])

			# Sum the scores
			overlapGRanges$score <- sum(c(mcols(grA)[queryHits(overlaps)[i], "score"], 
							   mcols(grB)[subjectHits(overlaps)[i], "score"]))

			# Append to the new GRanges object
			summedGR <- c(summedGR, overlapGRanges)
		}

		# Combine with non-overlapping ranges
		nonOverlapsA <- setdiff(grA, summedGR)
		nonOverlapsA$score = NA
		overlapsA = findOverlaps(nonOverlapsA, grA)
		if (length(overlapsA) > 0) {
			for (i in seq_along(queryHits(overlapsA))) {
				nonOverlapsA$score[queryHits(overlapsA)[i]] = grA$score[subjectHits(overlapsA)[i]]
			}
		}
		nonOverlapsB <- setdiff(grB, summedGR)
		nonOverlapsB$score = NA
		overlapsB = findOverlaps(nonOverlapsB, grB)
		if (length(overlapsB) > 0) {
			for (i in seq_along(queryHits(overlapsB))) {
				nonOverlapsB$score[queryHits(overlapsB)[i]] = grB$score[subjectHits(overlapsB)[i]]
			}
		}
		combinedGR <- c(summedGR, nonOverlapsA, nonOverlapsB)
	}

	# Sort and return
	return(yyx_sort_GRanges(combinedGR))
}

# Example usage
# Assuming grA and grB are defined GRanges objects with a 'score' column
# grA = GRanges("chr", IRanges(seq(1,7,by=3), seq(3,9,by=3)))
# grA$score = 1:3
# grB = GRanges("chr", IRanges(c(2,6), c(5,10)))
# grB$score = 2:1
# sumScores(grA, grB)


#anno_filename = anno_filenames[1]
anno_raw_granges = NULL
anno_coverage_granges = NULL
if(grepl(".bed$", anno_filename, ignore.case=TRUE)){
	cat("Now read anno bed " %.% anno_filename %.% " ...\n")
	anno_raw_granges = yyx_read_bed_or_bdg_into_GRanges(anno_filename)
#	if(length(anno_filenames) > 1){
#		for(k in 2:length(anno_filenames)){
#			anno_filename = anno_filenames[k]
#			if(grepl(".bed$", anno_filename, ignore.case=TRUE)){
#				cat("Now read more anno bed " %.% anno_filename %.% " ...\n")
#				more_raw_granges = yyx_read_bed_or_bdg_into_GRanges(anno_filename)
#				anno_raw_granges = c(anno_raw_granges, more_raw_granges)
#			}else{
#				stop("Unmatched format according to file extension of " %.% anno_filename)
#			}
#		}
#	}
	anno_coverage_granges = as(coverage(anno_raw_granges), "GRanges")
	anno_coverage_granges = anno_coverage_granges[anno_coverage_granges$score > 0]
}else if(grepl(".bdg$", anno_filename, ignore.case=TRUE) || grepl(".bedgraph$", anno_filename, ignore.case=TRUE)){
	cat("Now read anno bdg " %.% anno_filename %.% " ...\n")
	anno_coverage_granges = yyx_read_bed_or_bdg_into_GRanges(anno_filename)
	if(should_abs_score){
		anno_coverage_granges$score = abs(anno_coverage_granges$score)
	}
#	if(length(anno_filenames) > 1){
#		for(k in 2:length(anno_filenames)){
#			anno_filename = anno_filenames[k]
#			if(grepl(".bdg$", anno_filename, ignore.case=TRUE) || grepl(".bedgraph$", anno_filename, ignore.case=TRUE)){
#				cat("Now read more anno bdg " %.% anno_filename %.% " ...\n")
#				more_coverage_granges = yyx_read_bed_or_bdg_into_GRanges(anno_filename)
#				if(should_abs_score){
#					more_coverage_granges$score = abs(more_coverage_granges$score)
#				}
#				anno_coverage_granges = sumScores(anno_coverage_granges, more_coverage_granges)
#			}else{
#				stop("Unmatched format according to file extension of " %.% anno_filename)
#			}
#		}
#	}
}else if(grepl(".bw$", anno_filename, ignore.case=TRUE) || grepl(".bigwig$", anno_filename, ignore.case=TRUE)){
	require(rtracklayer)
	cat("Now read anno bw " %.% anno_filename %.% " ...\n")
	anno_coverage_granges = import(anno_filename)
	if(should_abs_score){
		anno_coverage_granges$score = abs(anno_coverage_granges$score)
	}
#	if(length(anno_filenames) > 1){
#		for(k in 2:length(anno_filenames)){
#			anno_filename = anno_filenames[k]
#			if(grepl(".bw$", anno_filename, ignore.case=TRUE) || grepl(".bigwig$", anno_filename, ignore.case=TRUE)){
#				cat("Now read more anno bw " %.% anno_filename %.% " ...\n")
#				more_coverage_granges = import(anno_filename)
#				if(should_abs_score){
#					more_coverage_granges$score = abs(more_coverage_granges$score)
#				}
#				anno_coverage_granges = sumScores(anno_coverage_granges, more_coverage_granges)
#			}else{
#				stop("Unmatched format according to file extension of " %.% anno_filename)
#			}
#		}
#	}
}else{
	stop("Unrecognized format according to file extension of " %.% anno_filename)
}

anno_coverage_granges = yyx_sort_GRanges(anno_coverage_granges)


cat("Now read query " %.% query_bed_filename %.% " ...\n")
query_granges = yyx_read_query_bed_into_GRanges(query_bed_filename)
if(grepl(".bed$", anno_filename, ignore.case=TRUE)){
	query_granges$anno_count = countOverlaps(query_granges, anno_raw_granges)
	query_granges$anno_density = query_granges$anno_count / width(query_granges)
}
query_granges$anno_max = NA
query_granges$anno_avg = NA

overlaps <- findOverlaps(query_granges, anno_coverage_granges)
# Iterate over each range in queryGR
for (i in seq_along(query_granges)) {
	if(sum(queryHits(overlaps) == i) > 0){
		overlap_anno_coverage_granges = anno_coverage_granges[subjectHits(overlaps[queryHits(overlaps) == i])]
#		intersect_granges = NULL
#		for(j in seq_along(overlap_anno_coverage_granges)){   # very slow ...
#			if(is.null(intersect_granges)){
#				intersect_granges = intersect(query_granges[i], overlap_anno_coverage_granges[j])
#			}else{
#				intersect_granges = c(intersect_granges, intersect(query_granges[i], overlap_anno_coverage_granges[j]))
#			}
#		}
#		stopifnot(length(intersect_granges)==length(overlap_anno_coverage_granges))
#		overlap_lengths = width(intersect_granges)
		# the code above is very slow
		# Because query_granges[i] is a continuous block, anno_coverage_granges is sorted
		# so only the first and last overlap_anno_coverage_granges may need to shrink width
		if(start(overlap_anno_coverage_granges[1]) < start(query_granges[i])){
			start(overlap_anno_coverage_granges[1]) = start(query_granges[i])
		}
		if(end(overlap_anno_coverage_granges[length(overlap_anno_coverage_granges)]) > end(query_granges[i])){
			end(overlap_anno_coverage_granges[length(overlap_anno_coverage_granges)]) = end(query_granges[i])
		}
		overlap_lengths = width(overlap_anno_coverage_granges)
		query_granges$anno_max[i] = max(overlap_anno_coverage_granges$score)
		query_granges$anno_avg[i] = sum(overlap_anno_coverage_granges$score * overlap_lengths) / width(query_granges)
	}
}

exportToBEDWithMetadata(query_granges, output_prefix %.% ".bed")


