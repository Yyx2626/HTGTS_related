#### Usage: Rscript this.r <input.tlx> <output_preifx> <ref.fa>
####   [motif (default: CATG)] [max_dist (default: 10)] [motif_sites.bed]

#### 2023-01-02 morning, Fred asked to develope a new peak calling program that can identify peaks out of local background noise level, not thresholding on signal height.

args = commandArgs(TRUE)
if(length(args) < 3){
	stop("not enough command-line parameters
Usage: Rscript this.r <input.tlx> <output_preifx> <ref.fa>
  [motif (default: CATG)] [max_dist (default: 10)] [motif_sites.bed]")
}

input_tlx_filename = args[1]
output_prefix = args[2]
ref_fa_filename = args[3]
motif = "CATG"
if(length(args) >= 4){
	motif = args[4]
}
max_dist = 10
if(length(args) >= 5){
	max_dist = as.numeric(args[5])
}
motif_sites_bed_filename = NULL
if(length(args) >= 6){
	motif_sites_bed_filename = args[6]
}



options(stringsAsFactors=FALSE)

library(tidyverse)
library(GenomicRanges)
library(Biostrings)



`%.%` = function(x,y) paste0(x,y)

# Function to read .fai or chromSize file and create a seqinfo object for GenomicRanges
readSeqInfo <- function(file_path, genome_name = NA) {
  # Read the file
  data <- read.table(file_path, header = FALSE)
  
  # Check the file type (.fai or chromSize) to determine the position of columns
  if (ncol(data) >= 2) {
    # .fai files usually have at least 5 columns, chromSize files have 2 columns
    seqnames_col <- 1
    seqlengths_col <- ifelse(ncol(data) >= 5, 2, 2)
  } else {
    stop("Invalid file format. The file should be a .fai or chromSize file.")
  }

  # Create Seqinfo object
  seqinfo <- Seqinfo(seqnames = data[, seqnames_col],
                     seqlengths = data[, seqlengths_col],
                     isCircular = rep(FALSE, length(data[, seqnames_col])),
                     genome = genome_name
                )
  
  return(seqinfo)
}

# Function to find motif sites using regular expressions and save to GenomicRanges
findMotifSites <- function(fasta_file, motif) {
  # Read the FASTA file
  sequences <- readDNAStringSet(fasta_file)

  # Initialize a list to store GRanges for each sequence
  all_sites <- list()

  # Iterate over each sequence
  for (i in seq_along(sequences)) {
    seqname <- names(sequences)[i]
    sequence <- sequences[[i]]

    # Find positions of the motif using pattern matching
    match_positions <- matchPattern(DNAString(motif), sequence)

    # Create GRanges object and add to the list
    all_sites[[seqname]] <- GRanges(
      seqnames = Rle(seqname),
      ranges = ranges(match_positions),
      strand = Rle("*")
    )
  }

  # Combine all GRanges objects into one
  return(yyx_unlist_a_list_of_GRanges(all_sites))
}


yyx_unlist_a_list_of_GRanges = function(GRanges_list){
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


yyx_read_tlx_junction_into_GRanges = function(filename){
	input = read_tsv(filename)
	GRanges(
		seqnames = Rle(input$Rname),
		ranges = IRanges(start=input$Junction, end=input$Junction),
		strand = Rle(ifelse(input$Strand > 0, "+" , "-"))
	)
}

yyx_read_bed_into_GRanges = function(filename){
	input = read_tsv(filename, col_names=FALSE)
	ans = GRanges(
		seqnames = Rle(input[[1]]),
		ranges = IRanges(start=input[[2]]+1, end=input[[3]]),
		strand = "*"
	)
	if(ncol(input) >= 4){
		ans$name = input[[4]]
	}
	if(ncol(input) >= 5){
		ans$score = input[[5]]
	}
	if(ncol(input) >= 6){
		strand(ans) = sub("[.]", "*", input[[6]])
	}
	ans
}


# Function to reduce overlapping ranges but not adjacent ranges
#reduceOverlappingOnly <- function(granges) {
#  if (!inherits(granges, "GenomicRanges")) {
#    stop("Input must be a GenomicRanges object.")
#  }
#
#  ans <- lapply(split(granges, seqnames(granges)), function(now_chr_granges) {
#   reduced_list = lapply(split(now_chr_granges, strand(now_chr_granges)), function(x) {
#    if (length(x) <= 1) return(x)
#    
#    x = sort(x)
#
#    start_positions <- start(x)
#    end_positions <- end(x)
#    keep <- rep(TRUE, length(x))
#  
#    for (i in rev(seq_along(x)[-1])) {
#      if (start_positions[i] <= end_positions[i - 1]) {
#        end_positions[i - 1] <- max(end_positions[i - 1], end_positions[i])
#        keep[i] <- FALSE
#      }
#    }
#  
#    reduced <- GRanges(seqnames(x)[keep],
#                       IRanges(start_positions[keep], end_positions[keep]),
#                       strand = strand(x)[keep])
#
#    return(reduced)
#   })
#   return(yyx_unlist_a_list_of_GRanges(reduced_list))
#  })
#
#  return(yyx_unlist_a_list_of_GRanges(ans))
#}
reduceOverlappingOnly <- function(granges) {
	reduce(granges, min.gapwidth=0)
}


# Function to summarize a GenomicRanges object with summed scores or frequencies
summarizeGRangesByScoresOrFreq <- function(granges) {
  # Check if the input is a GenomicRanges object
  if (!inherits(granges, "GenomicRanges")) {
    stop("Input must be a GenomicRanges object.")
  }

  # Check if 'score' column exists, if not, create one and set all values to 1
  if (!"score" %in% names(mcols(granges))) {
    granges$score <- 1
  }

  # Reduce the GenomicRanges object

  # Sum the scores or count frequencies of the reduced ranges
  ans <- lapply(split(granges, seqnames(granges)), function(now_chr_granges) {
   reduced_list = lapply(split(now_chr_granges, strand(now_chr_granges)), function(x) {
    if (length(x) <= 1) return(x)
    reduced_granges <- reduceOverlappingOnly(x)
    
    # Add the summed scores to the reduced GRanges object
    mcols(reduced_granges)$score <- sapply(seq_along(reduced_granges), function(i) {
      sum(x$score[which(start(x) <= end(reduced_granges)[i] & 
                              end(x) >= start(reduced_granges)[i])])
    })
    return(reduced_granges)
   })
   return(yyx_unlist_a_list_of_GRanges(reduced_list))
  })

  return(yyx_unlist_a_list_of_GRanges(ans))
}

# Function to find the nearest subject block for each query block
findNearestSubjectBlock <- function(query, subject) {
  # Check if the inputs are GenomicRanges objects
  if (!inherits(query, "GenomicRanges")) {
    stop("Query must be a GenomicRanges object.")
  }
  if (!inherits(subject, "GenomicRanges")) {
    stop("Subject must be a GenomicRanges object.")
  }

  # Find the nearest subject block for each query block
  nearest_hits <- distanceToNearest(query, subject)

  # Add the subject info and distance to the query object as metadata columns
  mcols(query)$nearest_subject_seqname <- seqnames(subject)[subjectHits(nearest_hits)]
  mcols(query)$nearest_subject_start <- start(subject)[subjectHits(nearest_hits)]
  mcols(query)$nearest_subject_end <- end(subject)[subjectHits(nearest_hits)]
  mcols(query)$distance <- mcols(nearest_hits)$distance

  return(query)
}




seqinfo_obj <- readSeqInfo(ref_fa_filename %.% ".fai")

if(is.null(motif_sites_bed_filename) || !file.exists(motif_sites_bed_filename)){
	cat("Now searching for motif " %.% motif %.% " in ref " %.% ref_fa_filename %.% " ...\n")
	print(system.time({
	CATG_sites <- findMotifSites(ref_fa_filename, motif)
	}))
	CATG_sites = yyx_unlist_a_list_of_GRanges(CATG_sites)
	CATG_sites = yyx_sort_GRanges(CATG_sites)
	if(!is.null(motif_sites_bed_filename)){
		cat("Now output motif sites to " %.% motif_sites_bed_filename %.% " ...\n")
		print(system.time({
		exportToBEDWithMetadata(CATG_sites, motif_sites_bed_filename)
		}))
	}
}else{
	cat("Now read motif sites from " %.% motif_sites_bed_filename %.% " ...\n")
	print(system.time({
	CATG_sites = yyx_read_bed_into_GRanges(motif_sites_bed_filename)
	}))
}

cat("Now read tlx " %.% input_tlx_filename %.% " ...\n")
print(system.time({
tlx_granges = yyx_read_tlx_junction_into_GRanges(input_tlx_filename)
}))

cat("Now count frequency on bp ...\n")
print(system.time({
summarized_tlx_granges <- summarizeGRangesByScoresOrFreq(tlx_granges)
}))

cat("Now find nearest motif ...\n")
print(system.time({
nearest_CATG_granges <- findNearestSubjectBlock(summarized_tlx_granges, CATG_sites)
}))

cat("Now filter on motif distance <= " %.% max_dist %.% " ...\n")
print(system.time({
filtered_CATG_granges = nearest_CATG_granges[abs(nearest_CATG_granges$distance) <= max_dist]
}))
#seqnames(filtered_CATG_granges) = filtered_CATG_granges$nearest_subject_seqname
ranges(filtered_CATG_granges) = IRanges(start=filtered_CATG_granges$nearest_subject_start, end=filtered_CATG_granges$nearest_subject_end)
# Remove the nearest subject information from mcols
mcols(filtered_CATG_granges)$nearest_subject_seqname <- NULL
mcols(filtered_CATG_granges)$nearest_subject_start <- NULL
mcols(filtered_CATG_granges)$nearest_subject_end <- NULL

filtered_CATG_granges = yyx_sort_GRanges(filtered_CATG_granges)

cat("Now sum up frequency (stranded) and output to " %.% output_prefix %.% ".bed ...\n")
print(system.time({
summarized_filtered_CATG_granges = summarizeGRangesByScoresOrFreq(filtered_CATG_granges)
}))   # 5s

exportToBEDWithMetadata(yyx_sort_GRanges(summarized_filtered_CATG_granges), output_prefix %.% ".bed")

cat("Now sum up frequency (unstranded) and output to " %.% output_prefix %.% ".bdg ...\n")
unstranded_filtered_CATG_granges = summarized_filtered_CATG_granges
strand(unstranded_filtered_CATG_granges) = "*"
print(system.time({
unstranded_filtered_CATG_granges = summarizeGRangesByScoresOrFreq(unstranded_filtered_CATG_granges)
}))   # 5s

exportToBedgraph(unstranded_filtered_CATG_granges, output_prefix %.% ".bdg")



