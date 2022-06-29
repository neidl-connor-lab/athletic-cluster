#!/usr/bin/env Rscript

## setup -------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(argparse))
options(stringsAsFactors=FALSE)

# read in command-line arguments
args <- ArgumentParser()
args$add_argument("-r", "--ref", required=TRUE, help="reference FASTA")
args$add_argument("-v", "--vcf", required=TRUE, help="LoFreq VCF")
args$add_argument("-c", "--cov", required=TRUE, help="coverage TSV")
args$add_argument("-o", "--ofile", required=TRUE, help="GZipped FASTA")
args <- args$parse_args()

# helper functions
mesg <- function(x) {
  paste("[MSG]", x) %>%
    paste0("echo '", ., "'") %>%
    system()
}
err <- function(x) {
  paste("[ERR]", x) %>%
    paste0("echo '", ., "'") %>%
    system()
  quit(status=1)
}

## check inputs and load files ---------------------------------
# SNVs
if(file.exists(args$vcf)) {
  mesg(paste("Valid VCF file:", args$vcf))
  snvs <- read.delim(args$vcf, comment.char="#", header=FALSE,
                     col.names=c("Segment", "Position", "ID", 
                                 "Ref", "Alt", "Quality.score", 
                                 "FILTER", "INFO"),
                     colClasses=c("character", "integer", "character",
                                  "character", "character", "integer", 
                                  "character", "character")) %>%
          # get frequencies
          mutate(Frequency=as.numeric(str_extract(INFO, "(?<=AF=)[0-9\\.]+"))) %>%
          # quality filter & consensus filter
          filter(FILTER=="PASS",
                 Frequency >= 0.5) %>%
          # keep only necessary info
          select(Segment, Position, Ref, Alt)
} else {
  err(paste("Invalid VCF file:", args$vcf))
}

# coverage
if(file.exists(args$cov)) {
  mesg(paste("Valid coverage file:", args$cov))
  covs <- read.delim(args$cov, comment.char="#", header=FALSE,
                     col.names=c("Segment", "Position", "Depth"))
} else {
  err(paste("Invalid coverage file:", args$cov))
}

# ref seq
if(file.exists(args$ref)) {
  mesg(paste("Valid reference:", args$ref))
  refs <- readDNAStringSet(args$ref)
  # make data frame of segments, positions, and ref nucleotides
  refs <- names(refs) %>%
          lapply(function(i) {
            y <- str_extract(i, "[^ ]+")
            x <- refs[i] %>%
                 as.character() %>%
                 strsplit("") %>%
                 unlist()
            data.frame(Position=1:length(x),
                       Segment=y,
                       RefSeq=x,
                       row.names=paste0(y, 1:length(x)))
          }) %>%
          do.call(rbind, .)
}

# output directory and file
odir <- dirname(args$ofile)
if(dir.exists(odir)) {
  paste("Valid output directory:", odir) %>%
    mesg()
} else {
  paste("Creating output directory:", odir) %>%
    mesg()
  dir.create(odir, recursive=TRUE)
}
rm(odir)
mesg(paste("Output file:", args$ofile))

## build consensus -------------------------------------------------------------
# add coverage and refseq to SNV data frame
snvs <- covs %>%
        left_join(snvs, by=c("Segment", "Position")) %>%
        right_join(refs, by=c("Segment", "Position")) %>%
        mutate(Consensus=as.character(NA))
rm(refs, covs)

# set nucleotides WITHOUT Alt alleles to ref
x <- which(is.na(snvs$Alt))
snvs[x, "Consensus"] <- snvs[x, "RefSeq"]

# set nucleotides WITH Alt alleles to Alt
# in case of insertion, Alt with have >1 character but that's okay
# in case of deletion, Alt will only have reference, so no change
x <- which(!is.na(snvs$Alt))
snvs[x, "Consensus"] <- snvs[x, "Alt"]

# set nucleotides with <10X coverage to N
# SNVs were filtered to a minimum depth of 10X by LoFreq
# so no risk of overwriting later when handling deletions
x <- which(snvs$Depth < 10)
snvs[x, "Consensus"] <- "N"
rm(x)

# pull out preliminary consensus 
# includes everything but deletions
# segments & positions are already ordered from coverage table
csus <- unique(snvs$Segment) %>%
        lapply(function(i) {
          snvs %>%
            filter(Segment==i) %>%
            select(Consensus) %>%
            unlist()
        })
names(csus) <- unique(snvs$Segment)

# find deletions and incorporate (if any)
snvs <- filter(snvs, nchar(Ref) > nchar(Alt))
if(dim(snvs)[1] > 0) {
  for(i in 1:dim(snvs)[1]) {
    seg <- as.character(snvs[i, "Segment"])
    pos <- as.integer(snvs[i, "Position"])
    del <- as.character(snvs[i, "Ref"])
    # first character in deletion is the "anchor"
    # so only need to remove nucleotides 2-n
    # step 1: calculate number of nucleotides to remove
    del <- nchar(del)-1
    # step 2: get positions to delete
    pos <- (pos+1):(pos+del)
    # step 3: set positions to NA
    csus[[seg]][pos] <- NA
  }
}
rm(i, seg, pos, del)

# collapse segments and remove 3' and 5' low-coverage ends
csus <- csus %>%
        lapply(function(i) {
          # removed deleted bases
          na.omit(i) %>%
            # collapse into a single string
            paste(collapse="") %>%
            # clip N from ends
            str_extract("(?<=^|N)[ACGT].+[ACGT](?=N|$)")
        })

# add sample ID to segment names
# extracting sample ID from output file name
names(csus) <- basename(args$ofile) %>%
               str_remove("\\.fa\\.gz$") %>%
               paste(names(csus))

# format as DNAString as save
csus <- DNAStringSet(unlist(csus))

# remove deletions and clip 3' and 5' Ns
csus <- csus %>%
        lapply(function(i) {
          na.omit(i) %>%
            paste0(collapse="") %>%
            str_extract("(?<=^|N)[ACGT].+[ACGT](?=N|$)")
        })
# format as DNA string
names(csus) <- paste(args$id, names(csus))
csus <- DNAStringSet(unlist(csus))
writeXStringSet(csus, args$ofile, compress=TRUE)
