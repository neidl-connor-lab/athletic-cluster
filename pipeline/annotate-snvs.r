#!/usr/bin/env Rscript

## setup -------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(argparse))
options(stringsAsFactors=FALSE)

# annotation functions
coding <- function(r, a) {
  # check INDEL type
  if(r$INDEL) {
    # get number of bases added/removed
    # could be in REF or ALT depending if insertion or deletion
    # don't count the original, unchanging base
    n <- max(nchar(r$NT.REF), nchar(r$NT.ALT))-1
    # if n is divisible by 3, no frame shift
    if(n%%3==0) {
      r$INDEL.frameshift <- FALSE
    } else {
      r$INDEL.frameshift <- TRUE
    }
  } else {
    r$INDEL.frameshift <- FALSE
  }
  # get gene start and end
  x <- min(c(a$Start, a$End))
  y <- max(c(a$Start, a$End))
  # pull out nucleotide sequence
  f <- a$NucSeq %>%
       strsplit("") %>%
       unlist()
  # set nucleotide indices
  names(f) <- x:y
  
  # substitute nucleotide
  # no special handling needed if substitution or insertion
  if(nchar(r$NT.REF)==1) {
    x <- as.character(r$Position)
    y <- as.character(r$NT.ALT)
    f[x] <- y
  } else { # deletion
    # NT.ID format is POS-Xyy-X
    # where "y" is the deleted nucleotide(s)
    # get number of y 
    n <- nchar(r$NT.REF)-1
    # get positions to delete
    # pos+1 to pos+n
    # format as character to access the proper positions by name
    x <- as.character((r$Position+1):(r$Position+n))
    # replace deleted bases with NA
    # then remove those positions
    f[x] <- NA
    f <- f[!is.na(f)]
  }
  # format substituted string
  alt <- paste(f, collapse="") %>%
         DNAString()
  # change to reverse complement if on (-) strand
  if(a$Strand=="-") {
    alt <- reverseComplement(alt)
  }
  
  # translate nucleotide sequence
  # then format for comparison
  alt <- translate(alt) %>%
         as.character() %>%
         strsplit("") %>%
         unlist()
  # format reference AA sequence
  ref <- a$AASeq %>%
         strsplit("") %>%
         unlist()
  
  # find first AA change, if any
  x <- which(ref != alt)
  # annotate as synonymous if nothing in x
  if(length(x)==0) {
    r$Type <- "Synonymous"
    r$AA.ID <- paste0(r$Gene, " ", r$NT.REF, r$Position, r$NT.ALT)
    # fill in missing columns with NA
    x$AA.Position <- as.integer(NA)
    x$AA.REF <- as.character(NA)
    x$AA.ALT <- as.character(NA)
  } else { # annotate as nonsynonymous if ≥ 1 AA changes
    # get AA position
    x <- x[1]
    # fill in columns
    r$Type <- "Nonsynonymous"
    r$AA.Position <- x
    r$AA.REF <- as.character(ref[x])
    r$AA.ALT <- as.character(alt[x])
    r$AA.ID <- paste0(r$Gene, " ", r$AA.REF, r$AA.Position, r$AA.ALT)
  }
  
  # return final annotation
  return(r)
}
annot <- function(v, a) {
  # initiate return data.frame
  r <- data.frame(Segment=as.character(v["CHROM"]),
                  Position=as.integer(v["Position"]),
                  NT.ID=as.character(v["ID"]),
                  NT.REF=as.character(v["REF"]),
                  NT.ALT=as.character(v["ALT"]),
                  Quality.score=as.integer(v["Quality.score"]),
                  INDEL=as.logical(v["INDEL"]),
                  Read.depth=as.integer(v["Read.depth"]),
                  Frequency=as.numeric(v["Frequency"]),
                  Strand.bias=as.integer(v["Strand.bias"]),
                  Percent=as.numeric(v["Percent"]))
  # set INDEL type
  if(nchar(r$NT.REF) > 1) {
    r$INDEL.type <- "Deletion"
  } else if(nchar(r$NT.ALT) > 1) {
    r$INDEL.type <- "Insertion"
  } else {
    r$INDEL.type <- "None"
    r$INDEL.frameshift <- FALSE
  }
  # get gene if coding
  a <- a %>%
       filter(Genome == r$Segment,
              Start <= r$Position,
              End >= r$Position)
  # check for noncoding/coding/overlap
  # annotate if coding
  if(dim(a)[1]==0) { # noncoding
    r$Type <- "Noncoding"
    r$Gene <- "None"
    r$AA.ID <- paste0("Noncoding ", r$NT.REF, r$Position, r$NT.ALT)
    # fill in missing return values
    r$AA.Position <- as.integer(NA)
    r$AA.REF <- as.character(NA)
    r$AA.ALT <- as.character(NA)
  } else if(dim(a)[1]==1) { # coding but no ORF overlap
    r$Type <- "Coding"
    r$Gene <- as.character(a$Gene)
    # annotate
    r <- coding(r, a)
  } else { # coding, overlapping ORFs
    r$Type <- "Coding"
    # loop over ORFs to annotate
    r <- 1:dim(a)[1] %>%
         lapply(function(i) {
           x <- r
           y <- a[i, ]
           x$Gene <- y$Gene
           coding(x, y)
         }) %>%
         plyr::join_all(type="full")
  }
  return(r)
}
mesg <- function(...) {
  x <- paste("echo -e '[MSG]", ..., "'")
  system(x)
}
err <- function(...) {
  x <- paste("echo -e '[ERR]", ..., "'")
  system(x)
  q(status=1)
}

## inputs ------------------------------------------------------
args <- ArgumentParser()
args$add_argument("-c", "--csv", required=TRUE, 
                  help="CSV annotation file")
args$add_argument("-v", "--vcf", required=TRUE, 
                  help="LoFreq VCF file")
args$add_argument("-o", "--ofile", required=TRUE, 
                  help="output CSV file")
args <- args$parse_args()

# read in annotation CSV
if (file.exists(args$csv)) {
  ann <- read.csv(args$csv)
} else {
  err("Invalid annotation file:", args$csv)
}

# read in VCF
if (file.exists(args$vcf)) {
  vcf <- read.delim(args$vcf, comment.char="#", header=FALSE,
                    col.names=c("CHROM", "Position", "ID", 
                                "REF", "ALT", "Quality.score", 
                                "FILTER", "INFO"),
                    colClasses=c("character", "integer", "character",
                                 "character", "character", "integer", 
                                 "character", "character"))
} else {
  err("Invalid VCF file:", args$vcf)
}

# check VCF is not empty
if (dim(vcf)[1]==0) {
  mesg("No SNVs found. Removing empty VCF file.")
  file.remove(args$vcf)
  q(status=0)
}

# create output directory if necessary
x <- dirname(args$ofile)
if (!dir.exists(x)) {
  mesg("Creating annotation output directory:", x)
  dir.create(x, recursive=TRUE, mode="0740")
}
rm("x")

# format VCF file
vcf <- vcf %>%
       # remove failing SNVs
       filter(FILTER=="PASS") %>% 
       # extract info
       mutate(ID=paste0(Position, "-", REF, "-", ALT), 
              INDEL=str_detect(INFO, "INDEL"),
              Read.depth=str_extract(INFO, "(?<=DP=)[0-9]+(?=;)"),
              Frequency=str_extract(INFO, "(?<=AF=)[0-9\\.]+(?=;)"),
              Strand.bias=str_extract(INFO, "(?<=SB=)[0-9]+(?=;)")) %>%
       # format new columns
       mutate(Read.depth=as.numeric(Read.depth),
              Frequency=as.numeric(Frequency),
              Strand.bias=as.numeric(Strand.bias),
              Percent=100*Frequency)

## annotate SNVs -----------------------------------------------
options(warn=-1) # suppressing warnings
vcf <- apply(vcf, 1, annot, a=ann)
options(warn=0) # turning warnings back on
       
# join all annotations
snvs <- vcf[[1]]
if (length(vcf) > 1) {
  for(i in 2:length(vcf)) {
    x <- vcf[[i]]
    y <- intersect(colnames(snvs), colnames(x))
    snvs <- full_join(snvs, x, by=y)
  }
  rm("x", "i", "y")
}

# show table of SNVS
table(snvs$Type)

# write csv
write.csv(snvs, args$ofile, row.names=FALSE)
