#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
setwd("/rprojectnb/lasvchal/Jacquelyn/papers/athletic-cluster")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))
theme_set(theme_classic())
options(stringsAsFactors=FALSE)

# helper functions
mesg <- function(...) {
  x <- paste("echo -e '[MSG]", ..., "'")
  system(x)
}
err <- function(...) {
  x <- paste("echo -e '[ERR]", ..., "'")
  system(x)
  q(status=1)
}
read.coverage <- function(id, dir=".", prefix="") {
  paste0(dir, "/", prefix, id, ".tsv") %>%
    read.delim(col.names=c("Segment", "Position", "Depth")) %>%
    mutate(Barcode=id) %>%
    select(Barcode, Position, Depth)
}
read.vcf <- function(id, dir=".", prefix="") {
  paste0(dir, "/", prefix, id, ".vcf") %>%
    read.delim(comment.char="#",
               col.names=c("Segment", "Position", "ID", "Ref", "Alt", 
                           "Qual", "Filter", "Info"),
               colClasses=c(Ref="character", Alt="character")) %>%
    filter(Filter=="PASS") %>%
    mutate(Frequency=as.numeric(str_extract(Info, "(?<=AF=)[0-9\\.]+")),
           ID=paste0(Position, "-", Ref, "-", Alt),
           Barcode=id) %>%
    select(Barcode, ID, Frequency)
}
calculate.differences <- function(m) {
  d <- matrix(nrow=dim(m)[2], 
              ncol=dim(m)[2], 
              dimnames=list(colnames(m),
                            colnames(m)))
  for(i in 1:dim(m)[2]) {
    # get first vector
    v1 <- m[,i]
    for(j in i:dim(m)[2]) {
      # get second vector
      v2 <- m[,j]
      # calculate differences
      x <- sum(v1 != v2)
      # insert in both halves of difference matrix
      d[i,j] <- x
      d[j,i] <- x
    }
  }
  return(d)
}
assign.clusters <- function(d, r=0) {
  # set upper diagonal to NA, that way we only consider "earlier" samples
  # when assigning clusters
  diag(d) <- NA
  for(i in 1:dim(d)[1]) {
    d[i, i:dim(d)[1]] <- NA
  }
  
  # init cluster vector
  cluster <- rep(NA, dim(d)[1])
  
  # assign the first sample to cluster 1
  cluster[1] <- 1
  
  # loop over rows and assign clusters
  for(i in 2:dim(d)[1]) {
    # what other sample columns are in sample #i's cluster?
    # since upper triangle is NA, we automatically don't consider "later"
    # samples or matches to oneself
    parent <- which(d[i,] %in% r)
    # if there are parent matches, find the "oldest" parent 
    #    and attach its cluster ID to sample #i
    # otherwise (no parent matches), get the max cluster ID, increment
    #    by one, and assign to sample #i
    if(length(parent) > 0) {
      parent <- min(parent)
      cluster[i] <- cluster[parent]
    } else {
      cluster[i] <- max(cluster, na.rm=TRUE) + 1
    }
  }
  return(cluster)
}

# arguments
args <- ArgumentParser()
args$add_argument("-i", "--idir", default="data", help="input directory")
args$add_argument("-o", "--odir", default="clustering", help="output directory")
args$add_argument("-m", "--meta", default="sampleinfo.csv", help="sample data")
args$add_argument("-p", "--pfix", default="bu", help="filename prefix")
args$add_argument("-t", "--thld", default=7000, help="clustering threshold")
args <- args$parse_args()

# check input directory
if(dir.exists(args$idir)) {
  mesg("Valid input directory:", args$idir)
} else {
  err("Input directory not found:", args$idir)
}

# check output directory
if(dir.exists(args$odir)) {
  mesg("Valid output directory:", args$odir)
} else {
  mesg("Creating output directory:", args$odir)
  dir.create(args$odir, recursive=TRUE)
}

# check and load metadata (barcode, collection date, lineage)
if(file.exists(args$meta)) {
  mesg("Valid metadata file:", args$meta)
  meta <- read.csv(args$meta, 
                   na.strings="",
                   colClasses=c(Collection.date="Date")) %>%
          select(Barcode, Collection.date, Lineage) %>%
          mutate(QC.sequencing=(!is.na(Lineage)))
  mesg("Metadata contains", dim(meta)[1], "samples, of which", 
       sum(meta$QC.sequencing), "pass sequencing QC")
} else {
  err("Invalid metadata file:", args$meta)
}

# print out file prefix
mesg("Using file prefix:", args$pfix)

## sequencing QC ---------------------------------------------------------------
# coverage: all TSV files in input directory
covs <- meta$Barcode %>%
        lapply(read.coverage, dir=args$idir, prefix=args$pfix) %>%
        do.call(rbind, .)

# remove 3' and 5' UTRs
covs <- filter(covs,
               Position > 265,
               Position < 29675)

# add coverage QC to metadata
meta <- covs %>%
        group_by(Barcode) %>%
        summarise(Median.depth=median(Depth),
                  Missing.bases=sum(Depth < 10),
                  .groups="drop") %>%
        mutate(Coverage.percent=100*(29903-Missing.bases)/29903,
               QC.clustering=(Missing.bases < args$thld)) %>%
        right_join(meta, by="Barcode")

# save quality stats
f <- paste0(args$odir, "/qc.csv")
mesg("Writing quality:", f)
meta %>%
  select(Barcode, QC.sequencing, QC.clustering, 
         Median.depth, Missing.bases, Coverage.percent) %>%
  write.csv(f, row.names=FALSE)
rm(f)

## pairwise comparison ---------------------------------------------------------
mesg("Calculating pairwise differences...")

# filter to only samples passing clustering QC
meta <- filter(meta, QC.clustering)

# import VCF files
snvs <- meta$Barcode %>%
        lapply(read.vcf, dir=args$idir, prefix=args$pfix) %>%
        do.call(rbind, .) %>%
        # only considering ≥50% for cluster
        filter(Frequency >= 0.5) %>%
        # add position for coverage filter
        mutate(Position=as.integer(str_extract(ID, "^[0-9]+")))

# add in coverage to remove UTRs and low-coverage bases
snvs <- covs %>%
        filter(Depth >= 10) %>%
        select(Barcode, Position) %>%
        inner_join(snvs, by=c("Barcode", "Position")) %>%
        select(Barcode, ID, Frequency)

# make SNV matrix; everything at >0.5 round up to one; reduce others to zero
matx <- snvs %>%
        reshape2::dcast(ID ~ Barcode, 
                        fill=0, 
                        value.var="Frequency") %>%
        column_to_rownames("ID") %>%
        as.matrix()
matx[matx >= 0.5] <- 1
matx[matx < 0.5] <- 0
rm(snvs)

# make difference matrix and save
difs <- calculate.differences(matx)
f <- paste0(args$odir, "/difs.csv")
mesg("Writing difference matrix:", f)
write.csv(difs, f)
rm(matx, f)

## assign clusters -------------------------------------------------------------
mesg("Assigning clusters...")
# reorder metadata and difs primarily by date and then by barcode
meta <- meta %>%
        arrange(Barcode) %>%
        arrange(Collection.date)
difs <- difs[as.character(meta$Barcode), as.character(meta$Barcode)]

# assign clusters by 0, 1, and 2 difference thresholds
meta$Cluster.dif0 <- assign.clusters(difs, r=0)
meta$Cluster.dif1 <- assign.clusters(difs, r=0:1)
meta$Cluster.dif2 <- assign.clusters(difs, r=0:2)

# write cluster output
f <- paste0(args$odir, "/clusters.csv")
mesg("Writing clusters:", f)
meta %>%
  select(Barcode, starts_with("Cluster")) %>%
  write.csv(f, row.names=FALSE)
rm(f)

## some courtesy plots ---------------------------------------------------------
mesg("Making basic plots of results...")

# difference distribution
# remember all values appear twice so divide count by 2
d <- difs
diag(d) <- NA
d %>%
  as.data.frame() %>%
  reshape::melt(id.vars=c()) %>%
  na.omit() %>% # remove diagonal
  group_by(value) %>%
  summarise(Count=n(),
            .groups="drop") %>%
  mutate(Count=Count/2) %>%
  ggplot(aes(value, Count)) +
  geom_col(col="black", fill="black") +
  labs(x="Mutation differences",
       y="Comparisons",
       title="Difference distribution")
ggsave(paste0(args$odir, "/difs.png"), 
       width=10, height=10, units="cm")
rm(d)

# clusters with n > 1 duration
meta %>%
  group_by(Cluster.dif2) %>%
  summarise(Start=min(Collection.date),
            End=max(Collection.date),
            Cases=n(),
            .groups="drop") %>%
  # remove clusters with only one case
  filter(Cases > 1) %>%
  mutate(Cluster=factor(Cluster.dif2, levels=unique(Cluster.dif2))) %>%
  ggplot() +
  geom_segment(aes(x=Start, xend=End, y=Cluster, yend=Cluster)) +
  geom_point(aes(x=Start, y=Cluster)) +
  geom_point(aes(x=End, y=Cluster)) +
  labs(x="Swab date",
       y="Cluster ID",
       title="Cluster duration (n > 1)")
ggsave(paste0(args$odir, "/clusters-duration.png"), 
       units="cm", width=12, height=12)

# clusters with n > 3 accumulation
x <- meta %>%
     group_by(Cluster.dif2) %>%
     summarise(Cases=n(),
               .groups="drop") %>%
     filter(Cases > 3) %>%
     select(Cluster.dif2) %>%
     unlist()
seq(min(meta$Collection.date),
    max(meta$Collection.date),
    by=1) %>%
  lapply(function(i) {
    meta %>%
      rename(Cluster=Cluster.dif2) %>%
      filter(Cluster %in% x,
             Collection.date <= i) %>%
      group_by(Cluster) %>%
      summarise(Cases=n(),
                .groups="drop") %>%
      right_join(expand.grid(Cluster=x, Date=i),
                 by="Cluster") %>%
      replace_na(list(Cases=0))
  }) %>%
  do.call(rbind, .) %>%
  mutate(Cluster=factor(Cluster, levels=x)) %>%
  ggplot(aes(Date, Cases, col=Cluster)) +
  geom_line(size=2) +
  labs(x="Swab date",
       y="Total cases",
       title="Cluster accumulation (n > 3)") +
  scale_color_brewer(palette="Paired")
ggsave(paste0(args$odir, "/clusters-accumulation.png"), 
       units="cm", width=12, height=12)

## all done! -------------------------------------------------------------------
mesg("Clustering complete! Outputs can be found in the", args$odir, "directory")
sessionInfo()
