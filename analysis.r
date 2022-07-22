#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
setwd("/rprojectnb/lasvchal/Jacquelyn/papers/athletic-cluster")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ComplexHeatmap))
options(stringsAsFactors=FALSE)
theme_set(theme_classic())

# helper variables
colfun <- circlize::colorRamp2(0:3, c("#de2d26", "#000000", "#bdbdbd", "#ffffff"))
colsource <- c("Genomics"="#377eb8", "Contact tracing"="#4daf4a")

## inputs ----------------------------------------------------------------------
# cluster info and labels
meta <- read.csv("analysis/clusterinfo.csv")
# add in dates
meta <- read.csv("sampleinfo.csv",
                 colClasses=c(Collection.date="Date")) %>%
        select(Barcode, Collection.date) %>%
        right_join(meta, by="Barcode")
# add in cluster IDs
meta <- read.csv("clustering/clusters.csv") %>%
        right_join(meta, by="Barcode")

# load difference matrix and subset
difs <- read.csv("clustering/difs.csv", row.names=1) %>%
        as.matrix()
colnames(difs) <- rownames(difs)

## figure 1 --------------------------------------------------------------------
# figure 1A: heatmap of all BU
png("analysis/heatmap-all.png", 
    units="cm", res=300, width=12, height=10)
ComplexHeatmap::Heatmap(difs, name="# difs", 
                        col=colfun,
                        row_dend_width=unit(1, "cm"),
                        border=TRUE,
                        show_row_names=FALSE,
                        show_column_dend=FALSE,
                        show_column_names=FALSE)
dev.off()
tiff("analysis/heatmap-all.tiff", 
     units="cm", res=300, width=12, height=10)
ComplexHeatmap::Heatmap(difs, name="# difs", 
                        use_raster=TRUE,
                        col=colfun,
                        row_dend_width=unit(1, "cm"),
                        border=TRUE,
                        show_row_names=FALSE,
                        show_column_dend=FALSE,
                        show_column_names=FALSE)
dev.off()

## figure 2A -------------------------------------------------------------------
# heatmap of epi only
m <- filter(meta, Source=="Contact tracing")
d <- difs[as.character(m$Barcode), as.character(m$Barcode)]
png("analysis/heatmap-epi.png", 
    units="cm", res=300, width=10, height=8)
ComplexHeatmap::Heatmap(d, name="# difs", 
                        col=colfun,
                        row_dend_width=unit(1, "cm"),
                        border=TRUE,
                        show_row_names=FALSE,
                        show_column_dend=FALSE,
                        show_column_names=FALSE)
dev.off()
tiff("analysis/heatmap-epi.tiff", 
     units="cm", res=300, width=10, height=8)
ComplexHeatmap::Heatmap(difs, name="# difs", 
                        use_raster=TRUE,
                        col=colfun,
                        row_dend_width=unit(1, "cm"),
                        border=TRUE,
                        show_row_names=FALSE,
                        show_column_dend=FALSE,
                        show_column_names=FALSE)
dev.off()
rm(m, d)

## figure 3 --------------------------------------------------------------------
# figure 3A: genomics 
# heatmap of epi only
m <- meta %>%
     filter(Label!="Eliminated") %>%
     arrange(Source)
d <- difs[as.character(m$Barcode), as.character(m$Barcode)]
a <- columnAnnotation(Source=m$Source,
                      border=TRUE,
                      col=list(Source=colsource))
png("analysis/heatmap-genomics.png", 
    units="cm", res=300, width=12, height=8)
ComplexHeatmap::Heatmap(d, name="# difs", 
                        col=colfun,
                        row_dend_width=unit(1, "cm"),
                        top_annotation=a,
                        border=TRUE,
                        show_row_names=FALSE,
                        show_column_dend=FALSE,
                        show_column_names=FALSE)
dev.off()
tiff("analysis/heatmap-genomics.tiff", 
     units="cm", res=300, width=12, height=8)
ComplexHeatmap::Heatmap(d, name="# difs", 
                        col=colfun,
                        row_dend_width=unit(1, "cm"),
                        top_annotation=a,
                        border=TRUE,
                        show_row_names=FALSE,
                        show_column_dend=FALSE,
                        show_column_names=FALSE)
dev.off()
rm(m, d, a)

# figure 3B: accumulation
meta %>%
  filter(Label!="Eliminated") %>%
  arrange(Collection.date, Source) %>%
  rownames_to_column("Cases") %>%
  mutate(Days=Collection.date-min(Collection.date),
         Cases=as.integer(Cases)) %>%
  ggplot(aes(Days, Cases, fill=Source)) +
  geom_point(pch=21, size=3) +
  scale_fill_manual(values=colsource) +
  labs(x="Days since first case",
       y="Case count",
       title="Accumulation of cases over time") +
  theme(legend.position="bottom")
ggsave("analysis/accumulation.png", units="cm", width=10, height=10)
ggsave("analysis/accumulation.tiff", units="cm", width=10, height=10)
