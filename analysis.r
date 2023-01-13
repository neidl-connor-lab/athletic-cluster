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
# all samples
meta <- read.csv("sampleinfo.csv",
                 na.strings="",
                 colClasses=c(Collection.date="Date"))

# load difference matrix and subset
difs <- read.csv("clustering/difs.csv", row.names=1) %>%
        as.matrix()
colnames(difs) <- rownames(difs)

## figure 1 --------------------------------------------------------------------
# figure 1A: treemap of lineages
meta %>%
  na.omit() %>% 
  # count total of each lineage
  group_by(Lineage) %>%
  summarise(Cases=n(),
            .groups="drop") %>%
  # keep top 11 for coloring purposes
  arrange(desc(Cases)) %>%
  mutate(Lineage=factor(Lineage, levels=c(unique(Lineage)[1:11], "Other"))) %>%
  replace_na(list(Lineage="Other")) %>%
  # group "other" into a single line
  group_by(Lineage) %>%
  summarise(Cases=sum(Cases),
            .groups="drop") %>%
  ggplot(aes(fill=Lineage, area=Cases, label=Lineage)) +
  treemapify::geom_treemap(col="black") +
  treemapify::geom_treemap_text(col="black", place="bottomleft") +
  scale_fill_brewer(palette="Paired") +
  theme(legend.position="none")
ggsave("analysis/treemap.png", units="cm", width=10, height=10)
ggsave("analysis/treemap.tiff", units="cm", width=10, height=10)

# figure 1B: heatmap of all AY.3
m <- filter(meta, Lineage=="AY.3")
d <- difs[as.character(m$Barcode), as.character(m$Barcode)]
png("analysis/heatmap-ay3.png", 
    units="cm", res=300, width=12, height=10)
ComplexHeatmap::Heatmap(d, name="# difs", 
                        col=colfun,
                        row_dend_width=unit(1, "cm"),
                        border=TRUE,
                        show_row_names=FALSE,
                        show_column_dend=FALSE,
                        show_column_names=FALSE)
dev.off()
tiff("analysis/heatmap-ay3.tiff", 
     units="cm", res=300, width=12, height=10)
ComplexHeatmap::Heatmap(d, name="# difs", 
                        use_raster=TRUE,
                        col=colfun,
                        row_dend_width=unit(1, "cm"),
                        border=TRUE,
                        show_row_names=FALSE,
                        show_column_dend=FALSE,
                        show_column_names=FALSE)
dev.off()
rm(m, d)

## figure 2A -------------------------------------------------------------------
# add cluster info and labels to sample info
meta <- read.csv("analysis/clusterinfo.csv") %>%
        left_join(meta, by="Barcode")
# add in cluster IDs
meta <- read.csv("clustering/clusters.csv") %>%
        right_join(meta, by="Barcode")

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
ComplexHeatmap::Heatmap(d, name="# difs", 
                        use_raster=TRUE,
                        col=colfun,
                        row_dend_width=unit(1, "cm"),
                        border=TRUE,
                        show_row_names=FALSE,
                        show_column_dend=FALSE,
                        show_column_names=FALSE)
dev.off()
rm(m, d)

# figure 2B: genomics 
m <- meta %>%
     filter(Label!="Eliminated") %>%
     arrange(Source)
d <- difs[as.character(m$Barcode), as.character(m$Barcode)]
a <- columnAnnotation(Source=m$Source,
                      border=TRUE,
                      col=list(Source=colsource))
png("analysis/heatmap-genomics.png", 
    units="cm", res=300, width=17, height=8)
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
     units="cm", res=300, width=17, height=8)
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

# figure 2C: accumulation
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


## misc. stats -----------------------------------------------------------------
# nucleotide change separating A & B identical clusters
snvs.a <- read.csv("data/bu1762.csv") %>%
          filter(Position > 265,
                 Position < 29675,
                 Percent > 50) %>%
          select(NT.ID, AA.ID)
snvs.b <- read.csv("data/bu1780.csv") %>%
          filter(Position > 265,
                 Position < 29675,
                 Percent > 50) %>%
          select(NT.ID, AA.ID)
x <- which(!(snvs.b$NT.ID %in% snvs.a$NT.ID))
snvs.b[x, ]
rm(snvs.a, snvs.b, x)

# differences between eliminated samples and others (1797 & 1830)
x <- meta %>%
     filter(Label!="Eliminated") %>%
     select(Barcode) %>%
     unlist() %>%
     as.character()
min(difs["1797", x])
min(difs["1830", x])
