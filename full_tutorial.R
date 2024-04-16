library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyr)
library(knitr)
library(formatR)
library(reshape2)
library(DoubletFinder)
library(M3Drop)
library(cluster)
library(igraph)
library(edgeR)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(gprofiler2)
library(goseq)
library(CellChat)

## Preparing Files

files <- list.files(path = "your_path", recursive = F, full.name = T)
files <- as.data.frame(files)

## Name Your Conditions/Samples 

files$"Conditions/Samples" <- c("your_conditions/samples")

## Creating a Seurat Object


## Quality Control 


