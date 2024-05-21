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
library(progress)
library(Matrix)
library(glmGamPoi)

## Preparing Files

files <- list.files(path = "your_path", recursive = F, full.name = T)
files <- list.files(path = "/home/bryanl/scratch/MBI4850G/finalproject",recursive = F, full.name = T)
files <- list.files(path = "/home/bryanl/scratch/testdatasets",recursive = F, full.name = T)

files <- as.data.frame(files)
files$files <- files

## Name Your Conditions/Samples 

files$Condition <- c("your_conditions_names")
files$Samples <- c("your_sample_names")

## Creating a Seurat Object

YourVariable <- getAllSeuratObject(files)
YourVariable <- unlist(YourVariable)

## Quality Control 


