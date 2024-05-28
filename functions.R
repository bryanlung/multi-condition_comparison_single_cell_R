## Data Manipulation

getAllSeuratObject <- function(files, min.cells = 3, min.features = 200,
                              meta.data = NULL) {
        output_list <- list()
        pb <- progress_bar$new(
                format = "  Importing Your Data [:bar] :percent in :elapsed",
                total = length(files$files), clear = FALSE, width= 60)
        for (i in files$files) {
                pb$tick()
                if (grepl("txt.gz$", i) | grepl("txt.gz$", i) == T) {
                        j <- files$Condition[files == i] 
                        k <- read.delim(paste(i, sep = ""), row.names = 1)
                        k <-  as.sparse(k)
                        l <- files$Samples[files == i]
                        DatasetName <- paste(j,l,i, sep = "_")
                        output_list[[DatasetName]] <- k        
                        }  
                else if (grepl(".rds$", i) | grepl(".RDS$", i) == T) {
                        j <- files$Condition[files == i] 
                        k <- readRDS(paste(i, sep = ""))
                        l <- files$Samples[files == i]
                        DatasetName <- paste(j,l,i, sep = "_")
                        output_list[[DatasetName]] <- k
                }
                else if (grepl(".mtx$", i) | grepl("filtered_feature_bc_matrix$", i) == T) {
                        j <- files$Condition[files == i] 
                        k <- read10X(paste(i, sep = ""))
                        k <- as.sparse(k)
                        l <- files$Samples[files == i]
                        DatasetName <- paste(j,l,i, sep = "_")
                        output_list[[DatasetName]] <- k
                }
                else if (grepl(".csv$", i) | grepl(".csv.gz$", i) == T) {
                        j <- files$Condition[files == i] 
                        k <- read.csv(paste(i, sep = ""), row.names = 1)
                        k <- as.sparse(k)
                        l <- files$Samples[files == i]
                        DatasetName <- paste(j,l,i, sep = "_")
                        output_list[[DatasetName]] <- k
                }
                else if (grepl(".tsv$", i) | grepl(".tsv.gz$", i) == T) {
                        j <- files$Condition[files == i] 
                        k <- read.csv(paste(i, sep = "\t"), row.names = 1)
                        k <- as.sparse(k)
                        l <- files$Samples[files == i]
                        DatasetName <- paste(j,l,i, sep = "_")
                        output_list[[DatasetName]] <- k
                } else {
                        stop("Your file type is currently not supported.")
                }
                Sys.sleep(1/100)
        }
        split_DatasetName <- strsplit(names(output_list), "_")
        Condition <- sapply(split_DatasetName, function(x){x[[1]]})
        Samples <- sapply(split_DatasetName, function(x){x[[2]]})
        Sample.Ident <- sapply(split_DatasetName, function(x){paste(x[c(3)],
               collapse= "_")})
        Seurat_list <- list()
        pb1 <- progress_bar$new(
                format = "  Creating Your Seurat Object [:bar] :percent in :elapsed",
                total = length(files$files), clear = FALSE, width= 60)
        for (i in seq(from = 1, to = length(files$files), by = 1)) {
                j <- i
                pb1$tick()
                if (class(output_list[[i]]) != "SeuratObject") {
                        j <- i         
                        k <- CreateSeuratObject(counts = output_list[i], 
                                project = Sample.Ident[j], min.cells = min.cells, 
                                min.features = min.features, meta.data = meta.data)
                        k@meta.data$Condition <- Condition[j]
                        k@meta.data$Sample.Ident <- Sample.Ident[j]
                        k@meta.data$Samples <- Samples[j]
                        SeuratObjectName <- paste(names(output_list)[j], "SeuratObject", sep = "_")
                        Seurat_list[[SeuratObjectName]] <- k
                }
                Sys.sleep(1/100)
        }
        data_list <- list(output_list = output_list, Seurat_list = Seurat_list)
        print(paste("Importing data is now completed.", length(files$files), 
                    "Seurat objects were created."))
        return(data_list)
}

SeuratMerge <- function(SeurObj) {
        Merged_list <- list()
        N = length(SeurObj$Seurat_list)
        if (N == 2) { 
                GlobalMerged <- merge(SeurObj$Seurat_list[[1]], SeurObj$Seurat_list[[2]])
                CompleteMerge <- paste("CompleteMerge")
                Merged_list[[CompleteMerge]] <- GlobalMerged
                print("Dataset merging is now completed.")
                return(Merged_list)
        } else if (N > 2) {
                GlobalMerged <- merge(SeurObj$Seurat_list[[1]], SeurObj$Seurat_list[[2]])
                for(i in seq(from = 3, to = length(SeurObj$Seurat_list), by = 1)) {
                        GlobalMerged <- merge(GlobalMerged, SeurObj$Seurat_list[[i]])
                        CompleteMerge <- paste("CompleteMerge")
                        Merged_list[[CompleteMerge]] <- GlobalMerged
                }
        print("Dataset merging is now completed.")
        return(Merged_list)        
        }
}
      
recSeuratMerge <- function(SeurObj) {
        N <- length(SeurObj) + 0.1
        if (N == 1.1) {
                GlobalMerged <- SeurObj[[1]]
                return(GlobalMerged)
        }
        if (N == 2.1) {
                GlobalMerged <- merge(SeurObj[[1]], SeurObj[[2]])
                return(GlobalMerged) 
        } else { 
                set1 = 1:floor(N/2)
                set2 = ceiling(N/2):N
                a = recSeuratMerge(SeurObj[set1]) 
                b = recSeuratMerge(SeurObj[set2])
                GlobalMerged <- merge(a,b)
                return(GlobalMerged)
        }
}

getRecursiveMerge <- function(SeurObj) {
        A <- recSeuratMerge(SeurObj)
        print(paste("Dataset merging is now completed.", length(files$files), 
                "Seurat objects were merged."))
        return(A)
}

## Quality Control

getQCViolinPlot <- function(seurobj, SavePlots = c("FALSE", "TRUE")) {
        if (length(grep( "^mt-", rownames(seurobj), value = T)) > 1) {
                print("Mouse dataset detected.")
                seurobj[["percent.mt"]] <- PercentageFeatureSet(seurobj, pattern = "^mt")
        }
        else if (length(grep( "^mt.", rownames(seurobj), value = T)) > 1) { 
                print("Mouse dataset detected.")
                seurobj[["percent.mt"]] <- PercentageFeatureSet(seurobj, pattern = "^mt.")
        }
        else if (length(grep( "^MT-", rownames(seurobj), value = T)) > 1) { 
                print("Human dataset detected.")
                seurobj[["percent.mt"]] <- PercentageFeatureSet(seurobj, pattern = "^MT-")
        } else {
                stop("Mitochondrial pattern not found. Please find the mitochondrial pattern associated with 
                your species.")
        }
        plot1 <- FeatureScatter(seurobj, feature1 = "nCount_RNA", 
                feature2 = "percent.mt", group.by = "Condition")
        plot2 <- FeatureScatter(seurobj, feature1 = "nCount_RNA", 
                feature2 = "nFeature_RNA", group.by = "Condition")
        plot <- print(plot1 + plot2) 
        ViolinPlot <- print(VlnPlot(seurobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                ncol = 3, group.by = "Condition"))
        if (SavePlots[1] == "TRUE") { 
                seurobj@misc$QCplots <- plot
                seurobj@misc$QCViolinPlot <- ViolinPlot
                print("Step 1 of quality control is completed. Please proceed to data subsetting.")
                return(seurobj)
        
        }
        if (SavePlots[1] == "FALSE") { 
                print("Step 1 of quality control is completed. Please proceed to data subsetting.")
                return(seurobj)
        }
}
          
getSubsetThresholds <- function(seurobj, nFeature_RNA_bot, nFeature_RNA_top,
        nCount_RNA_bot, nCount_RNA_top, percent.mt.thresh, scale_factor = 10000, 
        n_features = 3000, normalization_method = c("FALSE","LogNormalize", "CLR"), 
        selection_method = c("FALSE","vst", "mean.var.plot", "dispersion"), 
        SavePlots = c("FALSE", "TRUE")) { 
                seurobj <- subset(seurobj, subset = nFeature_RNA > nFeature_RNA_bot 
                        & nFeature_RNA < nFeature_RNA_top & percent.mt < percent.mt.thresh & 
                        nCount_RNA > nCount_RNA_bot & nCount_RNA < nCount_RNA_top)
                ViolinPlot <- print(VlnPlot(seurobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,
                        group.by = "Condition"))
                if (normalization_method[1] == "LogNormalize") { 
                        seurobj <- NormalizeData(seurobj, normalization.method = "LogNormalize", 
                                scale.factor = scale_factor)
                        if (selection_method[1] == "vst") {
                                seurobj <- FindVariableFeatures(seurobj, selection.method = "vst", nfeatures = n_features)
                                print(top20 <- head(VariableFeatures(seurobj), 20))
                                plot1 <- VariableFeaturePlot(seurobj)
                                plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
                                plot <- print(plot1 + plot2)
                                seurobj <- ScaleData(seurobj)
                                if (SavePlots[1] == "TRUE") {
                                        seurobj@misc$FinalViolinPlot <- ViolinPlot
                                        seurobj@misc$HVGPlot <- plot
                                        print("Normalization is now complete. Please proceed to clustering.")
                                        return(seurobj)
                                }
                                if (SavePlots[1] == "FALSE") { 
                                        print("Normalization is now complete. Please proceed to clustering.")
                                        return(seurobj)
                                }
                        }
                        else if (selection_method[1] == "mean.var.plot") {
                                seurobj <- FindVariableFeatures(seurobj, selection.method = "mean.var.plot", nfeatures = n_features)  
                                print(top20 <- head(VariableFeatures(seurobj), 20))
                                plot1 <- VariableFeaturePlot(seurobj)
                                plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
                                plot <- print(plot1 + plot2)
                                seurobj <- ScaleData(seurobj)
                                if (SavePlots[1] == "TRUE") {
                                        seurobj@misc$FinalViolinPlot <- ViolinPlot
                                        seurobj@misc$HVGPlot <- plot
                                        print("Normalization is now complete. Please proceed to clustering.")
                                        return(seurobj)
                                }
                                if (SavePlots[1] == "FALSE") { 
                                        print("Normalization is now complete. Please proceed to clustering.")
                                        return(seurobj)
                                }
                       }
                       else if (selection_method[1] == "dispersion") { 
                                seurobj <- FindVariableFeatures(seurobj, selection.method = "dispersion", nfeatures = n_features)  
                                print(top20 <- head(VariableFeatures(seurobj), 20))
                                plot1 <- VariableFeaturePlot(seurobj)
                                plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
                                plot <- print(plot1 + plot2)
                                seurobj <- ScaleData(seurobj)
                                if (SavePlots[1] == "TRUE") {
                                        seurobj@misc$FinalViolinPlot <- ViolinPlot
                                        seurobj@misc$HVGPlot <- plot
                                        print("Normalization is now complete. Please proceed to clustering.")
                                        return(seurobj)
                                }
                                if (SavePlots[1] == "FALSE") { 
                                        print("Normalization is now complete. Please proceed to clustering.")
                                        return(seurobj)
                                }
                        }
                }
                else if (normalization_method[1] == "CLR") { 
                        seurobj <- NormalizeData(seurobj, normalization.method = "CLR", 
                                scale.factor = scale_factor)
                        if (selection_method[1] == "vst") {
                                seurobj <- FindVariableFeatures(seurobj, selection.method = "vst", nfeatures = n_features)
                                print(top20 <- head(VariableFeatures(seurobj), 20))
                                plot1 <- VariableFeaturePlot(seurobj)
                                plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
                                plot <- print(plot1 + plot2)
                                seurobj <- ScaleData(seurobj)
                                if (SavePlots[1] == "TRUE") {
                                        seurobj@misc$FinalViolinPlot <- ViolinPlot
                                        seurobj@misc$HVGPlot <- plot
                                        print("Normalization is now complete. Please proceed to clustering.")
                                        return(seurobj)
                                }
                                if (SavePlots[1] == "FALSE") { 
                                        print("Normalization is now complete. Please proceed to clustering.")
                                        return(seurobj)
                                }
                        }
                        else if (selection_method[1] == "mean.var.plot") {
                                seurobj <- FindVariableFeatures(seurobj, selection.method = "mean.var.plot", nfeatures = n_features)
                                print(top20 <- head(VariableFeatures(seurobj), 20))
                                plot1 <- VariableFeaturePlot(seurobj)
                                plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
                                plot <- print(plot1 + plot2)
                                seurobj <- ScaleData(seurobj)
                                if (SavePlots[1] == "TRUE") {
                                        seurobj@misc$FinalViolinPlot <- ViolinPlot
                                        seurobj@misc$HVGPlot <- plot
                                        print("Normalization is now complete. Please proceed to clustering.")
                                        return(seurobj)
                                }
                                if (SavePlots[1] == "FALSE") { 
                                        print("Normalization is now complete. Please proceed to clustering.")
                                        return(seurobj)
                                }
                       }
                        else if (selection_method[1] == "dispersion") { 
                                seurobj <- FindVariableFeatures(seurobj, selection.method = "dispersion", nfeatures = n_features)
                                print(top20 <- head(VariableFeatures(seurobj), 20))
                                plot1 <- VariableFeaturePlot(seurobj)
                                plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
                                plot <- print(plot1 + plot2)
                                seurobj <- ScaleData(seurobj)
                                if (SavePlots[1] == "TRUE") {
                                        seurobj@misc$FinalViolinPlot <- ViolinPlot
                                        seurobj@misc$HVGPlot <- plot
                                        print("Normalization is now complete. Please proceed to clustering.")
                                        return(seurobj)
                                }
                                if (SavePlots[1] == "FALSE") { 
                                        print("Normalization is now complete. Please proceed to clustering.")
                                        return(seurobj)
                                }
                        }
                }
                else if (normalization_method[1] == "FALSE") {
                        if (selection_method[1] == "FALSE") {
                                seurobj <- SCTransform(seurobj, vars.to.regress = "percent.mt", verbose = FALSE)
                                print(top20 <- head(VariableFeatures(seurobj), 20))
                                Var1 <- seurobj@assays$SCT@SCTModel.list[[1]]@feature.attributes
                                Var1$Genes <- rownames(Var1)
                                plot <- print(ggplot(Var1, aes(gmean, residual_variance, label = rownames(Var1))) + 
                                        geom_point() +
                                        scale_y_log10() + scale_x_log10() + geom_text(aes(label = Genes), 
                                        data= Var1[Var1$Genes %in% top20,], hjust=0, vjust=0, color = "red") +
                                        theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),
                                        panel.background = element_blank(),
                                        axis.line = element_line(colour = "black")))
                                if (SavePlots[1] == "TRUE") {
                                        seurobj@misc$FinalViolinPlot <- ViolinPlot
                                        seurobj@misc$HVGPlot <- plot
                                        print("Normalization is now complete. Please proceed to clustering.")
                                        return(seurobj)
                                }
                                if (SavePlots[1] == "FALSE") { 
                                        print("Normalization is now complete. Please proceed to clustering.")
                                        return(seurobj)
                                }
                        }
                }
}
                    
##subtest <- getSubsetThresholds(seurobj,200,2000,200,5000,2, normalization_method= "FALSE", selection_method= "FALSE", SavePlots = "TRUE")                             

## Clustering and Doublet Removal  

getPCs <- function(seurobj, SavePlots = c("FALSE","TRUE")) {
        all.genes <- rownames(seurobj)
        seurobj <- ScaleData(seurobj, features = all.genes)
        seurobj <- ScaleData(seurobj, vars.to.regress = "percent.mt")
        seurobj <- RunPCA(seurobj, features = VariableFeatures(object = seurobj))
        Var1 <- print(DimPlot(seurobj, reduction = "pca", group.by = "Condition"))
        Elbow <- print(ElbowPlot(seurobj, ndims = 50))
        data <- Elbow$data
        x <- data$dims
        y <- data$stdev
        reg <- lm(tail(y, 30) ~ tail(x, 30))
        PCtable <- matrix(nrow = length(x), ncol = 1)
        j <- 1
        for (i in y) {
                tmp <- (i - (reg$coefficients[2]*j + reg$coefficients[1]))
                PCtable[j,] <- tmp
                j <- j + 1
        }
        advisedPCs <- which(PCtable < 0)[1] - 1
        Var2 <- print(ggplot(data, aes(dims, stdev)) + geom_point() + geom_smooth(data = tail(data, 30), 
                method = "lm" ,color = "red",fullrange=TRUE) +
                theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black")))
        seurobj@misc$advisedPCs <- advisedPCs
        seurobj@misc$PClinearregression <- reg
        if (SavePlots[1] == "TRUE") {
                seurobj@misc$PCAPlot <- Var1
                seurobj@misc$OrigElbowPlot <- Elbow
                seurobj@misc$LinearRegressionFittedElbow <- Var2
                print(paste("Principal component analysis complete.", advisedPCs, "principal components are suggested for downstream analysis.")) 
                return(seurobj)
        }
        if (SavePlots[1] == "FALSE") {        
                print(paste("Principal component analysis complete.", advisedPCs, "principal components are suggested for downstream analysis."))
                return(seurat_output)
        }
}
        

getClusters <- function(seurobj, dim = advisedPCs, sequence = 0.25, 
        start.res = 0.25, end.res = 2, SCT = c("TRUE", "FALSE")) {   
        seurobj <- FindNeighbors(seurobj, dims = 1:dim)
        for (i in seq(from = start.res, to = end.res, by = sequence)) {
                seurobj <- FindClusters(seurobj, resolution = i)
        }               
        tmp <-  1
        ARI <- c()
        names1 <- c()
        nSample <- round(5000/length(unique(seurobj@meta.data$seurat_clusters)))
        TotalSampledCells <- list()
        pb <- progress_bar$new(
                format = "  Calculating silhouette score [:bar] :percent in :elapsed",
                total = length(unique(seurobj@meta.data$seurat_clusters)), 
                        clear = FALSE, width= 60)
        for (i in 0:(length(unique(seurobj@meta.data$seurat_clusters))-1)) {
                pb$tick()
                getCol <- colnames(seurobj[, seurobj@meta.data$seurat_clusters == i])
                SampledCells <- sample(getCol, size = min(nSample, length(getCol)), replace = F)
                TotalSampledCells[[tmp]] <- SampledCells
                tmp <- tmp + 1
                Sys.sleep(1/100)
        }
        tmp <- 1
        downsampled.obj <- seurobj[, unlist(TotalSampledCells)]
        avg_cluster_sil_scores <- c()
        for(i in seq(from = start.res, to = end.res - sequence, by = sequence)) {
                j = i + sequence
                if (SCT[1] == "TRUE") { 
                        Var1 <- seurobj@meta.data[, paste0("SCT_snn_res.",i)]
                        Var2 <- seurobj@meta.data[, paste0("SCT_snn_res.",j)]
                        Var3 <- noquote(paste0("SCT_snn_res.",j))
                }
                if (SCT[1] == "FALSE") { 
                        Var1 <- seurobj@meta.data[, paste0("RNA_snn_res.",i)]
                        Var2 <- seurobj@meta.data[, paste0("RNA_snn_res.",j)]
                        Var3 <- noquote(paste0("RNA_snn_res.",j))
                }                        
                ARI[tmp] <- compare(Var1, Var2, method = "adjusted.rand")
                names1[tmp] <- Var3
                if (SCT[1] == "TRUE") {
                        Idents(downsampled.obj) <- downsampled.obj@meta.data[, paste0("SCT_snn_res.",j)]
                }
                if (SCT[1] == "FALSE") { 
                        Idents(downsampled.obj) <- downsampled.obj@meta.data[, paste0("RNA_snn_res.",j)]
                }
                downsampled.obj@meta.data$seurat_clusters <- Idents(downsampled.obj)
                dist.matrix <- dist(x = Embeddings(object = downsampled.obj[["pca"]])[, 1:dim])
                sil <- silhouette(as.numeric(as.character(downsampled.obj@meta.data$seurat_clusters)), dist=dist.matrix)
                cluster_sil_scores <- aggregate(sil[,3], by=list(sil[,1]), mean)    
                avg_cluster_sil_scores[tmp] <- mean(cluster_sil_scores$x)
                tmp <- tmp + 1 
        }       
        names(ARI) <- names1
        names(avg_cluster_sil_scores) <- names1
        print("ARI")
        print(ARI)
        print("Silhouette Score")
        print(avg_cluster_sil_scores)
        Var4 <- ARI + avg_cluster_sil_scores
        Var5 <- names(Var4[which(max(Var4) == Var4)])
        Idents(seurobj) <- seurobj@meta.data[, Var5]
        seurobj@meta.data$seurat_clusters <- Idents(seurobj) 
        seurobj@misc$suggested_res <- Var5
        seurobj@misc$ARIs <- ARI
        seurobj@misc$Silhouette_Scores <- avg_cluster_sil_scores
        if (SCT[1] == "TRUE") {
                Var5 <- as.numeric(noquote(sub("SCT_snn_res.", "", Var5)))
        }
        if (SCT[1] == "FALSE") {
                Var5 <- as.numeric(noquote(sub("RNA_snn_res.", "", Var5)))
        }
        print(paste("A resolution of", Var5, "was used for clustering. Please proceed to doublet removal."))
        return(seurobj)
}

## Doublet Removal at UMAP Generation

findDoublets <- function(seurobj, dim = advisedPCs, SCT = c("FALSE", "TRUE")) { 
        doublet_list <- list()
        for (i in unique(seurobj@meta.data$Samples)) { 
                type.freq <- table(seurobj@meta.data$seurat_clusters[seurobj@meta.data$Samples == i])/ncol(seurobj[,seurobj@meta.data$Samples == i])
                homotypic.prop <- sum(type.freq^2)
                nEXP = 0.009*(ncol(seurobj[,seurobj@meta.data$Samples == i])/1000)*(1-homotypic.prop)*ncol(seurobj[,seurobj@meta.data$Samples == i])
                pN <- 0.25
                PC <- dim
                if (SCT[1] == "TRUE") {
                        sweep.out <- paramSweep(seurobj[,seurobj@meta.data$Samples == i], PCs=1:PC, sct = T)
                }
                if (SCT[1] == "FALSE") {
                        sweep.out <- paramSweep(seurobj[,seurobj@meta.data$Samples == i], PCs=1:PC, sct = F)
                }
                sweep.stats <- summarizeSweep(sweep.out)
                maxBCreal <- data.frame(which(sweep.stats == max(sweep.stats$BCreal), arr.ind=TRUE))
                DoubletParameters <- print(sweep.stats[maxBCreal$row,])
                pK <- as.numeric(as.character(DoubletParameters$pK))
                pN <- as.numeric(as.character(DoubletParameters$pN))
                if (SCT[1] == "TRUE") {
                        tmpseurobj <- doubletFinder(seurobj[,seurobj@meta.data$Samples == i], PCs=1:PC, pN=pN, pK=pK, nExp=nEXP, sct = T)
                        X <- paste("DF.classifications", pN, pK, nEXP, sep="_")
                        Var1 <- tmpseurobj@meta.data[,X]
                        doublet_list[[i]] <- Var1
                }
                if (SCT[1] == "FALSE") {
                        tmpseurobj <- doubletFinder(seurobj[,seurobj@meta.data$Samples == i], PCs=1:PC, pN=pN, pK=pK, nExp=nEXP, sct = F) 
                        X <- paste("DF.classifications", pN, pK, nEXP, sep="_")
                        Var1 <- tmpseurobj@meta.data[,X]
                        doublet_list[[i]] <- Var1
                }
        }
        seurobj@meta.data$DFCLASSIFICATIONS <- unlist(doublet_list)
        print("Doublet removal is now completed. Please proceed to UMAP generation and cell annotation.")
        return(seurobj)
}

getMarkers <- function(seurobj, dim = advisedPCs, SavePlots = c("FALSE", "TRUE"), only_pos = TRUE, 
        min_pct = 0.01, logfc_threshold = -Inf, SCT = c("FALSE", "TRUE")) {
                seurobj <- RunUMAP(seurobj, dims=1:dim)
                UMAP <- print(DimPlot(seurobj, reduction = "umap"))
                UMAP_Condition <- print(DimPlot(seurobj, reduction = "umap", group.by = "Condition"))
                UMAP_DFCLASSFICATIONS <- print(DimPlot(seurobj, reduction="umap", group.by= "DFCLASSIFICATIONS"))
                if (SCT[1] == "TRUE") {
                         Var1 <- PrepSCTFindMarkers(seurobj)
                         Var1 <- FindAllMarkers(Var1, only.pos = only_pos, min.pct = min_pct, logfc.threshold = logfc_threshold,
                                 assay = "SCT")
                }
                if (SCT[1] == "FALSE") {
                        Var1 <- FindAllMarkers(seurobj, only.pos = only_pos, min.pct = min_pct, logfc.threshold = logfc_threshold)
                }
                Var1 %>%
                group_by(cluster) %>%
                dplyr::filter(avg_log2FC > 1)
                seurobj@misc$annotation_markers <- Var1
                if (SavePlots[1] == "TRUE") {
                        seurobj@misc$UMAP <- UMAP
                        seurobj@misc$UMAP_Condition <- UMAP_Condition
                        seurobj@misc$UMAP_DFCLASSFICATIONS  <- UMAP_DFCLASSFICATIONS 
                        print("UMAP is now generated. Please use marker genes to annotate clusters.")
                        return(seurobj)
                }
                if (SavePlots[1] == "FALSE") { 
                        print("UMAP is now generated. Please use marker genes to annotate clusters.")
                        return(seurobj)
                }
}

## Simpson Index with respect to each condition

getSimpson <- function(seurobj, SCT = c("TRUE", "FALSE")) {
        simpson_index <- function(sample, cluster) {
                tab <- table(sample, cluster)
                tab <- t(t(tab)/colSums(tab))
                return(colSums(tab^2))
        }

        simpson_index_optimal <- function(samples) {
        tab <- table(samples)
        tab <- tab/sum(tab)
        return(sum(tab^2))
        }
        tmp <- 1
        totalsimpson <- c()
        totalsimpson.optimal <- c()
        names1 <- c()
        output_list <- list()
        for (i in unique(seurobj@meta.data$Condition)) { 
                tmpseurobj <- seurobj[, seurobj@meta.data$Condition == i]
                simpson <- round(simpson_index(tmpseurobj@meta.data[, "Samples"],
                        tmpseurobj@meta.data$seurat_clusters), digits=2)
                simpson.optimal <- round(simpson_index_optimal(
                tmpseurobj@meta.data[, "Samples"]), digits=2)
                simpson[is.na(simpson)] = 0
                simpson.optimal[is.na(simpson.optimal)] = 0 
                totalsimpson[tmp] <- round(mean(simpson), digits=2)
                totalsimpson.optimal[tmp] <- simpson.optimal
                names1[tmp] <- i
                names(totalsimpson) <- names1
                names(totalsimpson.optimal) <- names1
                output_list[[i]] <- simpson
                tmp <- tmp + 1
        }
        simpson <- as.data.frame(output_list)
        simpson$cluster <- rownames(simpson)
        simpson <- melt(simpson)
        Var1 <- length(unique(simpson$cluster)) - 1
        simpson$clusters <- factor(simpson$cluster, levels = c(0, 1:Var1))
        totalsimpson.optimal <- as.data.frame(totalsimpson.optimal)
        totalsimpson.optimal$Condition <- rownames(totalsimpson.optimal)
        totalsimpson.optimal$Placeholder <- rep("A", length(rownames(totalsimpson.optimal)))
        p <- ggplot(simpson, aes(x = clusters , y = variable , fill = value)) + geom_tile() +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "white")) +
                labs(x= "Clusters" , y = NULL) +  scale_fill_gradientn(colours = c('#ffffcc','#a1dab4','#41b6c4','#2c7fb8','#253494'),
                values = c(0, 0.25, 0.5, 0.75, 1), name= "Simpson") + theme(plot.title = element_text(hjust = 0.5)) + 
                ggtitle("Simpson per Cluster")
        q <- ggplot(totalsimpson.optimal, aes(x = Placeholder ,y = Condition, fill = totalsimpson.optimal)) + 
                geom_tile() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "white"), 
                axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
                scale_fill_gradientn(colours = c('#ffffcc','#a1dab4','#41b6c4','#2c7fb8','#253494'),
                limits = c(0, 1), values = c(0, 0.25, 0.5, 0.75, 1), name= "Simpson") + labs(x = NULL, y = NULL) + 
                ggtitle("Expected Simpson") + theme(plot.title = element_text(hjust = 0.5)) + 
                geom_text(aes(label= totalsimpson.optimal))
        SimpsonPlot <- print(ggarrange(p, q, common.legend = TRUE, legend="bottom", widths = c(1,0.5),  align = "h"))
        Var2 <- print(round(mean(totalsimpson), digits=2))
        Var3 <- print(round(mean(totalsimpson.optimal$totalsimpson.optimal), digits=2))
        seurobj@misc$simpson <- Var2
        seurobj@misc$simpson.optimal <- Var3
        seurobj@misc$SimpsonPlot <- SimpsonPlot 
        if (Var1 > 2* Var2) {
                print("Integration is recommended")
        } else {
                print("Integration is not recommended")
        }
        return(seurobj)
}   

## Automatic Annotator

autoAnnotate <- function(seurobj) { 
        mmusculus.database <- MouseRNAseqData()
        tmpseurobj <- seurobj
        tmpseurobj[["RNA"]] <- as(tmpseurobj[["RNA"]], Class="Assay")
        tmpseurobj <- as.SingleCellExperiment(tmpseurobj)
        singleR_annot <- SingleR(test = tmpseurobj, ref = mmusculus.database, assay.type.test=1, 
                labels = mmusculus.database$label.main)
        seurobj@meta.data$singleR <- singleR_annot$pruned.labels

        
## Pseudobulk Expression Matrix

group_rowmeans <- function(MAT, group_labs, type=c("mean","sum")) {
        d <- split(seq(ncol(MAT)), factor(group_labs));
        if (type[1] == "mean") {
                mus <- sapply(d, function(group) my_rowMeans(MAT[,group]))
        } else {
                mus <- sapply(d, function(group) my_rowSums(MAT[,group]))
        }
        return(mus);
}

my_rowSums <- function(x) {
        if (!is.null(ncol(x))) {
                if (ncol(x) > 1) {
                        return(Matrix::rowSums(x))
                }
        }
        return(x);
}

get_pseudobulk <- function(mat, clusters, donors, nexclude=10) {
          tab <- table(clusters, donors)
          exclude = tab < nexclude & tab > 0
          exclude <- which(exclude, arr.ind=T)

          tmp <- paste(clusters, donors, sep="_")
          exclude_tmp <- paste(rownames(tab)[exclude[,1]], colnames(tab)[exclude[,2]], sep="_")
          remove <- tmp %in% exclude_tmp
          mat <- mat[,!remove]
          clusters <- factor(clusters[!remove])
          donors <- factor(donors[!remove])
          c <- split(seq(ncol(mat)), clusters);
          donor_freqs <- table(donors)/length(donors)
          # avg expression per donor in this cluster
          clust_expr <- sapply(c, function(clust) {
                d_expr <- group_rowmeans(mat[,clust], donors[clust], type="sum");
                if(is.null(dim(d_expr))) {
                        l <- sapply(d_expr, length)
                        keep <- which(l == nrow(mat))
                        d_expr <- matrix(d_expr[[keep]], ncol=length(keep), byrow=FALSE);
                        rownames(d_expr) <- rownames(mat);
                        colnames(d_expr) <- paste(clusters[clust[1]], levels(donors)[keep], sep="_")
                } else {
                        colnames(d_expr) <- paste(clusters[clust[1]], colnames(d_expr), sep="_")
                }
                return(d_expr);
          })
          out <- clust_expr[[1]];
          for (i in 2:length(clust_expr)) {
                c_names <- c(colnames(out), colnames(clust_expr[[i]]))
                out <- cbind(out, clust_expr[[i]]);
                if (is.null(dim(out))){
                        out <- matrix(out, ncol=1)
                        rownames(out) <- rownames(mat)
                }
                colnames(out) <- c_names
         }
         return(out)
}
