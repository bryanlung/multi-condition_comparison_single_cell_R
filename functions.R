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
                        k <- as.sparse(k)
                        DatasetName <- paste(j,i, sep = "_")
                        output_list[[DatasetName]] <- k        
                        }  
                else if (grepl(".rds$", i) | grepl(".RDS$", i) == T) {
                        j <- files$Condition[files == i] 
                        k <- readRDS(paste(i, sep = ""))
                        DatasetName <- paste(j,i, sep = "_")
                        output_list[[DatasetName]] <- k
                }
                else if (grepl(".mtx$", i) | grepl("filtered_feature_bc_matrix$", i) == T) {
                        j <- files$Condition[files == i] 
                        k <- read10X(paste(i, sep = ""))
                        k <- as.sparse(k)
                        DatasetName <- paste(j,i, sep = "_")
                        output_list[[DatasetName]] <- k
                }
                else if (grepl(".csv$", i) | grepl(".csv.gz$", i) == T) {
                        j <- files$Condition[files == i] 
                        k <- read.csv(paste(i, sep = ""), row.names = 1)
                        k <- as.sparse(k)
                        DatasetName <- paste(j,i, sep = "_")
                        output_list[[DatasetName]] <- k
                }
                else if (grepl(".tsv$", i) | grepl(".tsv.gz$", i) == T) {
                        j <- files$Condition[files == i] 
                        k <- read.csv(paste(i, sep = "\t"), row.names = 1)
                        k <- as.sparse(k)
                        DatasetName <- paste(j,i, sep = "_")
                        output_list[[DatasetName]] <- k
                } else {
                        stop("Your file type is currently not supported.")
                }
                Sys.sleep(1/100)
        }
        split_DatasetName <- strsplit(names(output_list), "_")
        Condition <- sapply(split_DatasetName, function(x){x[[1]]})
        Sample.Ident <- sapply(split_DatasetName, function(x){paste(x[c(2)],
               collapse= "_")})
        Seurat_list <- list()
        pb1 <- progress_bar$new(
                format = "  Creating Your Seurat Object [:bar] :percent in :elapsed",
                total = length(files$files), clear = FALSE, width= 60)
        for (i in seq(from = 1, to = length(files$files), by = 1)) {
                j = i
                pb1$tick()
                if (class(output_list[[i]]) != "SeuratObject") {
                        j = i         
                        k <- CreateSeuratObject(counts = output_list[i], project = 
                                Sample.Ident[j], min.cells = min.cells, 
                                min.features = min.features, meta.data = meta.data)
                        k@meta.data$Condition <- Condition[j]
                        k@meta.data$Sample.Ident <- Sample.Ident[j]
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

getQCViolinPlot <- function(seurobj) {
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
                feature2 = "percent.mt")
        plot2 <- FeatureScatter(seurobj, feature1 = "nCount_RNA", 
                feature2 = "nFeature_RNA")
        plot <- print(plot1 + plot2) 
        ViolinPlot <- print(VlnPlot(seurobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                ncol = 3))
        seurat_output <- list(plot, ViolinPlot, seurobj)
        
        print("Step 1 of quality control is completed. Please proceed to data subsetting.")
        return(seurat_output)
}
          
getSubsetThresholds <- function(seurobj, nFeature_RNA_bot, nFeature_RNA_top,
        nCount_RNA_bot, nCount_RNA_top, percent.mt.thresh, scale_factor = 10000, 
        n_features = 3000, normalization_method = c("FALSE","LogNormalize", "CLR", "RC"), 
        selection_method = c("vst", "mean.var.plot", "dispersion")) { 
                seurobj[[3]] <- subset(seurobj[[3]], subset = nFeature_RNA > nFeature_RNA_bot 
                        & nFeature_RNA < nFeature_RNA_top & percent.mt < percent.mt.thresh & 
                        nCount_RNA > nCount_RNA_bot & nCount_RNA < nCount_RNA_top)
                ViolinPlot <- print(VlnPlot(seurobj[[3]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
                seurobj[[3]] <- SCTransform(seurobj[[3]], vars.to.regress = "percent.mt", verbose = FALSE)
                if (normalization_method[2]) { 
                        if (selection_method[1]) {
                                seurobj[[3]] <- NormalizeData(seurobj[[3]], normalization.method = "LogNormalize", 
                                        scale.factor = scale_factor)
                                seurobj[[3]] <- NormalizeData(seurobj[[3]])
                                seurobj <- FindVariableFeatures(seurobj, selection.method = "vst", nfeatures = n_features)
                                print(top20 <- head(VariableFeatures(seurobj), 20))
                                plot1 <- VariableFeaturePlot(seurobj)
                                plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
                                plot <- print(plot1 + plot2)
                                seurobj <- ScaleData(seurobj)
                                seurat_output <- list(ViolinPlot, plot, seurobj[[3]])
                                return(seurat_output)
                        }
                        else if (selection_method[2]) {
                                seurobj[[3]] <- NormalizeData(seurobj[[3]], normalization.method = "LogNormalize", 
                                        scale.factor = scale_factor)
                                seurobj[[3]] <- NormalizeData(seurobj[[3]])
                                seurobj <- FindVariableFeatures(seurobj, selection.method = "mean.var.plot", nfeatures = n_features)  
                                print(top20 <- head(VariableFeatures(seurobj), 20))
                                plot1 <- VariableFeaturePlot(seurobj)
                                plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
                                plot <- print(plot1 + plot2)
                                seurobj <- ScaleData(seurobj)
                                seurat_output <- list(ViolinPlot, plot, seurobj[[3]])
                                return(seurat_output)
                       }
                        else if (selection_method[3]) { 
                                seurobj[[3]] <- NormalizeData(seurobj[[3]], normalization.method = "LogNormalize", 
                                        scale.factor = scale_factor)
                                seurobj[[3]] <- NormalizeData(seurobj[[3]])
                                seurobj <- FindVariableFeatures(seurobj, selection.method = "dispersion", nfeatures = n_features)  
                                print(top20 <- head(VariableFeatures(seurobj), 20))
                                plot1 <- VariableFeaturePlot(seurobj)
                                plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
                                plot <- print(plot1 + plot2)
                                seurobj <- ScaleData(seurobj)
                                seurat_output <- list(ViolinPlot, plot, seurobj[[3]])
                                return(seurat_output)
                        }
                }
                else if (normalization_method[3]) { 
                        if (selection_method[1]) {
                                seurobj[[3]] <- NormalizeData(seurobj[[3]], normalization.method = "CLR", 
                                        scale.factor = scale_factor)
                                seurobj[[3]] <- NormalizeData(seurobj[[3]])
                                seurobj <- FindVariableFeatures(seurobj, selection.method = "vst", nfeatures = n_features)
                                print(top20 <- head(VariableFeatures(seurobj), 20))
                                plot1 <- VariableFeaturePlot(seurobj)
                                plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
                                plot <- print(plot1 + plot2)
                                seurobj <- ScaleData(seurobj)
                                seurat_output <- list(ViolinPlot, plot, seurobj[[3]])
                                return(seurat_output)
                        }
                        else if (selection_method[2]) {
                                seurobj[[3]] <- NormalizeData(seurobj[[3]], normalization.method = "CLR", 
                                        scale.factor = scale_factor)
                                seurobj[[3]] <- NormalizeData(seurobj[[3]])
                                seurobj <- FindVariableFeatures(seurobj, selection.method = "mean.var.plot", nfeatures = n_features)
                                print(top20 <- head(VariableFeatures(seurobj), 20))
                                plot1 <- VariableFeaturePlot(seurobj)
                                plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
                                plot <- print(plot1 + plot2)
                                seurobj <- ScaleData(seurobj)
                                seurat_output <- list(ViolinPlot, plot, seurobj[[3]])
                                return(seurat_output)
                       }
                        else if (selection_method[3]) { 
                                seurobj[[3]] <- NormalizeData(seurobj[[3]], normalization.method = "CLR", 
                                        scale.factor = scale_factor)
                                seurobj[[3]] <- NormalizeData(seurobj[[3]])
                                seurobj <- FindVariableFeatures(seurobj, selection.method = "dispersion", nfeatures = n_features)
                                print(top20 <- head(VariableFeatures(seurobj), 20))
                                plot1 <- VariableFeaturePlot(seurobj)
                                plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
                                plot <- print(plot1 + plot2)
                                seurobj <- ScaleData(seurobj)
                                seurat_output <- list(ViolinPlot, plot, seurobj[[3]])
                                return(seurat_output)
                        }
                }
                else if (normalization_method[4]) { 
                        if (selection_method[1]) {
                                seurobj[[3]] <- NormalizeData(seurobj[[3]], normalization.method = "RC", 
                                        scale.factor = scale_factor)
                                seurobj[[3]] <- NormalizeData(seurobj[[3]])
                                seurobj <- FindVariableFeatures(seurobj, selection.method = "vst", nfeatures = n_features)
                                print(top20 <- head(VariableFeatures(seurobj), 20))
                                plot1 <- VariableFeaturePlot(seurobj)
                                plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
                                plot <- print(plot1 + plot2)
                                seurobj <- ScaleData(seurobj)
                                seurat_output <- list(ViolinPlot, plot, seurobj[[3]])
                                return(seurat_output)
                        }
                        else if (selection_method[2]) {
                                seurobj[[3]] <- NormalizeData(seurobj[[3]], normalization.method = "RC", 
                                        scale.factor = scale_factor)
                                seurobj[[3]] <- NormalizeData(seurobj[[3]])
                                seurobj <- FindVariableFeatures(seurobj, selection.method = "mean.var.plot", nfeatures = n_features) 
                                print(top20 <- head(VariableFeatures(seurobj), 20))
                                plot1 <- VariableFeaturePlot(seurobj)
                                plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
                                plot <- print(plot1 + plot2)
                                seurobj <- ScaleData(seurobj)
                                seurat_output <- list(ViolinPlot, plot, seurobj[[3]])
                                return(seurat_output)
                       }
                        else if (selection_method[3]) { 
                                seurobj[[3]] <- NormalizeData(seurobj[[3]], normalization.method = "RC", 
                                        scale.factor = scale_factor)
                                seurobj[[3]] <- NormalizeData(seurobj[[3]])
                                seurobj <- FindVariableFeatures(seurobj, selection.method = "dispersion", nfeatures = n_features)
                                print(top20 <- head(VariableFeatures(seurobj), 20))
                                plot1 <- VariableFeaturePlot(seurobj)
                                plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
                                plot <- print(plot1 + plot2)
                                seurobj <- ScaleData(seurobj)
                                seurat_output <- list(ViolinPlot, plot, seurobj[[3]])
                                return(seurat_output)
                        }
                }
                print("Normalization is now complete. Please proceed to clustering and doublet removal.")
}
                                
                                



## Adjusted Rand Index

getARI <- function(seurobj, sequence = 0.25, start.res = 0.25, end.res = 2) {
        for (i in seq(from = start.res, to = end.res, by = sequence)) {
                seurobj <- FindClusters(seurobj, resolution = i)
}               
        tmp <-  1
        ARI <- c()
        names1 <- c()
        for(i in seq(from = start.res, to = end.res - 0.25, by = sequence)) {
                j = i + sequence
                Var1 <- seurobj@meta.data[, paste0("RNA_snn_res.",i)]
                Var2 <- seurobj@meta.data[, paste0("RNA_snn_res.",j)]
                Var3 <- noquote(paste0("RNA_snn_res.",j))
                ARI[tmp] <- compare(Var1, Var2, method = "adjusted.rand")
                names1[tmp] <- Var3
                tmp <- tmp + 1 
        }       
        names(ARI) <- names1
        Var4 <- names(ARI[which(max(ARI)== ARI)])
        Var4 <- as.numeric(noquote(sub("RNA_snn_res.", "", Var4)))
        seurobj <- FindClusters(seurobj, resolution = Var4)
        return(seurobj)
}

## Simpson Index

getSimpson <- function(seurobj, samples, resolution) {
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

        cluster_comp <- table(seurobj@meta.data[, samples],
                seurobj@meta.data[, (paste0("RNA_snn_res.",resolution))])
        simpson <- round(simpson_index(seurobj@meta.data[, samples],
                seurobj@meta.data$seurat_clusters), digits=2)
        simpson.optimal <- round(simpson_index_optimal(
                seurobj@meta.data[, samples]), digits=2)
        Var1 <- max(colSums(cluster_comp/ncol(seurobj))) * 1.1
        b_loc <- barplot(cluster_comp/ncol(seurobj),
                 col=RColorBrewer::brewer.pal(nrow(cluster_comp), "Set3"),
                 xlab = "Clusters", ylab = "Proportion",
                 ylim = c(0, Var1))
                 text(x=b_loc, y=colSums(cluster_comp/ncol(seurobj)),
                 labels=simpson, pos=3, cex = 0.7)
                 legend("topright",
                 c(paste("Expected Simpson =", simpson.optimal),
                 paste("Average Simpson =", round(mean(simpson), digits=2)),
                 paste("Maximum Simpson =", max(simpson))), bty="n")
        if (round(mean(simpson) > 2* simpson.optimal)) {
                print("Integration is recommended")
        } else {
                print("Integration is not recommended")
        }
}

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
