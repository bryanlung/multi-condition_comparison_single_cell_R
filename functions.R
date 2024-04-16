## Quality Control

getSeuratObject <- function(seurobj, 
        object.data <- Read10X(data.dir = paste("/home/bryanl/scratch/PradoData", path, "filtered_feature_bc_matrix", sep = "/"))
object <- CreateSeuratObject(counts = object.data, project = paste(object, sep ="") , min.cells = 3, min.features = 200)

getQC <- function(seurobj, object1, path, species = "mmusculus") {
        if (species == "mmusculus") {
                selection <- "^mt-"
        } else if (species == "hsapiens")
                {
                selection <- "^MT-"
        } else (species)
                {
                selection <- "Please find the mitochondrial pattern associated with 
                your species."
        }
                return(selection)
        }
        
        object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = species)
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
                print(simpson.optimal)
        b_loc <- barplot(cluster_comp/ncol(seurobj),
                 col=RColorBrewer::brewer.pal(nrow(cluster_comp), "Set3"),
                 xlab = "Clusters", ylab = "Proportion")
                 text(x=b_loc, y=colSums(cluster_comp/ncol(seurobj)),
                 labels=simpson, pos=3, cex = 0.7)
                 legend("topright",
                 c(paste("Expected Simpson =", simpson.optimal),
                 paste("Average Simpson =", round(mean(simpson), digits=2)),
                 paste("Maximum Simpson =", max(simpson))), bty="n")
        if (round(mean(simpson) > 2* simpson.optimal)) {
                print("Integration is not recommended")
        } else {
                print("Integration is recommended")
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
