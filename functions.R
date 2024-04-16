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
        cluster_comp <- table(seurobj@meta.data[, samples],
                seurobj@meta.data[, (paste0("RNA_snn_res.",resolution))])
        simpson <- round(simpson_index(seurobj@meta.data[, samples],
                seurobj@meta.data$seurat_clusters), digits=2)
        b_loc <- barplot(cluster_comp/ncol(seurobj),
                col=RColorBrewer::brewer.pal(nrow(cluster_comp), "Set3"),
                text(x=b_loc, y=colSums(cluster_comp/ncol(seurobj)),
                labels=simpson, pos=3),legend("topright",
                c(paste("Avg. Simpson =", round(mean(simpson), digits=2)),
                paste("Max. Simpson =", max(simpson)) ), bty="n"))
        return(b_loc)
}
