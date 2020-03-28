get_cluster_composition <- function(object, cluster_size, norm=T, color = NULL, map = T, leg=T, final=F, order=NULL, seurat=F, symbols=NULL){
  if ( seurat ) {
    types <- unique(sub("\\_.+", "", names(object$seurat_clusters)))
    types <- types[order(types)]
    cluster <- as.numeric(names(table(as.numeric(object$seurat_clusters))[table(as.numeric(object$seurat_clusters)) >= cluster_size]))
    f <- as.numeric(object$seurat_clusters) %in% cluster
    cpart <- object$seurat_clusters[f]
  }
  else {
    if ( final==F) {
      object@cpart <- object@cluster$kpart
    }
    if (!is.null(symbols)) {
      colnames(object@ndata) <- symbols
    }
    if ( !identical(colnames(object@ndata), names(object@cpart))) {
      names(object@cpart) <- colnames(object@ndata)
    }  
    
    cluster <- as.numeric(names(table(object@cpart)[table(object@cpart) >= cluster_size]))
    types <- unique(sub("\\_.+", "", colnames(object@ndata)))
    types <- types[order(types)]
    
    cpart <- object@cpart[as.numeric(object@cpart) %in% cluster ]
  }
  if (is.null(color)) {
    color <- rainbow(length(types))
  }
  iels <- cbind(cluster=as.numeric(cpart),sample=sub("\\_.+", "", names(cpart)))
  rownames(iels) <- names(cpart)
  iels <- data.frame(iels)
  counts <- as.matrix(table( iels$sample,iels$cluster))
  #
  if ( norm == T){
    if (seurat) { counts <- counts/as.numeric(table(sub("\\_.+", "", colnames(object@assays$RNA))))}
    else { counts <- counts/as.numeric(table(sub("\\_.+", "", colnames(object@ndata))))}  # normalize
  }
  rel_counts <- t(t(counts)/apply(counts,2,sum))
  rel_counts <- rel_counts[,order(as.numeric(colnames(rel_counts)))]
  if ( !is.null(order)) {
    rel_counts <- rel_counts[,as.character(order)]
  }
  #rel_counts <- rel_counts[,cluster]
  
  if (map == F && leg==T) {
    #pdf(file.path(".", paste(name, ".pdf")))
    barplot(rel_counts, col = NA, border = NA, axes = FALSE, axisnames = F)
    legend( "topright", pch=20, bty="n", cex=1, legend=types, col=color) 
    #dev.off()
  }
  
  else {
    #pdf(file.path(".", paste(name, ".pdf")))
    barplot(rel_counts, main="Sample Contribution to RaceID3 cluster", ylab="% of cluster",
            xlab="", col=color, names.arg=as.character(colnames(rel_counts)), cex.names=1, las=2)
    if ( leg ==T) {
      legend( "topright", pch=20, bty="n", cex=1, legend=types, col=color) 
    }
    #dev.off()
  }
}

