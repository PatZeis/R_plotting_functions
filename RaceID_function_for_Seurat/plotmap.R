library(cluster)
compmedoids <- function (objectassay, part) 
  
{
  m <- c()
  for (i in sort(unique(as.numeric(part)))) {
    f <- names(part)[as.numeric(part) == i]
    if (length(f) == 1) {
      m <- append(m, f)
    }
    else {
      
      
      g <- apply(as.matrix(objectassay@data[, as.numeric(part) == 
                                                           i]) - as.vector(pam(t(objectassay@data[, as.numeric(part) == 
                                                                                                                 i]), 1)$medoids), 2, sum) == 0
      m <- append(m, names(part)[as.numeric(part) == i][g])
      
    }
  }
  m
}


plotmap <- function (object, final = TRUE, tp = 1, fr = FALSE, um = FALSE, 
                     cex = 0.5, seurat=F, medoids, reduction="umap") 
{ 
  
  if (seurat) {
    if (reduction=="umap") {d <- object@reductions$umap@cell.embeddings}
    if (reduction=="pca") { d <- object@reductions$pca@cell.embeddings}
    if (reduction=="tsne") {d <- object@reductions$tsne@cell.embeddings}
    part <- object$seurat_clusters
  }
  else {
    if (length(object@tsne) == 0 & length(object@fr) == 0 & length(object@umap) == 
        0) 
      stop("run comptsne/compfr/compumap before plotlabelsmap")
    if (final & length(object@cpart) == 0) 
      stop("run findoutliers before plotmap")
    if (!final & length(object@cluster$kpart) == 0) 
      stop("run clustexp before plotmap")
    if (!is.numeric(tp) | (is.numeric(tp) & tp > 1 | tp < 0)) 
      stop("tp has to be a number between 0 and 1 (transparency)")
    if (!is.logical(fr)) 
      stop("fr has to be TRUE or FALSE")
    if (!is.logical(um)) 
      stop("um has to be TRUE or FALSE")
    if (fr == FALSE & um == FALSE & dim(object@tsne)[1] == 0) {
      if (dim(object@fr)[1] != 0) {
        fr <- TRUE
      }
      else if (dim(object@umap)[1] != 0) {
        um <- TRUE
      }
    }
    part <- if (final) 
      object@cpart
    else object@cluster$kpart
    if (fr) {
      d <- object@fr
    }
    else if (um) {
      d <- object@umap
    }
    else {
      d <- object@tsne
    }
  }
  if (is.null(object@fcol)) { fcol <- rainbow(length(unique(as.numeric(part))))}
  else { fcol <- object@fcol}
  row.names(d) <- names(part)
  plot(d, xlab = "", ylab = "", cex = 0, axes = FALSE)
  for (i in 1:max(as.numeric(part))) {
    if (sum(as.numeric(part) == i) > 0) 
      points(d[as.numeric(part) == i, 1], d[as.numeric(part) == i, 2], col = adjustcolor(fcol[i], 
                                                                                         tp), pch = 20, cex = cex)
  }
  for (i in 1:max(as.numeric(part))) {
    if (sum(as.numeric(part) == i) > 0) 
      points(d[medoids[i], 1], d[medoids[i], 
                                 2], col = adjustcolor(fcol[i], tp), pch = 20, 
             cex = 4)
    if (sum(as.numeric(part) == i) > 0) 
      points(d[medoids[i], 1], d[medoids[i], 
                                 2], col = adjustcolor("white", tp), pch = 20, 
             cex = 3)
    if (sum(as.numeric(part) == i) > 0) 
      text(d[medoids[i], 1], d[medoids[i], 
                               2], i, col = adjustcolor("black", tp), cex = 0.75, 
           font = 4)
  }
}

## example
## plotmap(restinglung.integrated2, seurat = T, medoids = medoids)
