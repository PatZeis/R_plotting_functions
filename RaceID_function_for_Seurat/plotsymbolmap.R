plotsymbolsmap <- function (object, types, subset = NULL, samples_col = NULL, cex = 0.5, 
                            fr = FALSE, um = FALSE, leg = TRUE, map = TRUE, seurat=T, reduction="umap") 
{
  if ( seurat == F)  {
    if (length(object@tsne) == 0 & length(object@fr) == 0 & length(object@umap) == 
        0) 
      stop("run comptsne/compfr/compumap before plotlabelsmap")
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
  }
  if (is.null(subset)) 
    subset <- unique(types)
  h <- sort(unique(types)) %in% subset
  if (!is.null(subset)) {
    fp <- rep(FALSE, length(types))
    fp[types %in% subset] <- TRUE
  }
  if (is.null(samples_col)) {
    samples_col <- rainbow(length(unique(types[fp])))
  }
  else {
    samples_col <- samples_col[h]
  }
  if (seurat) {
    if (reduction=="umap") {d <- object@reductions$umap@cell.embeddings}
    if (reduction=="pca") { d <- object@reductions$pca@cell.embeddings}
    if (reduction=="tsne") {d <- object@reductions$pca@cell.embeddings}
  }
  else {
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
  if (map) {
    plot(d, xlab = "", ylab = "", axes = FALSE, cex = cex, 
         pch = 20, col = "grey")
    for (i in 1:length(unique(types[fp]))) {
      f <- types == sort(unique(types[fp]))[i]
      points(d[f, 1], d[f, 2], col = samples_col[i], pch = 20, 
             cex = cex)
    }
  }
  else {
    plot(d, xlab = "", ylab = "", axes = FALSE, cex = 0, 
         pch = 20, col = "grey", xlim = c(min(d[, 1]), max(d[, 
                                                             1])), ylim = c(min(d[, 2]), max(d[, 2])))
  }
  if (leg) 
    legend("topleft", legend = sort(unique(types[fp])), col = samples_col, 
           pch = 20, cex = 0.75, bty = "n")
}