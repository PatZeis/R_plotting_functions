getRC_object <- function(object) {
            colsum <- Matrix::colSums(object@assays$RNA@data)
            object@assays[["RC"]] <- object@assays$RNA
            object@assays$RC@data <- t(t(as.matrix(object@assays$RNA@data))/colsum)
            return(object)
}  


plotexpmap <- function (object, g, n = NULL, logsc = FALSE, imputed = FALSE, 
                        fr = FALSE, um = FALSE, cells = NULL, cex = 1, map = TRUE, 
                        leg = TRUE, noise = FALSE, seurat=F, mintotal=NULL, reduction="umap" ) 
{
  if (seurat) {
    if ( is.null(mintotal)){
      stop("if seurat = T, mintotal as scaling factor has to be set")
    }
    DefaultAssay(object) <- "RC"

    
    l <- as.vector(object@assays$RC@data[g, ] * mintotal + 
                     0.1)
    if (is.null(cells)) {
      cells <- colnames(object@assays$RC@data)}
    if (is.null(n)) 
      n <- g[1]
    if (logsc) {
      f <- l == 0
      l <- log2(l)
      l[f] <- NA
    }
    h <- colnames(object@assays$RC@data) %in% cells
  }
  else {
    if (length(object@tsne) == 0 & length(object@fr) == 0 & length(object@umap) == 
        0) 
      stop("run comptsne/compfr/compumap before plotlabelsmap")
    if (!is.logical(fr)) 
      stop("fr has to be TRUE or FALSE")
    if (!is.logical(um)) 
      stop("um has to be TRUE or FALSE")
    if (length(intersect(g, rownames(object@ndata))) < length(unique(g))) 
      stop("second argument does not correspond to set of rownames slot ndata of SCseq object")
    if (!is.numeric(logsc) & !is.logical(logsc)) 
      stop("argument logsc has to be logical (TRUE/FALSE)")
    if (!is.null(cells)) {
      if (sum(!cells %in% colnames(object@ndata)) > 0) 
        stop("cells has to be a subset of cell ids, i.e. column names of slot ndata")
    }
    if (fr == FALSE & um == FALSE & dim(object@tsne)[1] == 0) {
      if (dim(object@fr)[1] != 0) {
        fr <- TRUE
      }
      else if (dim(object@umap)[1] != 0) {
        um <- TRUE
      }
    }
    if (imputed & length(object@imputed) == 0) 
      stop("imputing needs to be done by running compdist with knn > 0")
    if (is.null(n)) 
      n <- g[1]
    if (is.null(cells)) 
      cells <- colnames(object@ndata)
    knn <- object@imputed$knn
    if (!noise) {
      if (length(g) == 1) {
        l <- as.vector(object@ndata[g, ] * min(object@counts) + 
                         0.1)
      }
      else {
        l <- apply(as.data.frame(as.matrix(object@ndata)[g, 
                                                         ]) * min(object@counts), 2, sum) + 0.1
      }
      if (imputed) {
        l <- apply(rbind(object@imputed$nn, object@imputed$probs), 
                   2, function(y) {
                     ind <- y[1:(knn + 1)]
                     p <- y[(knn + 2):(2 * knn + 2)]
                     sum(l[ind] * p)
                   })
      }
    }
    else {
      if (is.null(object@noise)) 
        stop("run noise analysis first!")
      if (length(g) == 1) {
        l <- as.vector(object@noise[g, ] + 0.1)
      }
      else {
        l <- apply(as.data.frame(as.matrix(object@noise)[g, 
                                                         ]), 2, sum) + 0.1
      }
    }
    if (logsc) {
      f <- l == 0
      l <- log2(l)
      l[f] <- NA
    }
    h <- colnames(object@ndata) %in% cells
  }
  mi <- min(l, na.rm = TRUE)
  ma <- max(l, na.rm = TRUE)
  ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  ColorLevels <- seq(mi, ma, length = length(ColorRamp))
  v <- round((l - mi)/(ma - mi) * 99 + 1, 0)
  if ( seurat) {
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
  pardefault <- par()
  layout(matrix(data = c(1, 3, 2, 4), nrow = 2, ncol = 2), 
         widths = c(5, 1, 5, 1), heights = c(5, 1, 1, 1))
  par(mar = c(3, 5, 2.5, 2))
  if (!leg) 
    n <- NA
  plot(c(min(d[, 1]), max(d[, 1])), c(min(d[, 2]), max(d[, 
                                                         2])), xlab = NA, ylab = NA, main = n, pch = 20, cex = 0, 
       col = "lightgrey", axes = FALSE)
  if (map) {
    v <- v[h]
    d <- d[h, ]
    kk <- order(v, decreasing = F)
    points(d[kk, 1], d[kk, 2], col = ColorRamp[v[kk]], pch = 20, 
           cex = cex)
  }
  if (leg) {
    par(mar = c(10, 2.5, 2.5, 4))
    image(1, ColorLevels, matrix(data = ColorLevels, ncol = length(ColorLevels), 
                                 nrow = 1), col = ColorRamp, xlab = "", ylab = "", 
          xaxt = "n")
    layout(1)
    par(mar = pardefault$mar)
  }
}
