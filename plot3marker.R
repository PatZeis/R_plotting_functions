plot3marker <- function(object, logsc=T, cluster=NULL, all=F, gene1, gene2, gene3, leg=T, map=T, cor = F, cex=6, exclude=NULL){
  names(object@cpart) <- colnames(object@ndata)
  ma <- as.matrix(object@ndata* min(object@counts)) + 0.1
  ma <- data.frame(ma)
  if ( all == T && !is.null(cluster)) {
    stop("either specific cluster or all clusters")
  }
  if (!is.null(cluster)) {
    if ( length(cluster) > 1) { 
    cl <- names(object@cpart)[object@cpart %in% cluster]
    main <- paste("cluster", paste(cluster, collapse = "_"), sep = "_")}
    else {
    cl <- names(object@cpart)[object@cpart == cluster]
    main = paste("cluster", cluster, sep = "_")}}
  if ( all == T) {
    cl <- names(object@cpart)
    main <- "all clusters"
  }
  if (!is.null(exclude)) {
    cl <- cl[!grepl(exclude,cl)]
  }
  if ( cor ){
    corel <- cor(as.numeric(ma[gene1,cl]), as.numeric(ma[gene2,cl]))}
  else { corel <- ""}
  l <- as.numeric(ma[gene3, cl])
  m1 <- cbind(as.numeric(ma[gene1,cl]), as.numeric(ma[gene2,cl]))
  m1 <- data.frame(m1)
  if (logsc==T) {
    f <- l == 0
    l <- log2(l)
    l[f] <- NA
  }
  mi <- min(l, na.rm = TRUE)
  ma <- max(l, na.rm = TRUE)
  ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  ColorLevels <- seq(mi, ma, length = length(ColorRamp))
  v <- round((l - mi)/(ma - mi) * 99 + 1, 0)
  kk <- order(v, decreasing = F)
  pardefault <- par()
  layout(matrix(data = c(1, 3, 2, 4), nrow = 2, ncol = 2), 
         widths = c(5, 1, 5, 1), heights = c(5, 1, 1, 1))
  par(mar = c(4, 5, 2.5, 2))
  if (!leg) {
    n <- NA
    plot(c(min(m1[, 1]), max(m1[, 1])), c(min(m1[, 2]), max(m1[, 
                                                               2])), xlab = NA, ylab = NA, main = paste("expression in ", main, " corelation=", corel, sep=""), pch = 20, cex = 0, 
         col = "lightgrey", axes = FALSE)}
  else {
    plot(c(min(m1[, 1]), max(m1[, 1])), c(min(m1[, 2]), max(m1[, 
                                                               2])), xlab = gene1, ylab = gene2, main = paste("expression in ", main, " corelation=", corel, sep="") , pch = 20, cex = 0, 
         col = "lightgrey" )}
  
  if (map) {
    
    points(m1[kk,1], m1[kk,2], col = ColorRamp[v[kk]],pch =20, cex=cex)
    #points(m1[kk,1], m1[kk,2], col = "black",pch =1, cex=4)
  }
  if (leg) {
    par(mar = c(10, 2.5, 2.5, 4))
    image(1, ColorLevels, matrix(data = ColorLevels, ncol = length(ColorLevels), 
                                 nrow = 1), col = ColorRamp, xlab = gene3, ylab = "", 
          xaxt = "n")
    layout(1)
    par(mar = pardefault$mar)
  }
  
  
}