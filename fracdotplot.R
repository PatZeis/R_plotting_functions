getRC_object <- function(object) {
  colsum <- Matrix::colSums(object@assays$RNA@data)
  object@assays[["RC"]] <- object@assays$RNA
  object@assays$RC@data <- t(t(as.matrix(object@assays$RNA@data))/colsum)
  return(object)
}  

fracdotplot <- function( object, genes, cluster=NULL, sampclus = NULL, samples=NULL, limup, limdown, zsc=T, map=T, leg=T, mintotal=NULL, seurat=F) {
  library(ggplot2)
  library(RColorBrewer)
  if (seurat) {
    if ( is.null(object@assays$RC@data)) { 
      stop("compute relative counts matrix for seurat object")}
    ndata <- as.matrix(object@assays$RC@data * mintotal) + 0.1
    ndata <- data.frame(ndata)
  }
  else {
    if ( length(object@cpart) == 0) {
      object@cpart <- object@cluster$kpart
    }
    if ( ! identical(names(object@cpart), colnames(object@ndata)) ) { names(object@cpart) <- colnames(object@ndata)}
    ndata <- as.matrix(object@ndata * min(object@counts)) + 0.1
    ndata <- data.frame(ndata)}

  if ( !is.null(cluster) & !is.null(samples)) {
    stop("define either clusters OR samples")
  }
  if (is.null(cluster) & is.null(samples)) {
    stop("define either clusters OR samples")
  }
  if (!is.null(cluster)) {
  genevec <- c()
  clustervec <- c()
  fraction <- c()
  scaled_mean <- c()
  log2mean <- c()
  if (seurat) {
    cpart <- object$seurat_clusters
  }
  else {
    cpart <- object@cpart
  }
  for ( i in 1:length(genes)) {
    repgene <- rep(genes[i], length(cluster))
    
    if (!is.null(sampclus)) {
      meang <- mean(as.numeric(ndata[genes[i], grep(sampclus, colnames(ndata))]))
      sdg <- sd(ndata[genes[i],grep(sampclus, colnames(ndata))])
    }
    else {
      meang <- mean(as.numeric(ndata[genes[i],]))
      sdg <- sd(ndata[genes[i],]) }
    repclus <- c()
    frac <- c()
    cent_mean <- c()
    log2_mean <- c()
    for ( n in 1:length(cluster)) {
      if (!is.null(sampclus)) {
        clus <- names(cpart[as.numeric(cpart) == cluster[n]])
        clus <- clus[grep(sampclus, clus)]
      }
      else {
        clus <- names(cpart[as.numeric(cpart) == cluster[n]])}
      leng_clus <- length(clus)
      leng_gene_in_clus <- length(which(ndata[genes[i], clus] > 0.1))
      frac <- c(frac, leng_gene_in_clus/leng_clus)
      #repclus <- c(repclus, paste("cl",cluster[n], sep="_"))
      repclus <- c(repclus, cluster[n])
      if (zsc) { cent_mean <- c(cent_mean, (mean(as.numeric(ndata[genes[i], clus])) - meang)/sdg)}
      else { log2_mean <- c(log2_mean, log2(mean(as.numeric(ndata[genes[i], clus]))))}
    }
    genevec <- c(genevec, repgene) 
    clustervec <- c(clustervec, repclus)
    fraction <- c(fraction, frac)
    if (zsc) { scaled_mean <- c(scaled_mean, cent_mean)}
    else { log2mean <- c(log2mean, log2_mean) }
  }
  if ( zsc==T ) {
    data <- data.frame(Gene = factor(genevec, levels = genes) , Cluster = factor(clustervec, levels = cluster), Fraction = fraction, Expression = scaled_mean )
  }
  else {
    data <- data.frame(Gene = factor(genevec, levels = genes) , Cluster = factor(clustervec, levels = cluster), Fraction = fraction, Expression = log2mean )  
  }
  data[which(data$Expression > limup), "Expression"] <- limup
  data[which(data$Expression < limdown), "Expression"] <- limdown
  ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7,name = "RdYlBu")))(100)
  
  frac <- ggplot(data, aes(x = Gene, y = Cluster))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank()) + theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())                                                                    
  
  if ( map==T && leg==F) {
    print(frac + geom_point(aes(size = Fraction, color = Expression))  + scale_colour_gradientn(colours = ColorRamp) + theme(legend.position="none"))
  }
  if (leg==T && map==F){
    print(frac + theme(axis.line = element_line(colour = "black")) + theme(axis.title = element_text(color = "black"), axis.text = element_text(color = "black"),axis.ticks = element_line(color = "black"),axis.text.x = element_text(angle = 90, hjust = 1)))
  }
  if (map==T && leg==T) {
    print(frac + geom_point(aes(size = Fraction, color = Expression))  + scale_colour_gradientn(colours = ColorRamp) + theme(axis.line = element_line(colour = "black")) + theme(axis.title = element_text(color = "black"), axis.text = element_text(color = "black"),axis.ticks = element_line(color = "black"),axis.text.x = element_text(angle = 90, hjust = 1)) )
  }
}

  if (!is.null(samples)) {
  genevec <- c()
  samplevec <- c()
  fraction <- c()
  fraction_log <- c()
  scaled_mean <- c()
  log2mean <- c()
  for ( i in 1:length(genes)) {
    repgene <- rep(genes[i], length(samples))
    meang <- mean(as.numeric(ndata[genes[i],]))
    sdg <- sd(ndata[genes[i],])
    repsamp <- c()
    frac <- c()
    cent_mean <- c()
    log2_mean <- c()
    for ( n in 1:length(samples)) {
      samp <- colnames(ndata)[grep(samples[n], colnames(ndata))]
      leng_samp <- length(samp)
      leng_gene_in_samp <- length(which(ndata[genes[i], samp ]> 0))
      frac <- c(frac, leng_gene_in_samp/leng_samp)
      repsamp <- c(repsamp, samples[n])
      if ( zsc ) { cent_mean <- c(cent_mean, (mean(as.numeric(ndata[genes[i], samp])) - meang)/sdg) } 
      else { log2_mean <- c(log2_mean, log2(mean(as.numeric(ndata[genes[i], samp])))) } 
    }
    genevec <- c(genevec, repgene) 
    samplevec <- c(samplevec, repsamp)
    fraction <- c(fraction, frac)
    if ( zsc ) { scaled_mean <- c(scaled_mean, cent_mean) }
    else { log2mean <- c(log2mean, log2_mean) } 
  }
  if (zsc==T) {
    data <- data.frame(Gene = factor(genevec, levels = genes) , Sample = factor(samplevec, levels = samples), Fraction = fraction, Expression = scaled_mean )
  }
  else {
    data <- data.frame(Gene = factor(genevec, levels = genes) , Sample = factor(samplevec, levels = samples), Fraction = fraction, Expression = log2mean )  
  }
  data[which(data$Expression > limup), "Expression"] <- limup
  data[which(data$Expression < limdown), "Expression"] <- limdown
  ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7,name = "RdYlBu")))(100)
  
  frac <- ggplot(data, aes(x = Gene, y = Sample))  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank()) + theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())                                                                    
  
  if ( map==T && leg==F) {
    print(frac + geom_point(aes(size = Fraction, color = Expression))  + scale_colour_gradientn(colours = ColorRamp) + theme(legend.position="none"))
  }
  if (leg==T && map==F){
    print(frac + theme(axis.line = element_line(colour = "black")) + theme(axis.title = element_text(color = "black"), axis.text = element_text(color = "black"),axis.ticks = element_line(color = "black"),axis.text.x = element_text(angle = 90, hjust = 1)))
  }
  if (map==T && leg==T) {
    print(frac + geom_point(aes(size = Fraction, color = Expression))  + scale_colour_gradientn(colours = ColorRamp) + theme(axis.line = element_line(colour = "black")) + theme(axis.title = element_text(color = "black"), axis.text = element_text(color = "black"),axis.ticks = element_line(color = "black"),axis.text.x = element_text(angle = 90, hjust = 1)) )
  }
  }
}
#fracdotplot(restinglung.integrated2, genes, cluster = cluster, sampclus = "LocksSkin", limup = 1, limdown = -1, mintotal = 2000, seurat)
