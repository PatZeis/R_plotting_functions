### enrichment scores 
require(RColorBrewer)
require(stringr)
get_enrichment <- function(object, sample, clustsize = 20, alternat = "greater") { 
  sample <- sample[order(sample)]
if (is.null(object@cpart)) {
  object@cpart <- object@cluster$kpart
}
  if ( length(sample) > 9)
    stop( "too many samples, use x <= 9")
types <- unique(sub("\\_.+", "", names(object@cpart)))
if ( sum(sample %in% types) != length(sample))
  stop("samples must be element of colnames")
pop <- length(object@cpart)
cluster <- as.numeric(names(table(object@cpart)[table(object@cpart) > clustsize]))
pvals <- list()
for ( k in 1:length(sample)) {
  size_samp <- sum(grepl(sample[k], names(object@cpart)))
  cat(paste(size_samp, "\n",sep=""))
  pvalues <- c()
  for ( i in 1:length(cluster)) {
    clusteri <- object@cpart[object@cpart == cluster[i]]
    size_cluster <- length(clusteri)
    size_samp_clust <- sum(grepl(sample[k], names(clusteri)))
    fisher <- fisher.test(matrix(c(pop - size_samp, size_samp,size_cluster - size_samp_clust, size_samp_clust),ncol=2), alternative = alternat)
    pv <- fisher$p.value
    pvalues <- append(pvalues, pv)
    names(pvalues)[i] <- paste(k, "_cl", cluster[i], sep="")
    
  }
  pvals[[k]] <- pvalues
  
}
pvals2 <- pvals
names(pvals2) <- sample

if ( length(sample) >= 3 ){
marker_col <- brewer.pal(length(sample), "Set1")}

else {
  marker_col <- brewer.pal(3, "Set1")
  marker_col <- marker_col[1:length(sample)]
}

types <- sub("(\\_|\\.).+","", colnames(object@ndata))
types_num <- as.numeric(table(types))
types <- names(table(types))

rel_size <- types_num/pop
rel_size_abs <- rel_size *100
  
values <- list()


a <- cluster

values <- list()
for ( i in 1:length(a)){
  enrichment <- c()
  x <- table(as.character(sapply(names(object@cpart[object@cpart==a[i]]), function(x) str_split(x, pattern="[_]")[[1]][1])))
  
  lbls <- names(x)
  cat(paste(lbls, "\n", sep=""))
  y <- which(types %in% lbls)  ## index for color
  slices <- as.numeric(x)/types_num[y] ### normalize with number of total cells of sample
  pct <- round(slices/sum(slices)*100)
  lbls <- paste(pct, " % ", lbls,sep=" ")
  enrichment <- pct/rel_size_abs[y]
  names(enrichment) <- lbls
  values[[i]] <- enrichment
  names(values)[i] <- paste("cluster", a[i], sep="")
  pdf(file.path(".", paste("cluster_norm_", alternat, "_", a[i],"distribution",".pdf", sep="")))
  print(pie(slices, labels=lbls,cex.main=1.5, cex.lab=0.75,col=marker_col[y], main=paste("Cluster ", a[i], " distribution ", " N=", length(object@cpart[object@cpart == a[i]]), sep=" ")))
  print(barplot(enrichment,  names.arg=c(round(-log10(pvals2[[1]][i]), 2), round(-log10(pvals2[[2]][i]),2)), col = marker_col, ylim = range(0:2)))
  legend("topright", names(x), fill=marker_col[y], cex=1, bty="n")
  dev.off()
}
}