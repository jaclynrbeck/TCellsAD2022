
plotPCACurves <- function(scRNA) {
  pca <- scRNA[["pca"]]
  eigValues <- (pca@stdev)^2
  totalVar <- sum(matrixStats::rowVars(GetAssayData(scRNA, slot="scale.data")))
  varExplained <- eigValues / totalVar
  
  plt1 <- ggplot(data.frame(y=varExplained, x=1:length(eigValues)), aes(x,y)) + 
    geom_line() + labs(title = "Variance explained")
  
  plt2 <- ggplot(data.frame(y=cumsum(varExplained), x=1:length(eigValues)), aes(x,y)) + 
    geom_line() + labs(title = "Variance explained (cumulative)")
  
  thresh <- length(which(cumsum(varExplained) < 0.2)) + 1
  thresh2 <- length(which(cumsum(eigValues/sum(eigValues)) < 0.8)) + 1
  print(paste0("20% variance explained at x = ", thresh))
  print(paste0("If using ", length(eigValues), 
               " dimensions only, 80% variance explained at x = ", thresh2))
  plt1 / plt2
}
