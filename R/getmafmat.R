getmafmat <- function(mat) {
  nsnp <- dim(mat)[2]
  var_maf <- sapply(1:nsnp, function(i) min(getAlleleFrequencies(mat[, i])))
  return(var_maf)
}
