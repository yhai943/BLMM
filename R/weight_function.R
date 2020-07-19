weight_function <- function(data, N, type = "uw") {
  Z <- matrix(unlist(data), nrow = N, byrow = F)
  nsnp <- dim(Z)[2]
  var_maf <- sapply(1:nsnp, function(i) min(getAlleleFrequencies(Z[, i])))
  if (type == "uw") {
    w <- rep(1, length(var_maf))
  } else if (type == "beta") {
    w <- dbeta(var_maf, 1, 25)^2
  } else if (type == "wss") {
    w <- 1 / (var_maf * (1 - var_maf))
  }
  return(w)
}
