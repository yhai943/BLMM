getmaf <- function(snp) {
  if (dim(data.frame(snp))[2] != 1) {
    stop(
      "getAllelelFrequencies takes a [N x 1] vector of genotypes, but ",
      "provided genotype dimensions are: [", dim(snp)[1], " x ",
      dim(snp)[2], "]"
    )
  }

  pp <- length(which(snp == 0))
  pq2 <- length(which(snp == 1))
  qq <- length(which(snp == 2))

  if (pp + pq2 + qq != length(snp[!is.na(snp)])) {
    stop("SNP vector contains alleles not encoded as 0, 1 or 2")
  }
  p <- (2 * pp + pq2) / (2 * length(snp[!is.na(snp)]))
  return(min(p, 1 - p))
}
