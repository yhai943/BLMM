#' Compute allele frequencies.
#'
#' @param snp Alleles encoded as 0, 1 or 2 as a number
#' @description
#' Compute allele frequenciesGenotypic
#' @details
#' If AA, AB,and BB are the frequencies of the three genotypes at a locus with two alleles,
#' allele frequency is estimated as (2pp+pq) / (2N), where p is the frequency of A-allele;
#' q is the frequency of the B-allele in the sample; N is the number of samples.
#' @return pred a vector of predictions
#' @examples
#' getAlleleFrequencies(1)
#' @export

getAlleleFrequencies <- function(snp) {
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
  return(c(p, 1 - p))
}
