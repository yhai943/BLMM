#' Bayesian linear mixed model with multiple random effects
#'
#' @param fit \code{vb} fit from \code{vb_fit_rare_common}
#' @param epsilon a threshold for the increase of the variance
#' @param maf_beta a filtering threshold for the weight of parameter beta
#' @param maf_u a filtering threshold for the weight of parameter U
#' @usage vb_predictive(fit, epsilon = 1e-3, maf_beta = 0.5, maf_u = 0.5)
#' @examples
#' fit <- vb_fit_rarecommon(y = y_train, genotype = data, weight.type = "wss")
#' vb_predictive(fit)
#' @export
vb_predictive <- function(x, ...) {
  UseMethod("vb_predictive")
}
#' @export
vb_predictive.vb <- function(fit, epsilon = 1e-3, maf_beta = 0.5, maf_u = 0.5) {

  #------------------------------------------------------------------------------------------------------------
  # (0) initial

  y <- NULL

  B <- "Yes"

  index_test <- fit$index_test # index of the NA value in Y (index for testing samples)

  # index_fixed_marker = fit$index_fixed_marker # index of markers for fixed term

  N <- fit$N # total num of samples

  M <- fit$M # num of regions in the testing data

  Px <- fit$Px # num of SNPs in the testing data

  P <- fit$DIM # num of SNPs in the each region

  w_beta <- fit$E_w_beta # weight of beta

  w_u <- fit$E_w_u # weight of u

  geno_fix <- fit$genotype0
  # all data - fixed term

  geno_random <- lapply(1:M, function(m) attr(fit, "UD")[[m]]) # all data - random term

  geno_background <- fit$gamma_lambda # all data - background term

  X <- fit$X # data for training - fixed term

  Z <- fit$Z # data for training - random term

  m_u <- fit$m_u # mean of u

  m_beta <- fit$m_beta # mean of beta

  m_b <- fit$m_b # mean of background

  post_sigma2e <- fit$post_sigma2e # sigma2 of error

  post_sigma2b <- fit$post_sigma2b # sigma2 of error

  var_beta <- fit$post_sigma2beta # sigma2 of error

  var_u <- fit$post_sigma2u # sigma2 of error

  yc <- fit$y

  center <- fit$center

  scale <- fit$scale

  #------------------------------------------------------------------------------------------------------------
  # (1) sorting fixed effects by w_beta

  delta_sort_PVE_beta <- rep(-Inf, Px)

  # obtain the var_beta order index

  # (1) sorting random effects by w_u

  delta_sort_PVE_u <- rep(-Inf, M)

  # obtain the var_u order index

  lst_u <- sort(var_u, index.return = TRUE, decreasing = TRUE)

  #------------------------------------------------------------------------------------------------------------
  # (2) fixed effects - PVE and select top markers
  # PVE - the total proportion of variance in phenotype explained by the sparse effects and random effects terms together

  # calculate and store the PVE_beta
  PVE_beta <- sapply(1:Px, function(j) var(as.matrix(X[, j]) %*% m_beta[j]) / (post_sigma2e + var_beta + sum(var_u) + post_sigma2b))

  # finial index for fixed term # select top regions based on index_delta_PVE_beta
  fix_index <- which(w_beta >= maf_beta)

  # (2) random effects - PVE and select top regions -------------------------------------------------------------

  # calculate and store the PVE_u
  PVE_u <- sapply(1:M, function(m) var_u[[m]] / (post_sigma2e + var_beta + sum(var_u) + post_sigma2b))

  # sort the PVE_u by the pi_u (decreasing) (add 0 to be first for further calculation)

  sort_PVE_u <- append(0, sapply(1:M, function(m) sum(PVE_u[lst_u$ix][1:m])))

  # index of change
  if (length(sort_PVE_u) > 1) {
    delta_sort_PVE_u <- sapply(1:M, function(m) delta_sort_PVE_u[m] <- sort_PVE_u[m + 1] - sort_PVE_u[m])

    index_delta_PVE_u <- which(abs(delta_sort_PVE_u) > epsilon)

    len <- min(length(index_delta_PVE_u), ceiling(0.15 * M))
  } else if (length(sort_PVE_u) == 1) {
    delta_sort_PVE_u <- sort_PVE_u

    index_delta_PVE_u <- 1

    len <- 1
  }

  random_index <- lapply(lst_u, `[`, lst_u$x %in% head(lst_u$x, len))

  #------------------------------------------------------------------------------------------------------------
  # (3) calculate the PGE - the proportion of genetic variance explained by the sparse effects terms

  #------------------------------------------------------------------------------------------------------------
  # (4) fixed effects - calculate the predictive dist - fixed term mean

  # data for testing - fixed term
  # X0 <- geno_fix[index_test, index_fixed_marker]
  X0 <- geno_fix[index_test, ]

  # Predictive mean
  # all
  # beta_pred0 <- lapply(1:Px, function(j) as.matrix(X0[,j]) %*% w_beta[j] %*% m_beta[j] )

  beta_pred0 <- lapply(1:Px, function(j) as.matrix(X0[, j]) %*% m_beta[j])

  # extract the selected markers for pred
  beta_pred00 <- beta_pred0[fix_index]

  if (length(fix_index) != 0) {
    beta_pred <- Reduce(`+`, beta_pred00) + rep(0, length(index_test))
  } else {
    beta_pred <- as.vector(rep(0, length(index_test)))
  }

  # (4) random effects - calculate the predictive dist - random term mean

  # data for testing - random term
  # data for testing - random term
  Z0 <- lapply(1:M, function(m) geno_random[[m]][-index_test, ])

  Z1 <- lapply(1:M, function(m) geno_random[[m]][index_test, ])

  # Predictive mean
  # all

  # u0
  u_0 <- lapply(1:M, function(m) Z0[[m]] %*% m_u[[m]])

  # beta0 for u - prediction step 1

  u0_be <- lapply(1:M, function(m) ginv(crossprod(Z0[[m]])) %*% t(Z0[[m]]) %*% u_0[[m]])

  # u hat for u - prediction step 2

  u_pred0 <- lapply(1:M, function(m) Z1[[m]] %*% u0_be[[m]])

  # extract the selected regions for pred
  u_pred00 <- u_pred0[random_index$ix]

  if (length(random_index$ix) != 0) {
    u_pred <- Reduce(`+`, u_pred00) + rep(0, length(index_test))
  } else {
    u_pred <- as.vector(rep(0, length(index_test)))
  }

  # (4) background  - calculate the predictive dist - background term mean

  if (B == "Yes") {
    background_pred <- geno_background[index_test, ] %*% m_b
  } else {
    background_pred <- as.vector(rep(0, length(index_test)))
  }

  # Predictive variance
  #------------------------------------------------------------------------------------------------------------
  # (5) calculate the correlation - accuracy

  pred_value <- scale(beta_pred + u_pred + background_pred, T, T)

  m <- center
  s <- sqrt(scale)

  pred <- (drop(s) * pred_value) + m

  n_test <- length(index_test)

  # pred =  beta_pred + u_pred + background_pred

  return(as.numeric(pred))
}
