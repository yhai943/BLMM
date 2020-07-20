#' Fit Variational Bayes model
#'
#' @param y a phenotype vector of length n
#' @param genotype0 a list of genotype matrices of the target sequence
#' @param max_iter maximum number of iteration
#' @param weight.type type of weight function
#' @param maf.filter.snp a filtering threshold of minor allele frequency for the isolated predictors
#' @param epsilon_conv a numeric value giving the interval endpoint
#' @param Bsigma2beta a numeric value for sigma beta
#' @param theta_beta probability of causal variants
#' @param theta_u probability of causal regions
#' @param verbose informative iteration messages
#' @param kernel kernel type for covariance matrix
#' @description
#' Take in genotype, phenotype and analyses the target sequence by using BLMM.
#' @details
#' A hybrid model that includes a sparsity regression model and a LMM with multiple random effects.
#' The sparsity regression model part is designed to capture the strong predictive effects from isolated
#' predictors, whereas the LMM part is to capture the effects from a group of predictors located in
#' nearby genetic regions.
#' @examples
#' data("data", package = "BLMM")
#' # choose model type: 1. "uw" for BLMM-UW; 2. "beta" for BLMM-BETA; 3. "wss" for BLMM-WSS;
#' fit <- vb_fit_rarecommon(y = y_train, genotype = data, weight.type = "wss")
#' @importFrom MASS ginv
#' @importFrom data.table first
#' @importFrom boot logit inv.logit
#' @importFrom BeSS bess.one
#' @importFrom matrixcalc hadamard.prod
#' @importFrom psych tr
#'
#' @export

vb_fit_rarecommon <- function(y, genotype0, max_iter = 25000, weight.type = NULL, maf.filter.snp = 0.01, epsilon_conv = 1e-4,
                              Bsigma2beta = 1, theta_beta = 0.1, theta_u = 0.1, verbose = TRUE, kernel = "Lin") {

  #################################################
  #                                               #
  #               0.  Initial Input               #
  #                                               #
  #################################################
  if (verbose) {
    cat("[1] Initializing the algorithm:", "\n")
  }
  telliter <- 100

  #--------------------- A. import data --------------------------
  # options(digits=40)

  # Number of observations
  N <- NROW(y)

  # index of the NA value in Y (index of testing data)
  index_test <- as.vector(which(is.na(y)))


  #-------------------------- B. screening for random term and background --------------------------
  M0 <- length(genotype0)

  genotype0 <- lapply(1:M0, function(m) genotype0[[m]][, colSums(is.na(genotype0[[m]])) == 0])

  # y for training without NA
  yc <- scale(y, center = T, scale = T)[-index_test]

  center <- attr(scale(y, center = T, scale = T), "scaled:center")
  scale <- attr(scale(y, center = T, scale = T), "scaled:scale")
  #--------------------------- select the gene start------------------------------

  genotype <- genotype0

  #--------------------------- select the gene end------------------------------
  M <- length(genotype)

  regions_mat <- lapply(1:M, function(m) as.matrix(scale((genotype[[m]]), center = T, scale = F)))


  # combine all snps from all genes
  tot_rare_common_var0 <- do.call("cbind", genotype)

  # center the big matrix with all snps
  tot_rare_common_var <- as.matrix(scale(tot_rare_common_var0, center = T, scale = F))
  #-------------------------- B. screening for fixed term --------------------------
  index_var_snp <- lapply(1:M0, function(m) which(getmafmat(genotype0[[m]]) > maf.filter.snp))

  # obtain the common snps by using index

  var_snp <- lapply(1:M0, function(m) genotype0[[m]][, index_var_snp[[m]]])

  # combine all snps from all genes
  tot_snp_only0 <- do.call("cbind", var_snp)


  # screening ---------------------------------------
  Np_total <- dim(tot_snp_only0)[2]
  Np <- min(N, Np_total, 10)
  # Best subset selection
  screening_fix <- bess.one(tot_snp_only0[-index_test, ], yc, s = Np, family = "gaussian", normalize = TRUE)

  # obtain the top markers
  subset_fix <- screening_fix$bestmodel$model$xbest
  # obtain the index of top markers
  index_fixed_marker <- attr(subset_fix, "dimnames")[[2]]
  index_fixed_marker2 <- as.numeric(gsub("[^0-9.]", "", index_fixed_marker))

  tot_snp_only <- tot_snp_only0[, index_fixed_marker2]

  X <- tot_snp_only[-index_test, ]

  Px <- round(NCOL(X)) # P num of SNPs in each region

  ## -------------------------- add gene_add to genotype ----------------------------------

  ## ------------------------------------------------
  ## get the weight for each genetic region

  if (is.null(weight.type) == FALSE) {
    w0 <- lapply(1:M, function(m) weight_function(genotype[[m]], N = N, type = weight.type))

    if (weight.type == "uw") {
      w <- lapply(1:M, function(m) diag(as.numeric(scale(w0[[m]], F, scale = F))))
      pve <- 0.8
    } else if (weight.type == "beta") {
      w <- lapply(1:M, function(m) diag(as.numeric(scale(w0[[m]], F, scale = F))))
      pve <- 0.8
    } else if (weight.type == "wss") {
      w <- lapply(1:M, function(m) diag(as.numeric(scale(w0[[m]], F, scale = F))))
      pve <- 0.8
    }

    # data=genotype[[2]]; N = N; type = weight.type;
  }

  gc(verbose = FALSE)

  #-------------------------- C. reconstruct the covarience matrix for random term--------------------------

  # 1.  Z original SNP matrix

  NP <- lapply(1:M, function(m) dim(regions_mat[[m]])[2])
  ZList0 <- lapply(1:M, function(m) scale(regions_mat[[m]], F, F))

  # 3. select the kernel function - Lin, RBF, IBS, POLY2
  ## weight function - uw, wss, beta
  ## calculate kernel matrix

  if (kernel == "Lin" && is.null(weight.type) == TRUE) {
    Cov_List <- lapply(1:M, function(m) tcrossprod(ZList0[[m]])) # linear effects for SNPs
  } else if (kernel == "RBF" && is.null(weight.type) == TRUE) {
    Cov_List <- lapply(1:M, function(m) (kernelMatrix(rbfdot(sigma = 0.3), ZList0[[m]])@.Data)) # RBF effects for SNPs
  } else if (kernel == "POLY2" && is.null(weight.type) == TRUE) {
    Cov_List <- lapply(1:M, function(m) (tcrossprod(ZList0[[m]]) + 1)**2) # Poly 2 effects for SNPs
  } else if (kernel == "Lin" && is.null(weight.type) == FALSE) {
    Cov_List <- lapply(1:M, function(m) (ZList0[[m]] %*% w[[m]] %*% t(ZList0[[m]])) / (NP[[m]]))
    # Cov_List <- lapply(1:M, function(m) (ZList0[[m]] %*% w[[m]] %*% t(ZList0[[m]])) )
  } else if (kernel == "quadratic" && is.null(weight.type) == FALSE) {
    Cov_List <- lapply(1:M, function(m) ((ZList0[[m]] %*% w[[m]] %*% t(ZList0[[m]]) + 1)**2) / (NP[[m]]))
    #   Cov_List <- lapply(1:M, function(m) ((ZList0[[m]] %*% w[[m]] %*%t(ZList0[[m]]) + 1)**2) )
  } else {
    stop("currently only support Lin, RBF, POLY2 and IBS kernel and uw, beta, wss weight.type")
  }
  # a. Proportion of Variance
  #  temp <- lapply(1:M, function(m) princomp(ZList00[[m]],cor = TRUE, scores=TRUE) )
  temp <- lapply(1:M, function(m) prcomp(Cov_List[[m]]))
  zv <- lapply(1:M, function(m) temp[[m]]$sdev)
  pv <- lapply(1:M, function(m) zv[[m]]^2)
  cs <- lapply(1:M, function(m) pv[[m]] / sum(pv[[m]]))
  # b. Cumulative Proportion
  cm <- lapply(1:M, function(m) cumsum(cs[[m]]))
  # c. select number of pc explaining 80% of variance

  # 4. reconstruct the covariance matrix using eigenvalues and eigen vectors - random effects

  EigenList <- lapply(1:M, function(m) eigen(Cov_List[[m]])) # get eigen values and vectors of the covariance

  EigenvalueList0 <- lapply(1:M, function(m) ifelse(EigenList[[m]]$values < 0, 0, EigenList[[m]]$values)) # set negative eigen values to zero

  DIM <- lapply(1:M, function(m) max(min(length(which(EigenvalueList0[[m]] > 1e-30)), 45), 1))

  indexlist <- lapply(1:M, function(m) sort(EigenvalueList0[[m]], index.return = TRUE, decreasing = TRUE)) # get the eigen values sorted index

  vectorlist <- lapply(1:M, function(m) EigenList[[m]]$vectors[, indexlist[[m]]$ix[1:DIM[[m]]]]) # extract the eigen vectors with top eigen values

  EigenvalueList <- lapply(1:M, function(m) indexlist[[m]]$x[1:DIM[[m]]]) # get the top eigen values index (Top P)
  UListD <- lapply(1:M, function(m) (vectorlist[[m]] %*% diag(sqrt(EigenvalueList[[m]][1:DIM[[m]]]), DIM[[m]])) / 1000) # reconstruct X - (beta = U * D^(1/2))

  # reconstructed Z matrix

  Z <- lapply(1:M, function(m) UListD[[m]][-index_test, ])

  tUDUList <- lapply(1:M, function(m) crossprod(UListD[[m]][-index_test, ])) # t(U)*D*U

  # reconstructed covariance

  tZZList <- tUDUList

  #-----------------------------------------------------------------------------------------------------------------
  # 5. reconstruct the covariance matrix using eigenvalues and eigen vectors - background

  NP0 <- dim(tot_rare_common_var)[2]
  tot_var <- tot_rare_common_var0 # dropNA all variants

  if (is.null(weight.type) == FALSE) {
    weight.type <- "uw"

    w_g0 <- weight_function(tot_rare_common_var0, N = N, type = weight.type)

    if (weight.type == "uw") {
      w_g <- diag(as.numeric(scale(w_g0, F, scale = F)))
    } else {
      w_g <- diag(as.numeric(scale(w_g0, F, scale = T)))
    }
  }

  Cov_b <- (tot_var %*% w_g %*% t(tot_var)) / NP0 # all snps

  Eigen_b <- eigen(Cov_b) # get eigen values and vectors of the covariance

  Eigenvalue_b0 <- ifelse(Eigen_b$values < 0, 0, Eigen_b$values) # set negative eigen values to zero

  index_b <- sort(Eigenvalue_b0, index.return = TRUE, decreasing = TRUE) # get the eigen values sorted index

  # DIM0 <- length(which(Eigenvalue_b0!=0)) # remove eigenvalue 0

  # a. Proportion of Variance for g0 ----------------------------
  temp0 <- prcomp(Cov_b)
  zv0 <- temp0$sdev
  pv0 <- zv0^2
  cs0 <- pv0 / sum(pv0)
  # b. Cumulative Proportion
  cm0 <- cumsum(cs0)
  # c. select number of pc explaining 80% of variance
  DIM0 <- first(which(cm0 > 0.9))

  vector_b <- Eigen_b$vectors[, index_b$ix[1:DIM0]] # extract the eigen vectors with top eigen values

  Eigenvalue_b <- index_b$x[1:DIM0] # get the top eigen values index (Top P)

  # Gamma * (Lambda)^1/2

  gamma_lambda <- vector_b %*% diag(sqrt(Eigenvalue_b[1:DIM0])) # reconstruct X - (beta = U * D^(1/2))

  # reconstructed g matrix

  gamma_lambda_b <- gamma_lambda[-index_test, ] # gamma*(lambda)^(1/2)

  tglg <- crossprod(gamma_lambda_b) # t(gamma)*lambda*gamma

  # reconstructed covariance - background

  # tZZList <- tUDUList

  gc(verbose = FALSE)

  #--------------------------- D. initial the parameters for iterations ----------------------------------

  # s_beta

  # E_theta <- lapply(1:M, function(m) A_theta_u0[[m]]/(A_theta_u0[[m]] + B_theta_u0[[m]]) )

  # fixed theta and sigma of beta for fixed term
  #### fix ------------------------
  A_sigma2beta0 <- 0.1
  B_sigma2beta0 <- 0.1

  A_sigma2beta <- 0.1
  B_sigma2beta <- 0.1

  A_theta_beta0 <- 0.1
  B_theta_beta0 <- 0.1

  E_theta_beta <- rep(0.1, Px)

  B_sigma2beta <- 0.1 # fix

  # m_beta

  m_beta <- rep(1 / (Px * N), Px)

  E_w_beta <- rep(0.1, Px)

  XGammaBeta <- X %*% diag(E_w_beta) %*% m_beta

  logit_gamma_beta <- rep((0.05 / 0.95), Px)

  E_w_beta <- sapply(1:Px, function(j) inv.logit(logit_gamma_beta[j]))

  omega <- tcrossprod(E_w_beta) + hadamard.prod(diag(E_w_beta), (diag(1, Px) - diag(E_w_beta)))

  # B_sigma2beta <- rep(0.1, Px)

  #### random effects ------------------------

  A_theta_u0 <- 0.1
  B_theta_u0 <- 0.1

  A_sigma2u0 <- rep(0.1, M)
  B_sigma2u0 <- rep(0.1, M)

  A_sigma2u <- rep(A_sigma2u0 + N / 2, M)
  B_sigma2u <- rep(0.1, M)

  m_u <- lapply(1:M, function(m) rep(0.1, DIM[[m]]))

  logit_gamma_u <- rep(0.05 / 0.95, M)

  E_w_u <- lapply(1:M, function(m) inv.logit(logit_gamma_u[[m]]))

  E_theta_u <- rep(0.1, M)

  ZGammaU0 <- lapply(1:M, function(m) Z[[m]] %*% (drop(E_w_u[[m]]) * diag(1, DIM[[m]])) %*% m_u[[m]])

  ZGammaU <- do.call("cbind", ZGammaU0)

  # E_theta_u <- rep(E_theta_u00,M) ## random

  ####  background ------------------------
  m_b <- rep(0.1, DIM0)

  A_sigma2b <- 0.1
  B_sigma2b <- 0.1

  A_sigma2b0 <- 0.1
  B_sigma2b0 <- 0.1

  g <- gamma_lambda_b %*% m_b

  ####  error ------------------------

  A_sigma2e0 <- 0.1
  B_sigma2e0 <- 0.1

  A_sigma2e <- A_sigma2e0 + N / 2
  B_sigma2e <- 0.1

  #----------------------------------------------------------------------
  #----------------------------------------------------------------------
  #----------------------------------------------------------------------

  # Store the lower bounds L
  L <- rep(-Inf, max_iter)
  if (verbose) {
    cat("[1] Done!", "\n")
  }
  if (verbose) {
    cat("[2] Starting variational inference:", "\n")
  }

  # Iterate to find optimal parameters
  for (i in 2:max_iter) {

    #################################################
    #                                               #
    #      1.  Initial Paramters and Iteration      #
    #                                               #
    #################################################

    #----------------------# random effects #----------------------#

    B <- sapply(1:M, function(m) yc - XGammaBeta - rowSums(ZGammaU[, -m, drop = F]) - g)

    # Um

    s_u <- lapply(1:M, function(m) {
      solve(drop(((A_sigma2e0 + N / 2) / B_sigma2e)) * hadamard.prod(crossprod(Z[[m]]), (drop(E_w_u[[m]]) * diag(1, DIM[[m]]))) +
        drop(((A_sigma2u0[[m]] + DIM[[m]] / 2)) / B_sigma2u[[m]]) * diag(1, DIM[[m]]))
    })

    m_u <- lapply(1:M, function(m) (drop((A_sigma2e0 + N / 2) / B_sigma2e) * s_u[[m]]) %*% ((drop(E_w_u[[m]])) * diag(1, DIM[[m]])) %*% t(Z[[m]]) %*% B[, m, drop = F])

    # Gamma w_u

    logit_gamma_u <- lapply(1:M, function(m) {
      logit(E_theta_u[m]) + drop((A_sigma2e0 + N / 2) / B_sigma2e) * (
        t(B[, m, drop = F]) %*% Z[[m]] %*% m_u[[m]] - (1 / 2) * (
          t(m_u[[m]]) %*% t(Z[[m]]) %*% Z[[m]] %*% m_u[[m]] +
            tr(crossprod(Z[[m]]) %*% s_u[[m]])
        )
      )
    })

    E_w_u <- lapply(1:M, function(m) inv.logit(logit_gamma_u[[m]]))

    ZGammaU <- sapply(1:M, function(m) Z[[m]] %*% (drop(E_w_u[[m]]) * diag(1, DIM[[m]])) %*% m_u[[m]])

    # sigma B of U_m

    B_sigma2u <- lapply(1:M, function(m) {
      B_sigma2u0[[m]] + (1 / 2) * (crossprod(m_u[[m]]) + tr(s_u[[m]])
      )
    })

    #  B_sigma2u <- lapply(1:M, function(m) 0.1 )

    # theta for random effects

    A_theta_u <- sapply(1:M, function(m) E_w_u[[m]] + A_theta_u0)

    B_theta_u <- sapply(1:M, function(m) 1 - E_w_u[[m]] + B_theta_u0)

    E_theta_u <- sapply(1:M, function(m) A_theta_u[m] / (A_theta_u[m] + B_theta_u[m]))
    #  E_theta_u  <- rep(theta_u,M)

    #----------------------# fixed effects #----------------------#

    A <- yc - rowSums(ZGammaU[, , drop = F]) - g

    # Beta

    s_beta <- solve(drop(((A_sigma2e0 + N / 2) / B_sigma2e)) * hadamard.prod(crossprod(X), omega) +
      drop((A_sigma2beta0 + Px / 2) / B_sigma2beta) * diag(1, Px))

    m_beta <- drop(((A_sigma2e0 + N / 2) / B_sigma2e)) * s_beta %*% diag(E_w_beta) %*% t(X) %*% A


    # sigma2 for fixed effects

    A_sigma2beta <- (1 / 2) * N + A_sigma2beta0

    mtm_beta <- crossprod(m_beta)

    tr_s_beta <- tr(s_beta)

    B_sigma2beta <- (1 / 2) * (mtm_beta + tr_s_beta) + B_sigma2beta0

    # B_sigma2beta <- B_sigma2beta0

    # Gamma w_beta

    # use X_E_w_beta[,-j] = X[,-j] %*% diag(E_w_beta[-j]) to replace (for computation purpose)
    X_E_w_beta <- X %*% diag(E_w_beta)

    logit_gamma_beta <- sapply(1:Px, function(j) {
      logit(E_theta_beta[j]) -
        (1 / 2) * drop((A_sigma2e0 + N / 2) / B_sigma2e) * crossprod(X[, j]) %*% (crossprod(m_beta[j]) + s_beta[j, j]) +
        drop(((A_sigma2e0 + N / 2) / B_sigma2e)) * t(X[, j]) %*% (A * m_beta[j] -
          X_E_w_beta[, -j] %*% (m_beta[-j] * m_beta[j] + s_beta[-j, j])
        )
    })

    E_w_beta <- sapply(1:Px, function(j) inv.logit(logit_gamma_beta[j]))

    # theta for fixed effects

    A_theta_beta <- sapply(1:Px, function(j) E_w_beta[j] + A_theta_beta0)

    B_theta_beta <- sapply(1:Px, function(j) 1 - E_w_beta[j] + B_theta_beta0)

    E_theta_beta <- sapply(1:Px, function(j) A_theta_beta[j] / (A_theta_beta[j] + B_theta_beta[j]))

    #  E_theta_beta <- rep(theta_beta, Px)

    # Omega
    omega <- tcrossprod(E_w_beta) + hadamard.prod(diag(E_w_beta), (diag(1, Px) - diag(E_w_beta)))

    XGammaBeta <- X %*% diag(E_w_beta) %*% m_beta

    #----------------------# background #----------------------#

    B0 <- yc - XGammaBeta - rowSums(ZGammaU[, , drop = F])

    # Beta

    s_b <- solve(drop(((A_sigma2e0 + N / 2) / B_sigma2e)) * tglg +
      drop((A_sigma2b0 + DIM0 / 2) / B_sigma2b) * diag(1, DIM0))

    m_b <- drop(((A_sigma2e0 + N / 2) / B_sigma2e)) * s_b %*% t(gamma_lambda_b) %*% B0

    A_sigma2b <- (1 / 2) * N + A_sigma2b0

    mtm_b <- crossprod(m_b)

    tr_s_b <- tr(s_b)

    B_sigma2b <- (1 / 2) * (mtm_b + tr_s_b) + B_sigma2b0

    # background g
    g <- gamma_lambda_b %*% m_b


    #----------------------# error #----------------------#

    # sigma B of error

    #  CtC1  <- 2 * t(yc) %*% ( X %*% diag(E_w_beta) %*% m_beta
    #  )

    CtC1 <- X %*% diag(E_w_beta) %*% m_beta

    CtC2 <- sapply(1:M, function(m) drop(E_w_u[[m]]) * Z[[m]] %*% m_u[[m]])

    CtC3 <- g
    ####
    CtC4 <- tr(hadamard.prod(crossprod(X), omega) %*% s_beta)

    CtC5 <- sapply(1:M, function(m) tr(drop(E_w_u[[m]]) * crossprod(Z[[m]]) %*% s_u[[m]]))

    CtC6 <- tr(tglg %*% s_b)

    B_sigma2e <- B_sigma2e0 + (1 / 2) * (crossprod(yc - CtC1 - rowSums(CtC2) - CtC3) +
      sum(CtC4) +
      sum(CtC5) +
      sum(CtC6)
    )


    #---------------- calculating the sigma --------------------------

    ##########################################

    post_sigma2beta <- B_sigma2beta / (A_sigma2beta - 1)

    post_sigma2u <- sapply(1:M, function(m) B_sigma2u[[m]] / (A_sigma2u[m] - 1))

    post_sigma2b <- B_sigma2b / (A_sigma2b - 1)

    post_sigma2e <- B_sigma2e / (A_sigma2e - 1)

    ##########################################


    #################################################
    #                                               #
    #          2.  Iteration (lower Bound)          #
    #                                               #
    #################################################

    # Compute lower bound

    #--------------------p--------------------

    # lb_py <- -N/2 * log(2 * pi) - 1/2 * (A_sigma2e/B_sigma2e) * ( crossprod(yc - CtC1 - sum(CtC2) ) +
    #                                                                 sum(CtC3) + sum(CtC4))

    lb_py <- -N / 2 * log(2 * pi) - 1 / 2 * (A_sigma2e / B_sigma2e) * (crossprod(yc - rowSums(CtC2)) +
      sum(CtC3))

    # fix

    lb_pbeta <- N / 2 * log(2 * pi) - 1 / 2 * (A_sigma2beta / B_sigma2beta) * (crossprod(m_beta) + tr(s_beta))

    lb_psigma2_beta <- A_sigma2beta0 * log(B_sigma2beta0) - lgamma(A_sigma2beta0) -
      (A_sigma2beta0 + 1) * (B_sigma2beta / (A_sigma2beta - 1)) + B_sigma2beta0 * (A_sigma2beta / B_sigma2beta)

    lb_pgamma_beta <- sapply(1:Px, function(j) {
      E_w_beta[j] * (digamma(A_theta_beta[j]) - digamma(A_theta_beta[j] + B_theta_beta[j])) +
        (1 - E_w_beta[j]) * (digamma(B_theta_beta[j]) - digamma(A_theta_beta[j] + B_theta_beta[j]))
    })

    lb_ptheta_beta <- sapply(1:Px, function(j) {
      (A_theta_beta0 - 1) * (digamma(A_theta_beta[j]) - digamma(A_theta_beta[j] + B_theta_beta[j])) +
        (B_theta_beta0 - 1) * (digamma(B_theta_beta[j]) - digamma(A_theta_beta[j] + B_theta_beta[j]))
    })

    # random

    lb_pu <- sapply(1:M, function(m) -N / 2 * log(2 * pi) - N / 2 * (A_sigma2u[m] / B_sigma2u[[m]]) * (crossprod(m_u[[m]]) + tr(s_u[[m]])))

    lb_psigma2_u <- sapply(1:M, function(m) {
      A_sigma2u0[m] * log(B_sigma2u0[m]) - lgamma(A_sigma2u0[m]) -
        (A_sigma2u0[m] + 1) * (log(B_sigma2u[[m]]) - digamma(A_sigma2u[m])) +
        B_sigma2u0[m] * (A_sigma2u[m] / B_sigma2u[[m]])
    })

    lb_pgamma_u <- sapply(1:M, function(m) {
      E_w_u[[m]] * (digamma(A_theta_u[m]) - digamma(A_theta_u[m] + B_theta_u[m])) +
        (1 - E_w_u[[m]]) * (digamma(B_theta_u[m]) - digamma(A_theta_u[m] + B_theta_u[m]))
    })

    lb_ptheta_u <- sapply(1:M, function(m) {
      (A_theta_u0 - 1) * (digamma(A_theta_u[m]) - digamma(A_theta_u[m] + B_theta_u[m])) +
        (B_theta_u0 - 1) * (digamma(B_theta_u[m]) - digamma(A_theta_u[m] + B_theta_u[m]))
    })

    # background

    lb_pb <- N / 2 * log(2 * pi) - 1 / 2 * (A_sigma2b / B_sigma2b) * (crossprod(m_b) + tr(s_b))

    lb_psigma2_b <- A_sigma2b0 * log(B_sigma2b0) - lgamma(A_sigma2b0) -
      (A_sigma2b0 + 1) * (B_sigma2b / (A_sigma2b - 1)) + B_sigma2b0 * (A_sigma2b / B_sigma2b)

    # error

    lb_psigma2_e <- A_sigma2e0 * log(B_sigma2e0) - lgamma(A_sigma2e0) -
      (A_sigma2e0 + 1) * (log(abs(B_sigma2e)) - digamma(A_sigma2e)) -
      B_sigma2e0 * (A_sigma2e / B_sigma2e)

    #--------------------q--------------------

    # fix

    lb_qsigma2_beta <- -A_sigma2beta - log(B_sigma2beta) * lgamma(A_sigma2beta) +
      (1 + A_sigma2beta) * digamma(A_sigma2beta)

    #  lb_qbeta <- - 1/2 * log(exp(determinant(s_beta)$modulus[1]) + 1e-10) - 1/2 * log(2 * pi)

    lb_qbeta <- -1 / 2 * determinant(s_beta)$modulus[1] - 1 / 2 * log(2 * pi)

    lb_qgamma_beta <- sapply(1:Px, function(j) {
      E_w_beta[j] * log(E_w_beta[j] + 1e-10) +
        (1 - E_w_beta[j]) * log(1 - E_w_beta[j] + 1e-10)
    })

    lb_qtheta_beta <- sapply(1:Px, function(j) {
      lgamma(A_theta_beta[j] + B_theta_beta[j]) - lgamma(A_theta_beta[j]) - lgamma(B_theta_beta[j]) +
        (A_theta_beta[j] - 1) * (digamma(A_theta_beta[j]) - digamma(A_theta_beta[j] + B_theta_beta[j])) +
        (B_theta_beta[j] - 1) * (digamma(B_theta_beta[j]) - digamma(A_theta_beta[j] + B_theta_beta[j]))
    })

    # random

    lb_qsigma2_u <- sapply(1:M, function(m) {
      -A_sigma2u[m] - (log(B_sigma2u[[m]]) + log(lgamma(A_sigma2u[m]))) + (1 + A_sigma2u[m]) *
        digamma(A_sigma2u[m])
    })

    #  lb_qu <- sapply(1:M, function(m) - 1/2 * log(exp(determinant(s_u[[m]])$modulus[1]) + 1e-10) -1/2 * log(2 * pi) )

    lb_qu <- sapply(1:M, function(m) -1 / 2 * determinant(s_u[[m]])$modulus[1] - 1 / 2 * log(2 * pi))
    #  lb_qu <- sapply(1:M, function(m) - 1/2 * log(abs(determinant(s_u[[m]])$modulus[1]) + 1e-10) -1/2 * log(2 * pi) )

    #  lb_qu <- sapply(1:M, function(m) - 1/2 * log(det(s_u[[m]]) + 1e-10) -1/2 * log(2 * pi) )
    #  lb_qu <- sapply(1:M, function(m) - 1/2 * log(ifelse(is.infinite(det(s_u[[m]])), 1e20, det(s_u[[m]])) + 1e-10) -1/2 * log(2 * pi) )

    lb_qgamma_u <- sapply(1:M, function(m) {
      E_w_u[[m]] * log(E_w_u[[m]] + 1e-10) +
        (1 - E_w_u[[m]]) * log(1 - E_w_u[[m]] + 1e-10)
    })

    lb_qtheta_u <- sapply(1:M, function(m) {
      lgamma(A_theta_u[m] + B_theta_u[m]) - lgamma(A_theta_u[m]) - lgamma(B_theta_u[m]) +
        (A_theta_u[m] - 1) * (digamma(A_theta_u[m]) - digamma(A_theta_u[m] + B_theta_u[m])) +
        (B_theta_u[m] - 1) * (digamma(B_theta_u[m]) - digamma(A_theta_u[m] + B_theta_u[m]))
    })

    # background

    lb_qsigma2_b <- -A_sigma2b - log(B_sigma2b) * lgamma(A_sigma2b) +
      (1 + A_sigma2b) * digamma(A_sigma2b)

    #  lb_qb <- - 1/2 * log(exp(determinant(s_b)$modulus[1]) + 1e-10) - 1/2 * log(2 * pi)

    lb_qb <- -1 / 2 * determinant(s_b)$modulus[1] - 1 / 2 * log(2 * pi)

    # error
    lb_qsigma2_e <- -A_sigma2e - (log(abs(B_sigma2e)) + log(lgamma(A_sigma2e))) + (1 + A_sigma2e) * digamma(A_sigma2e)

    #---------------------L--------------------

    L[i] <- lb_py +
      lb_pbeta + lb_psigma2_beta + sum(lb_pgamma_beta) + sum(lb_ptheta_beta) +
      sum(lb_pu) + sum(lb_psigma2_u) + sum(lb_pgamma_u) + sum(lb_ptheta_u) +
      lb_pb + lb_psigma2_b +
      lb_psigma2_e -
      lb_qsigma2_beta - lb_qbeta - sum(lb_qgamma_beta) - sum(lb_qtheta_beta) -
      sum(lb_qu) - sum(lb_qsigma2_u) - sum(lb_qgamma_u) - sum(lb_qtheta_u) -
      lb_qsigma2_b - lb_qb -
      lb_qsigma2_e

    log_delta <- abs(L[i] - L[i - 1])
    # Brute force
    #  L0 <- L
    #  L0 <- L0[is.finite(L0)]

    #  if(length(unique(L0)) == length(L0) && i != max_iter){
    #    log_delta <- abs(L[i] - L[i - 1])
    #  }else if(length(unique(L0)) != length(L0) || i == max_iter){
    #    log_delta <- 1e-10
    #  }

    #################################################
    #                                               #
    #           3. Output Iteration Results         #
    #                                               #
    #################################################

    # Show VB difference
    if (verbose && i %% telliter == 0) {
      cat(
        "Iteration:\t", i, "\tLB:\t", L[i], "\tDelta:\t",
        log_delta, "\n"
      )
    }

    # Check for convergence
    if (log_delta < epsilon_conv) {
      if (verbose && i %% telliter != 0) {
        cat(
          "Iteration:\t", i, "\tLB:\t", L[i], "\tDelta:\t",
          log_delta, "\n"
        )
      }

      cat("[2] VB coverged!\n")
      break
    }
    # Check if VB converged in the given maximum iterations
    if (i == max_iter) {
      message("[2] VB did not converge!\n")
    }
  }

  if (log_delta < epsilon_conv) {
    obj <- structure(list(
      X = X, Z = Z, genotype0 = tot_snp_only, gamma_lambda = gamma_lambda, B_sigma2beta = B_sigma2beta,
      index_test = index_test, B_sigma2e = B_sigma2e, B_sigma2u = B_sigma2u, post_sigma2b = post_sigma2b,
      post_sigma2u = post_sigma2u, post_sigma2beta = post_sigma2beta, post_sigma2e = post_sigma2e,
      N = N, M = M, Px = Px, DIM = DIM, DIM0 = DIM0, Delta = log_delta, E_theta_u = E_theta_u,
      E_theta_beta = E_theta_beta, s_beta = s_beta, m_b = m_b, s_b = s_b,
      m_beta = m_beta, r = M, m_u = m_u, s_u = s_u, y = yc, scale = scale, center = center,
      E_w_beta = E_w_beta, E_w_u = E_w_u, Px = Px, index_fixed_marker = index_fixed_marker2,
      L = L[2:i]
    ), UD = UListD, class = "vb")
    return(obj)
  }
}
