#####################################
## Functions to generate pathways ###
#####################################

generate_geneset <- function(X,
                             Z,
                             N = length(X),
                             set_size,
                             num_null,
                             positive = TRUE,
                             vmax = 0.15,
                             bmax = 0.3,
                             emax = 0.2) {
  # intercept
  intercept = rnorm(set_size, sd = 1)

  # generate covariance matrix
  S <- matrix(nrow = set_size, ncol = set_size)
  diag(S) <- 1

  # effect size
  b <- runif(1, min = emax / 2, max = emax) 
  S[lower.tri(S)] <- b
  S <- t(S)
  S[lower.tri(S)] <- b

  # change direction of effect size if positive/negative
  if (!positive) {
    S[(set_size/2+1):set_size, 1:(set_size/2)] = -1 * S[(set_size/2+1):set_size, 1:(set_size/2)]
    S[1:(set_size/2), (set_size/2+1):set_size] = -1 * S[1:(set_size/2), (set_size/2+1):set_size]
  }

  # assign gene pairs with no differential co-expression
  if (num_null != 0) {
    num_cps <- choose(set_size, 2)
    df <- MESS::pairwise_combination_indices(set_size, self = FALSE)
    sid <- sample(1:set_size, replace = FALSE, size = num_null)
    sid <- apply(df, 1, FUN = function(x) sum(x %in% sid)) > 0
    df <- df[sid,, drop = F]
    S[df] <- 0
    S[cbind(df[, 2], df[, 1])] <- 0
  }

  # variance effect sizes
  evar <- runif(set_size, min = pmax(0, vmax - 0.1), max = vmax)

  # baseline correlation of genes
  bcor <- runif(n = N, min = bmax / 2, max = bmax)

  # generate gene set with mean, variance, and covariance effects
  Y <- matrix(nrow = N, ncol = set_size)
  for (i in 1:N) {
    M <- exp(1 + evar * X[i])
    M <- sqrt(tcrossprod(M))
    ICM <- (S * X[i] + bcor[i]) * M
    mu <- intercept  + Z[i] + X[i]
    diag(ICM) <- exp(1 + evar * X[i])
    Y[i,] <- MASS::mvrnorm(n = 1, mu = mu, Sigma = ICM)
  }

  return(list(Y = Y,
              X = X))
}

generate_data_nb <- function(X,
                             m = 100,
                             n = 200,
                             vqtl = FALSE,
                             baseline = rep(200, m),
                             disp = rep(0.5, m),
                             lib.size = sample(rep(c(5, 10), n/2) * 10^6),
                             pi0 = 0.0,
                             adj = NULL) {
  # initializations
  Xm <- t(matrix(X, nrow = n, ncol = m))
  effect.size <- rnorm(n = m, sd = .1)
  num.null <- round(pi0 * m)
  if (pi0 != 0) effect.size[1:num.null] <- 0

  # baseline counts + adjustment variable signal 
  baseline.counts <- matrix(baseline , ncol = 1) %*% matrix(rep(1, n), nrow = 1)
  adj.effect.size <- rnorm(n = m, sd = .1)
  adjm <- t(matrix(adj, nrow = n, ncol = m))

  # interaction to create vQTLs
  effect.size.int <- rnorm(n = m, sd = .025)
  W <-  t(matrix(rnorm(n*m), nrow = n, ncol = m))
  
  # expected counts
  mu <- exp(log(baseline.counts) + adj.effect.size *  (adjm) +  effect.size *  (Xm) + effect.size.int *  (Xm * W))

  # generate pathway using NB distribution
  phi <- matrix(rep(disp, n), ncol = n)
  counts <- matrix(MASS::rnegbin(m * n, mu, 1 / phi), nrow = m)
  
  # scale by library size
  counts <- round(t(t(counts) / colSums(counts) * lib.size))
  
  return(list(Y = counts,
              X = X,
              W = W))
}
