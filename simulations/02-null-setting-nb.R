########################################
########### Null simulations ###########
########################################

source("./00-functions.R")
library(kdca)
library(tidyverse)
library(digest)
library(parallel)

run_study <- function(N,
                      set_size,
                      add.covariate,
                      continuous,
                      prop_null,
                      vqtl,
                      id,
                      disp,
                      rep,
                      seed) {
  set.seed(seed)
  
  # generate primary variable
  if (continuous) {
    type_x <- "continuous"
    X <- runif(N, 0.0, 2.0)
  } else {
    type_x <- "categorical"
    X <- rbinom(n = N, size = 1, prob = 0.5)
  }
  
  # covariate
  if (add.covariate) {
    Z <- runif(N, 0.0, 2.0)
  }

  # generate data set of 10,000 genes from NB distribution
  sim <- generate_data_nb(X,
                          m = 10000,
                          n = N,
                          disp = rep(disp, 10000),
                          vqtl = vqtl,
                          lib.size = sample(rep(c(2, 10), N / 2) * 10^6),
                          baseline = exp(exp(rnorm(10000,
                                                   mean = 1.7,
                                                   sd = 0.15))),
                          pi0 = prop_null,
                          adj = Z)

  out <- NULL
  for (i in 1:1000) {
    set.seed(seed + i)
    
    # randomly create a pathway 
    id <- sample(1:10000, size = set_size, replace = FALSE)
    primary <- sim$X
    trans <- log2(sim$Y[id,] + 2)
    expr <- t(trans)

    # apply KDCA w/ library size as covariate
    p <- kdca(x = as.matrix(primary),
              mean_adjust = model.matrix(~Z + log2(colSums(sim$Y))),
              y = expr,
              type_x = type_x,
              perm.its = 1000)
    
    df <- data.frame(method = p$kernel,
                     id = i,
                     pvalues = c(p$p),
                     continuous = continuous,
                     rep = rep,
                     set_size = set_size,
                     seed = seed)

    out <- rbind(out,df)
  }
  
  saveRDS(out, file = paste0("./data/02-null-setting-nb/results", seed, ".rds"))
  rm(out)
  return(NULL)
}

# study design
design <- expand.grid(N = 300,
                      rep = 1:20,
                      id =  1000,
                      continuous = c(FALSE, TRUE),
                      add.covariate = TRUE,
                      vqtl = TRUE,
                      disp = 0.25,
                      set_size = c(10, 50),
                      prop_null = 0.0)

# study design
design <- design %>%
  dplyr::group_by(N,  id, rep, vqtl,  disp, add.covariate, continuous, prop_null, set_size) %>%
  dplyr::mutate(seed = readBin(digest(c(N, rep, disp, id, vqtl, add.covariate, continuous, prop_null, set_size), raw = TRUE), "integer"))

cl <- makeCluster(10, type = "PSOCK")
clusterExport(cl, varlist = c("design", "run_study", "generate_geneset_nb"))
out <- parLapply(cl, 1:nrow(design), function(ii) {
  library(tidyverse)
  library(ssizeRNA)
  library(kdca)
  library(limma)
  return(run_study(design[ii,]$N,
                   set_size = design[ii,]$set_size,
                   continuous = design[ii,]$continuous,
                   add.covariate = design[ii,]$add.covariate,
                   prop_null = design[ii,]$prop_null,
                   vqtl = design[ii,]$vqtl,
                   disp = design[ii,]$disp,
                   seed = design[ii,]$seed,
                   rep = design[ii,]$rep,
                   id = design[ii,]$id))
})
stopCluster(cl)
