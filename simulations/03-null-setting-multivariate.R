########################################
####### Alternative simulations ########
########################################

source("./00-functions.R")
library(kdca)
library(tidyverse)
library(digest)
library(parallel)

run_study <- function(N,
                      set_size,
                      num_null,
                      continuous,
                      positive,
                      vmax,
                      bmax,
                      emax,
                      id,
                      seed) {

  dir.create(file.path("./data/03-null-setting-multivariate/", seed), showWarnings = FALSE)
  set.seed(seed)

  # generate primary variable
  X <- matrix(rbinom(n = 3 * N, prob = 0.5, size = 1), ncol = 3)  
  
  # covariate
  Z <- runif(N, 0.0, 2.0)

  # generate gene set variable
  set.seed(seed + id)
  dat <- TRUE

  while (isTRUE(dat)) {
    # Ensure the simulated covariance is positive definitive
    dat <- tryCatch({
      generate_geneset(2 * (rowSums(X)) / 3,
                       Z,
                       N = N,
                       set_size = set_size,
                       num_null = num_null,
                       vmax = vmax,
                       bmax = bmax,
                       emax = emax,
                       positive = positive)
    }, error = function(e) {TRUE})
  }

  # apply KDCA
  df_perm <- kdca(x = as.matrix(X),
                  y = dat$Y,
                  mean_adjust = model.matrix(~scale(Z)),
                  type_x = ifelse(continuous, "continuous", "categorical"),
                  perm.its = 1000)
  df_perm$method = "KDCA"
  df_perm$seed = seed
  df_perm$id = id
  df_perm$set_size = set_size

  rm(dat)
  saveRDS(df_perm, file = paste0("./data/03-null-setting-multivariate/", seed, "/", id, ".rds"))

  return(NULL)
}

# study design
design <- expand.grid(N = c(300),
                      rep = c(1:20),
                      id = 1:1000,
                      continuous = TRUE,
                      positive = FALSE,
                      bmax = 0.5,
                      emax  = 0.20,
                      vmax = 0.1,
                      set_size = c(10, 50),
                      prop_null = 1)

design <- design %>%
  group_by(N,  bmax,  rep, vmax, positive, continuous, set_size, prop_null) %>%
  dplyr::mutate(seed = readBin(digest(c(N, rep, bmax, vmax,  positive, continuous, set_size, prop_null), raw = TRUE), "integer"))

cl <- makeCluster(10, type = "PSOCK")
clusterExport(cl, varlist = c("design", "run_study", "generate_geneset"))
out <- parLapply(cl, 1:nrow(design), function(ii) {
  library(tidyverse)
  library(kdca)
  return(run_study(design[ii,]$N,
                   set_size = design[ii,]$set_size,
                   num_null = round(design[ii,]$set_size * design[ii,]$prop_null),
                   continuous = design[ii,]$continuous,
                   bmax = design[ii,]$bmax,
                   emax = design[ii,]$emax,
                   positive = design[ii,]$positive,
                   vmax = design[ii,]$vmax,
                   seed = design[ii,]$seed,
                   id = design[ii,]$id))
})
stopCluster(cl)
