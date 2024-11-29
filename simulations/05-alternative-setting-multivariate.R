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
  # generate primary variable
  set.seed(seed)
  X <- matrix(rbinom(n = 3 * N, prob = 0.5, size = 1), ncol = 3)
  
  # covariate
  Z <- runif(N, 0.0, 2.0)

  # generate pathway
  set.seed(seed + id)
  dat <- TRUE
  while (isTRUE(dat)) {
    # Ensure the simulated covariance is positive definitive
    dat <- tryCatch({
      generate_geneset(2*(rowSums(X)) / 3, 
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
  df_perm$null = "permutation"

  # apply eigengene
  eigengene <- kdca:::eigengene(x = as.matrix(X),
                                y = dat$Y,
                                mean_adjust = model.matrix(~scale(Z)),
                                type_x = ifelse(continuous, "continuous", "categorical"),
                                perm.its = 1000)$pvalues

  eg <- data.frame(method = "Eigengene",
                   kernel = "Linear",
                   pvalues = eigengene,
                   null = "permutation")
  rm(dat)
  return(rbind(df_perm, eg))
}

# study design
design <- expand.grid(N = c(300),
                      rep = c(1:200),
                      id = 1:1,
                      continuous = TRUE,
                      positive = c(TRUE, FALSE),
                      bmax = 0.5,
                      emax  = 0.20,
                      vmax = 0.1,
                      set_size = c(10, 50),
                      prop_null = seq(0, 0.8, 0.2))

design <- design %>%
  group_by(N,  bmax, id, rep, vmax, positive, continuous, set_size, prop_null) %>%
  dplyr::mutate(seed = readBin(digest(c(N, rep, id, bmax, vmax,  positive, continuous, set_size, prop_null), raw = TRUE), "integer"))

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

out <- dplyr::bind_rows(out, .id = "column_label")
design$column_label <- as.character(1:nrow(design))
out <- design %>% right_join(out)

saveRDS(out, file = "./data/05-alternative-setting-multivariate.rds")
