########################################
########## Time comparisons ############
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
                      id,
                      seed) {
  # generate primary variable
  set.seed(seed)
  if (continuous) {
    X <- runif(N, 0.0, 2.0)
  } else {
    X <- rbinom(n = N, size = 1, prob = 0.5)
  }
  
  # covariate
  Z <- runif(N, 0.0, 2.0)
  
  # generate gene set variable
  set.seed(seed + id)
  dat <- TRUE
  while (isTRUE(dat)) {
    # Ensure the simulated covariance is positive definitive
    dat <- tryCatch({
      generate_geneset(X,
                       Z,
                       N = N,
                       set_size = set_size,
                       num_null = num_null,
                       vmax = vmax,
                       bmax = bmax,
                       positive = positive)
    }, error = function(e) {TRUE})
  }
  t1 <- proc.time()
  df_perm <- kdca(x = as.matrix(dat$X),
                  y = dat$Y,
                  mean_adjust = model.matrix(~scale(Z)),
                  type_x = ifelse(continuous, "continuous", "categorical"),
                  perm.its = 1000)
  t2 <- proc.time() - t1
  time <- t2[3]
  
  rm(dat)
  return(data.frame(time = time))
}

# study design
design <- expand.grid(N = c(100, 300, 500),
                      rep = c(1:50),
                      id = 1:1,
                      continuous = c(FALSE, TRUE),
                      positive = c(TRUE),
                      bmax = c(0.5),
                      vmax = c(0.1),
                      set_size = c(10, 50),
                      prop_null = 0)

design <- design %>%
  group_by(N,  bmax, id, rep, vmax, positive, continuous, set_size, prop_null) %>%
  dplyr::mutate(seed = readBin(digest(c(N, rep, id, bmax, vmax,  positive, continuous, set_size, prop_null), raw = TRUE), "integer"))

cl <- makeCluster(1, type = "PSOCK")
clusterExport(cl, varlist = c("design", "run_study", "generate_geneset"))
out <- parLapply(cl, 1:nrow(design), function(ii) {
  library(tidyverse)
  library(kdca)
  return(run_study(design[ii,]$N,
                   set_size = design[ii,]$set_size,
                   num_null = round(design[ii,]$set_size * design[ii,]$prop_null),
                   continuous = design[ii,]$continuous,
                   bmax = design[ii,]$bmax,
                   positive = design[ii,]$positive,
                   vmax = design[ii,]$vmax,
                   seed = design[ii,]$seed,
                   id = design[ii,]$id))
})
stopCluster(cl)

out <- dplyr::bind_rows(out, .id = "column_label")
design$column_label <- as.character(1:nrow(design))
out <- design %>% right_join(out)

saveRDS(out, file = "./data/06-time.rds")
