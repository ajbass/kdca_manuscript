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
                      num_null,
                      continuous,
                      positive,
                      vmax,
                      bmax,
                      id,
                      rep,
                      seed) {

  # generate primary variable
  set.seed(seed)
  if (continuous) {
    type_x <- "continuous"
    X <- runif(N, 0.0, 2.0)
  } else {
    type_x <- "categorical"
    X <- rbinom(n = N, size = 1, prob = 0.5)
  }

  # covariate
  Z <- runif(N, 0.0, 2.0)

  out <- NULL
  for (i in 1:id) {
    # generate pathway
    set.seed(seed + i)
    dat <- generate_geneset(X,
                            Z,
                            N = N,
                            set_size = set_size,
                            num_null = num_null,
                            vmax = vmax,
                            bmax = bmax,
                            positive = positive)
    
    # KDCA approach w/ permutation
    df_perm <- kdca(x = as.matrix(dat$X),
                    y = dat$Y,
                    mean_adjust = model.matrix(~scale(Z)),
                    type_x = ifelse(continuous, "continuous", "categorical"),
                    perm.its = 1000)
    df_perm$null <- "permutation"

    # KDCA approach w/ theoretical null 
    df_vtheo <- kdca:::run_dkat(x = as.matrix(dat$X),
                                y = dat$Y,
                                mean_adjust = model.matrix(~scale(Z)),
                                type_x = ifelse(continuous, "continuous", "categorical"))
    df_vtheo$null <- "theoretical_Var"
    rm(dat)

    df <- rbind(df_perm, df_vtheo)
    df$id <- i
    df$seed <- seed
    df$continuous <- continuous
    df$rep <- rep
    df$set_size <- set_size
    out <- rbind(out, df)
  }
  
  saveRDS(out, file = paste0("./data/01-null-setting/results", seed, ".rds"))
  return(out)
}

# study design
design <- expand.grid(N = 300,
                      rep = 1:20,
                      id = 1000,
                      continuous = c(TRUE, FALSE),
                      positive = TRUE,
                      bmax = 0.5,
                      vmax = 0.10,
                      set_size = c(10, 50),
                      prop_null = 1)

design <- design %>%
  group_by(N,  bmax, rep, vmax, positive, continuous, set_size, prop_null) %>%
  mutate(seed = readBin(digest(c(N, rep, bmax, vmax,  positive, continuous, set_size, prop_null), raw = TRUE), "integer"))

cl <- makeCluster(10, type = "PSOCK")
clusterExport(cl, varlist = c("design", "run_study", "generate_geneset"))
out <- parLapply(cl, 1:nrow(design), function(ii) {
  library(tidyverse)
  library(kdca)
  library(PearsonDS)
  return(run_study(design[ii,]$N,
                   set_size = design[ii,]$set_size,
                   num_null = round(design[ii,]$set_size * design[ii,]$prop_null),
                   continuous = design[ii,]$continuous,
                   bmax = design[ii,]$bmax,
                   positive = design[ii,]$positive,
                   vmax = design[ii,]$vmax,
                   seed = design[ii,]$seed,
                   id = design[ii,]$id,
                   rep = design[ii,]$rep))
})
stopCluster(cl)
