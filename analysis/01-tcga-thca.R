#####################################
######### Run TCGA analysis #########
#####################################

# Load packages
library(recount3)
library(tidyverse)
library(edgeR)
library(msigdbr)
library(sva)
library(parallel)
library(TCGAbiolinks)
options(recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release")

# download thca data using recount3
human_projects <- available_projects()
proj_info <- subset(human_projects,
                    project == "THCA" & project_type == "data_sources")
rse_coad <- create_rse(proj_info)
assay(rse_coad, "counts") <- transform_counts(rse_coad,
                                              by = "auc",
                                              round = FALSE)

# Remove overlaps
expr <- assay(rse_coad, "counts")
ov_gene_by_exon <- countOverlaps(rowRanges(rse_coad))
rdata <- rowData(rse_coad)
rdata <- rdata[ov_gene_by_exon==1, ]
expr <- expr[ov_gene_by_exon==1,]

# Select primary tumor samples, non-FFPE samples, age, BRAF
cdata <- colData(rse_coad)
TP <- cdata$tcga.cgc_sample_sample_type == "Primary Tumor"
age <- cdata$tcga.gdc_cases.diagnoses.age_at_diagnosis
sex <- cdata$tcga.gdc_cases.demographic.gender

ffpe <- cdata$tcga.cgc_sample_is_ffpe
race <- cdata$tcga.xml_race_list
subtype <- TCGAquery_subtype("THCA")
subtype <- data.frame(patient = cdata$tcga.xml_bcr_patient_barcode) %>%
  left_join(subtype)
subtype <- subtype$BRAF

id <- (!is.na(age)) & (!is.na(subtype)) & !(ffpe == "YES") & (TP) 
id[is.na(id)] <- FALSE
expr  <- expr[,id]
cdata <- cdata[id,]
age <- age[id]
subtype <- subtype[id]

# Calculate normalization factors and logCPM values
dataPrep <- expr
dge <- DGEList(counts = dataPrep)
dge <- calcNormFactors(dge, method = "TMM")
dge$genes <- dge$genes
dge$targets <- dge$samples
lib.size <- dge$samples$lib.size * dge$samples$norm.factors
counts <- dge$counts
cpm <- t(log2(t(counts + 0.5) / (lib.size + 1) * 1e+06))

# Filter low counts and select protein coding regions
prop_expressed = rowMeans(cpm >= 1)
id <- prop_expressed  >= 0.5
rdata <- rdata[id,]
dataPrep <- dataPrep[id,]
cpm <- cpm[id,]
keep <- rdata$gene_type %in% c("protein_coding")
dataPrep <- dataPrep[keep,]
cpm <- cpm[keep,]
rdata <- rdata[keep,]

# Select top 5,000 genes with highest variance
std.dev <- rowSds(cpm)
id <- rank(-(std.dev)) <= 5000
dataPrep <- dataPrep[id,]
cpm <- cpm[id,]
rdata <- rdata[id,]
cdata <- data.frame(external_id = colnames(cpm)) %>%
  left_join(data.frame(cdata))

# processing gene names
rownames(cpm) <- gsub('\\.[0-9]*', '', rownames(cpm))
rownames(cpm) <- rownames(cpm)
rownames(dataPrep) <- rownames(cpm)
rownames(rdata) <- rownames(cpm)

# Biocarta pathways
all_gene_sets = msigdbr(species = "Homo sapiens")
gene_set <- all_gene_sets %>% filter(gs_cat %in% c("C2")) %>%
  filter(gs_subcat %in% c("CP:BIOCARTA"))
genes <-  rownames(cpm)
pids <- unique(gene_set$gs_id)

# Estimate batch effects top PCs
expr <- log2(dataPrep + 2)
mod <-  matrix(1, nrow = dim(expr)[2], ncol = 1)
nsv <- num.sv(t(scale(t(expr), scale = F)), mod, method = "be")
dat <- scale(t(expr), scale = F)
ss <- svd(dat)

batch <- ss$u[, 1:nsv]
braf <- subtype
age <- scale(cdata$tcga.gdc_cases.diagnoses.age_at_diagnosis)

# run KDCA in parallel
run_kdca <- function(expr, pids, gene_set, batch, braf, age) {
  # select gene set
  gs <- gene_set %>% filter(gs_id %in% pids)
  id <- rownames(expr) %in% gs$ensembl_gene
  subset <- t(expr[id,, drop = F])
  set_size <- ncol(subset)

  # Focus on pathway sizes b/w [5, 150]
  if (set_size >= 5 & set_size <= 150) {
    # KDCA
    p_kdca <- kdca::kdca(x = as.matrix(braf),
                         y = scale(subset),
                         type_x = "categorical",
                         mean_adjust = model.matrix(~batch),
                         perm.its = 10000)$pvalues
    p_kdca2 <- kdca::kdca(x = as.matrix(age),
                          y = scale(subset),
                          type_x = "continuous",
                          mean_adjust = model.matrix(~batch),
                          perm.its = 10000)$pvalues

    # Eigengene
    p_eg <- kdca:::eigengene(x = as.matrix(braf),
                             y = scale(subset),
                             type_x = "categorical",
                             mean_adjust = model.matrix(~batch),
                             perm.its = 10000)$pvalues
    p_eg2 <- kdca:::eigengene(x = as.matrix(age),
                              y = scale(subset),
                              type_x = "continuous",
                              mean_adjust = model.matrix(~batch),
                              perm.its = 10000)$pvalues

    df <- data.frame(method = rep(c("KDCA (Combined)", "KDCA_Linear", "KDCA_Gaussian", "KDCA_Projection", "Eigengene"), 2),
                     pvalues = c(p_kdca, p_eg, p_kdca2, p_eg2),
                     type = rep(c("categorical", "continuous"), each = 5),
                     id = pids,
                     theoretical_size = dim(gs)[1],
                     tested_size = set_size,
                     cat = unique(gs$gs_cat),
                     description = gs$gs_name[1])
  } else {
    p_eg <- p_eg2 <- NA
    p_kdca <- p_kdca2 <- c(NA, NA, NA, NA)

    df <- data.frame(method = rep(c("KDCA (Combined)", "KDCA_Linear", "KDCA_Gaussian", "KDCA_Projection", "Eigengene"), 2),
                     pvalues = c(p_kdca, p_eg, p_kdca2, p_eg2),
                     type = rep(c("categorical", "continuous"), each = 5),
                     id = pids,
                     theoretical_size = dim(gs)[1],
                     tested_size = set_size,
                     cat = unique(gs$gs_cat),
                     description = gs$gs_name[1])
  }
  saveRDS(df, file = paste0("./data/", pids, ".rds"))
  return(df)
}

cl <- makeCluster(10, type = "PSOCK")
clusterExport(cl, varlist = c("gene_set", "batch", "run_kdca", "braf", "age", "expr", "pids"))
out <- parLapply(cl, 1:length(pids), function(ii) {
  library(tidyverse)
  library(kdca)
  set.seed(ii)
  return(run_kdca(expr,
                  pids[ii],
                  gene_set,
                  batch,
                  braf,
                  age))
})
stopCluster(cl)
