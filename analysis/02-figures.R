#####################################
##### Create manuscript figures #####
#####################################

# Calculate number of detections as function of FDR
fdr_discoveries <- function(pvalues, fdr = seq(0.0001, 1, 0.0001)) {
  d <- rep(NA, length(fdr))
  pvalues <- pvalues[!is.na(pvalues)]
  q <- p.adjust(pvalues, "BH")
  for (i in 1:length(fdr)) {
    d[i] <- sum(q <= fdr[i])
  }
  return(data.frame(detections = d,
                    fdr = fdr))
}

adjustedP <- function(pvalues, id, description, test, theo) {
  q <- p.adjust(pvalues, "BH")
  return(data.frame(bh = q, id = id,
                    description=description,
                    tested_size = test,
                    theoretical_size = theo))
}

# combined files
fnames <- list.files("./data/")
df <- NULL
for (f in fnames) {
  out <- readRDS(paste0("./data/03-kdca-thca/", f))
  df <- rbind(df, out)
}

# Reorder for plotting
df$method <- factor(df$method, labels = c("Eigengene",  "KDCA (Combined)", "KDCA (Gaussian)", "KDCA (Linear)", "KDCA (Projection)"))
df$method <- factor(df$method, levels = c("KDCA (Combined)", "KDCA (Linear)",  "KDCA (Gaussian)",  "KDCA (Projection)", "Eigengene"))
cbbPalette <- c("darkred", "#000000", "#009E73", "#0072B2",  "darkgrey", "#0072B2", "#D55E00", "#CC79A7")

# pathway name
name <- sapply(df$description,
               FUN = function(d) unlist(as.character(str_split(d, "_")[[1]][1])))
df$name <- name

# plotting prep
df0 <- df %>% 
  dplyr::select(description, type, theoretical_size, tested_size, pvalues, method) %>%
  filter(!is.na(pvalues)) %>%
  spread(method, pvalues) %>% arrange(`KDCA (Combined)`)
colnames(df0)[1:4] <- c("Pathway", "Risk factor", "Size", "Tested size")
df0$`Risk factor`[df0$`Risk factor` == "categorical"] <- "Age at diagnosis"
df0$`Risk factor`[df0$`Risk factor` == "continuous"] <- "BRAF mutation"
df0$Pathway <- str_split(df0$Pathway, "_", simplify = T)[,2]
df0 <- df0 %>%
  dplyr::select(-`KDCA (Linear)`, -`KDCA (Gaussian)`, -`KDCA (Projection)`) 

# table
library(xtable)
print(xtable(df0, digits = 4),  include.rownames=FALSE)

# discoveries figure
p1 <- df %>%
  group_by(method, type) %>%
  mutate(type2 = ifelse(type == "categorical", "BRAF mutation", "Age at Diagnosis")) %>%
  filter(!is.na(pvalues)) %>%
  group_by(method, type2) %>%
  dplyr:::do(fdr_discoveries(.$pvalues)) %>%
  filter(fdr < 1) %>%
  ggplot(aes(x = fdr, y = detections, color = method)) + theme_bw() + xlab("FDR threshold") +
  ylab("Discoveries") +
  geom_line() +
  scale_colour_manual(name = "Method", values = cbbPalette) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_grid(~type2)

ggsave(p1,width = 6,height = 2.5, filename = "./figures/04-power-comparison.pdf")

# KDCA combined vesus eigengene
p1 <- df %>% filter(method %in% c("KDCA (Combined)", "Eigengene")) %>%
  group_by(method, type) %>%
  mutate(type2= ifelse(type == "categorical", "BRAF mutation", "Age at Diagnosis")) %>%
  filter(!is.na(pvalues)) %>%
  group_by(method, type2) %>%
  dplyr:::do(f(.$pvalues, fdr)) %>%
  filter(fdr < 1) %>%
  ggplot(aes(x = fdr, y = detections, color = method)) + theme_bw() + xlab("FDR threshold") +
  ylab("Discoveries") +
  geom_line(size =.8) +
  scale_colour_manual(name = "Method", values = cbbPalette) +# ggtitle("COAD")  +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_grid(~ type2)

ggsave(p1,width = 6,height = 2.5, filename = "./figures/04-power-comparison-subset.pdf")

# Ranks the genes
out <- df %>%
  dplyr::filter(!is.na(pvalues)) %>% filter(method == "KDCA (Combined)", type == "continuous") %>%
  dplyr::group_by(method, type,name) %>%
  dplyr::do(adjustedP(.$pvalues, .$id, .$description, .$tested_size, .$theoretical_size)) %>% #data.frame(bh = p.adjust(.$pvalues, "BH"), id = .$id)) %>%
  filter(bh <= 0.15)


# Total discoveries
df %>%
  dplyr::filter(!is.na(pvalues)) %>%
  dplyr::group_by(method, type, name) %>%
  dplyr::do(adjustedP(.$pvalues, .$id, .$description, .$tested_size, .$theoretical_size)) %>% 
  filter(bh <= 0.15) %>%
  dplyr::group_by(method, type) %>%
  dplyr::summarise(length(bh))

# Average set size
df %>%
  dplyr::filter(!is.na(pvalues)) %>%
  dplyr::group_by(method, type, name) %>%
  dplyr::do(adjustedP(.$pvalues, .$id, .$description, .$tested_size, .$theoretical_size)) %>% #data.frame(bh = p.adjust(.$pvalues, "BH"), id = .$id)) %>%
  dplyr::group_by(method, type) %>%
  dplyr::summarise(mean(tested_size))

# calculate overlap for pathways
all_gene_sets = msigdbr(species = "human")
gene_set <- all_gene_sets %>% filter(gs_cat %in% c("C2")) %>%
  filter(gs_subcat %in% c( "CP:BIOCARTA"))
genes <-  rownames(cpm)

pids <- unique(gene_set$gs_id)
sets <- list()

# Get pathways to determine overlap
for (i in 1:length(pids)) {
  print(i)
  tmp_id <- pids[i]
  gs <- gene_set %>% filter(gs_id %in% tmp_id)
  id <- rownames(dataPrep) %in% gs$ensembl_gene
  subset <- t(dataPrep[id,, drop = F])

  set_size = ncol(subset)

  if (set_size >= 5 & set_size <= 150) {
    genes <- colnames(subset)
    df <- data.frame(genes = genes, geneset = tmp_id)
    sets <- rbind(sets, df)
  }
}

# Combine
id <- out$id[out$type == "continuous"] #out$type == "categorical"
tt <- gene_set %>% filter(gs_id %in% id)
tt <- tt[tt$human_ensembl_gene %in% rownames(cpm),]
ll <- list()
for (i in 1:length(id)) {
  ll[id[i]] <- list(tt$human_ensembl_gene[tt$gs_id== id[i]])
}

overlaps <- sapply(ll, function(g1)
  sapply(ll, function(g2) {round(length(intersect(g1, g2)) / length(g2) * 100)}))

# Genes that overlap the most with significant pathways
tt %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::summarise(freq = length(unique(gs_id))) %>%
  dplyr::arrange(-1*freq)
