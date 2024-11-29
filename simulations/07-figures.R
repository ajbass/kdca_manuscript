########################################
############### Figures ################
########################################

####### Null simulations ########
library(tidyverse)

# load null data simulations
ff <- list.files("./data/01-null-setting/")
cat <- NULL
for (f in ff) {
  o <- readRDS(paste0("./data/01-null-setting/", f))
  cat <- rbind(o, cat)
}

labs <- c("Continuous", "Categorical")
names(labs) <- c("TRUE", "FALSE")

# Type 1 error rate plot
p <- cat %>%
  filter(null == "permutation",
         kernel %in% c("Combined", "Gaussian", "Linear", "Projection")) %>%
  group_by(rep, kernel, continuous, null, set_size) %>%
  summarise(fpr = mean(pvalues <= 0.05)) %>%
  ggplot(aes(x = kernel, y = fpr)) +
  geom_boxplot(outlier.size = 0.4) +
  geom_abline(intercept = 0.05, slope = 0, linetype = "dashed", color = "red") +
  facet_grid(continuous ~ set_size, scales = "free",  labeller = labeller(continuous = labs)) +
  xlab("") +
  ylab("False positive rate") +
  theme_bw()  +
  theme(strip.background = element_rect(fill = NA, color = "black"), axis.line = element_line(color='black'),plot.title = element_text(hjust = 0.5),
        plot.background = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(size=12, vjust = -1.7), axis.text.x = element_text(color = "black"))

ggsave(p,  width = 6.0, height = 3.5, dpi = 700, filename = "./figures/01-kdca-null-typeierror.png")

p <- cat %>%
  filter(null == "theoretical_Var",
         kernel %in% c("Gaussian", "Linear", "Projection")) %>%
  group_by(rep, kernel, continuous, null, set_size) %>%
  summarise(fpr = mean(pvalues <= 0.05)) %>%
  ggplot(aes(x = kernel, y = fpr)) +
  geom_boxplot(outlier.size = 0.4) +
  geom_abline(intercept = 0.05, slope = 0, linetype = "dashed", color = "red") +
  facet_grid(continuous ~ set_size, scales = "free", labeller = labeller(continuous = labs)) +
  xlab("") +
  ylab("False positive rate") +
  theme_bw()  + ylim(0, 0.25) +
  theme(strip.background = element_rect(fill = NA, color = "black"), axis.line = element_line(color='black'),plot.title = element_text(hjust = 0.5),
        plot.background = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(size=12, vjust = -1.7), axis.text.x = element_text(color = "black"))

ggsave(p,  width = 6.0, height = 3.5, dpi = 700, filename = "./figures/01-kdca-null-typeierror2.png")

# load null NB data simulations
ff <- list.files("./data/01-null-setting-nb/")
cat <- NULL
for (f in ff) {
  o <- readRDS(paste0("./data/01-null-setting-nb/", f))
  o <- o %>% select(method, id, pvalues, continuous, rep, set_size, seed)
  cat <- rbind(o, cat)
}

cat$continuous <- factor( cat$continuous, labels = c("Categorical", "Continuous"))

# Type 1 error rate plot
p <- cat %>%
  filter( method %in% c("Combined", "Gaussian", "Linear", "Projection")) %>%
  group_by(rep, method, continuous,  set_size) %>%
  sample_n(1000) %>%
  summarise(fpr = mean(pvalues <= 0.05)) %>%
  ggplot(aes(x = method, y = fpr)) +
  geom_boxplot(outlier.size = 0.4) +
  geom_abline(intercept = 0.05, slope = 0, linetype = "dashed", color = "red") +
  facet_grid(continuous ~ set_size, scales = "free", labeller = labeller(continuous = labs)) +
  xlab("") +
  ylab("False positive rate") +
  theme_bw()  +
  theme(strip.background = element_rect(fill = NA, color = "black"), axis.line = element_line(color='black'),plot.title = element_text(hjust = 0.5),
        plot.background = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(size=12, vjust = -1.7), axis.text.x = element_text(color = "black"))

ggsave(p,  width = 6.0, height = 3.5, dpi = 700, filename = "./figures/02-kdca-null-typeierror-rnaseq.png")

# load null multivariate data simulations
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

f <- list.files("./data/03-null-setting-multivariate/", recursive = T, full.names = T)
com <- NULL
for (file in f) {
  oo <- readRDS(file)
  com <- rbind(oo, com)
}

cat <- design %>%
  left_join(com) %>%
  filter(!is.na(pvalues))

# Type 1 error rate plot
p <- cat %>%
   filter(kernel %in% c("Combined", "Gaussian", "Linear", "Projection")) %>%
  group_by(rep, method, kernel, set_size) %>%
  summarise(fpr = mean(pvalues <= 0.05)) %>%
  ggplot(aes(x = kernel, y = fpr)) +
  geom_boxplot(outlier.size = 0.4) +
  geom_abline(intercept = 0.05, slope = 0, linetype = "dashed", color = "red") +
  facet_grid(~ set_size, scales = "free") +
  xlab("") +
  ylab("False positive rate") +
  theme_bw()  +
  theme(strip.background = element_rect(fill = NA, color = "black"), axis.line = element_line(color='black'),plot.title = element_text(hjust = 0.5),
        plot.background = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(plot.title = element_text(size=12, vjust = -1.7), axis.text.x = element_text(color = "black"))

ggsave(p,  width = 6.0, height = 2.5, dpi = 700, filename = "./figures/03-kdca-null-typeierror-multi.png")


####### Alternative simulations ########
# load data
out <- readRDS("./data/04-alternative-setting.rds")

filt <- out %>% filter(method == "KDCA")
filt2 <- out %>% filter(method != "KDCA")
filt <- filt %>% mutate(method = paste0(method, " (", kernel, ")" ))
out <- rbind(filt,
             filt2)

# rename for plotting
out2 <- out %>% filter(method %in% c("Eigengene", "KDCA (Combined)", "KDCA (Linear)", "KDCA (Gaussian)", "KDCA (Projection)"))
out2$method <- factor(out2$method, labels = c("Eigengene", "KDCA (Combined)", "KDCA (Gaussian)", "KDCA (Linear)", "KDCA (Projection)"))

cbbPalette <- c("darkred", "#000000", "#009E73", "#0072B2",  "darkgrey", "#0072B2", "#D55E00", "#CC79A7")
out2$positive <- factor(out2$positive, labels = c("Positive/Negative", "Positive"))
out2$positive <- factor(out2$positive, levels = c("Positive", "Positive/Negative"))

p <- out2 %>%
  filter(null == "permutation",
         continuous == TRUE) %>%
  dplyr::group_by(prop_null,   positive, method, bmax, set_size) %>%
  dplyr::summarise(power = mean(pvalues <= 0.05)) %>%
  ggplot(aes(x = 1- prop_null, y = power, color = method)) +
  geom_point(size =1.25) +
  geom_line() +
  scale_colour_manual(name = "Method",
                      values = cbbPalette,
                      breaks = c("KDCA (Combined)", "KDCA (Linear)", "KDCA (Gaussian)", "KDCA (Projection)", "Eigengene")) +
  facet_grid(positive~set_size, scales = "free") +
  xlab("Proportion of non-null genes") +
  ylab("Power") +
  ylim(0,1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, color = "black"), axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank())

ggsave(p,  width = 5.75, height = 3.5, dpi = 700, filename = "./figures/04-kdca-alt-cont.png")

p <- out2 %>%
  filter(null == "permutation",
         continuous == FALSE) %>%
  dplyr::group_by(prop_null,   positive, method, bmax, set_size) %>%
  dplyr::summarise(power = mean(pvalues <= 0.05)) %>%
  ggplot(aes(x = 1- prop_null, y = power, color = method)) +
  geom_point(size =1.25) +
  geom_line() +
  scale_colour_manual(name = "Method",
                      values = cbbPalette,
                      breaks = c("KDCA (Combined)", "KDCA (Linear)", "KDCA (Gaussian)", "KDCA (Projection)", "Eigengene")) +
 facet_grid(positive~set_size, scales = "free") +
  xlab("Proportion of non-null genes") +
  ylab("Power") +
  ylim(0,1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, color = "black"), axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank())

ggsave(p, width = 5.75, height = 3.5, dpi = 700, filename = "./figures/04-kdca-alt-cat.png")

# multivariate setting
cont <- readRDS("./data/05-alternative-setting-multivariate.rds")

filt <- out %>% filter(method == "KDCA")
filt2 <- out %>% filter(method != "KDCA")
filt <- filt %>% mutate(method = paste0(method, " (", kernel, ")" ))
out <- rbind(filt,filt2)

# rename for plotting
out2 <- out %>% filter(method %in% c("Eigengene", "KDCA (Combined)", "KDCA (Linear)", "KDCA (Gaussian)", "KDCA (Projection)"))
out2$method <- factor(out2$method, labels = c("Eigengene", "KDCA (Combined)", "KDCA (Gaussian)", "KDCA (Linear)", "KDCA (Projection)"))

cbbPalette <- c("darkred", "#000000", "#009E73", "#0072B2",  "darkgrey", "#0072B2", "#D55E00", "#CC79A7")
out2$positive <- factor(out2$positive, labels = c("Positive/Negative", "Positive"))
out2$positive <- factor(out2$positive, levels = c("Positive", "Positive/Negative"))

p <- out2 %>%
  filter(null == "permutation",
         continuous == TRUE) %>%
  dplyr::group_by(prop_null,   positive, method, bmax, set_size) %>%
  dplyr::summarise(power = mean(pvalues <= 0.05)) %>%
  ggplot(aes(x = 1- prop_null, y = power, color = method)) +
  geom_point(size =1.25) +
  geom_line() +
  scale_colour_manual(name = "Method",
                      values = cbbPalette,
                      breaks = c("KDCA (Combined)", "KDCA (Linear)", "KDCA (Gaussian)", "KDCA (Projection)", "Eigengene")) +
  facet_grid(positive~set_size, scales = "free") +
  xlab("Proportion of non-null genes") +
  ylab("Power") +
  ylim(0,1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, color = "black"), axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank())


ggsave(p, width = 5.75, height = 3.5, dpi = 700, filename = "./figures/05-kdca-mult-alt.png")

####### Time comparisons ########
time <- readRDS("./data/06-time.rds")
time$continuous <- factor(time$continuous, labels = c("Categorical", "Continuous"))

p1 <- time %>% 
  ggplot(aes(x = as.factor(N), y = time/60, color = as.factor(set_size), group = interaction(set_size, N))) + 
  geom_boxplot() +
  facet_grid(~continuous, scales = "free") + theme_bw() +
  xlab("Sample size") + 
  ylab("Time (m)") + expand_limits(  y = 0)+ scale_color_manual(values=c("#69b3a2",  "black")) + labs(color = "Pathway size")

ggsave(p1, width = 4, height = 3, dpi = 700, filename = "./figures/05-time.png")
