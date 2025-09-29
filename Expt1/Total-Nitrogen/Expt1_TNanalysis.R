library(MCMCglmm)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(tibble)

# --------------------------
# Themes and variance function
# --------------------------
source("thesis_theme.R") # TPTN_theme() now includes cyanobacteria color scale
source("variance_explained.R") # ssbyvar()
source("config_paths.R") # data_folder path

# --------------------------
# Load data
# --------------------------
tn_dat <- read.csv(file.path(
  data_folder,
  "Total-Nitrogen",
  "TN_Master_Data.csv"
)) %>%
  mutate(actual_TN = TN_avg * 125) # account for dilution

metadata <- read.csv(file.path(
  data_folder,
  "Total-Nitrogen",
  "TNmetadata.csv"
)) %>%
  filter(sample %in% tn_dat$sample)

tn_dat <- tn_dat %>%
  left_join(metadata, by = "sample") %>%
  separate(
    treatment,
    into = c("geno", "cyano", "micro"),
    sep = "_",
    remove = FALSE
  ) %>%
  filter(
    !is.na(micro),
    actual_TN >= 0,
    sample != "36a2" # drop outlier
  ) %>%
  mutate(
    geno = factor(geno, levels = c("DR", "LR", "M", "TF", "UM", "W")),
    micro = factor(micro, levels = c("N", "H", "KF", "ODR"))
  )

geno_labels <- c(
  "DR" = "Durham \nReservoir",
  "LR" = "LaRoche \nPond",
  "M" = "Mill \nPond",
  "TF" = "Thompson \nFarm",
  "UM" = "Upper \nMill Pond",
  "W" = "Woodman \nRoad"
)
micro_labels <- c(
  "H" = "Home",
  "ODR" = "Dairy \nFarm",
  "N" = "None",
  "KF" = "Kingman \nFarm"
)

# --------------------------
# Normality check
# --------------------------
shapiro.test(tn_dat$actual_TN)

# --------------------------
# GLMM
# --------------------------
tn_mod1 <- MCMCglmm(
  actual_TN ~ geno:micro + micro:cyano + micro + cyano + geno,
  data = tn_dat,
  verbose = FALSE,
  nitt = 100000,
  thin = 500,
  burnin = 5000
)
summary(tn_mod1)

# --------------------------
# Plots
# --------------------------
# 1. Home microbiome by genotype
home_plot <- tn_dat %>%
  filter(micro == "H") %>%
  ggplot(aes(x = geno, y = actual_TN, color = cyano)) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 16,
    size = 2.5,
    position = position_dodge(0.6)
  ) +
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "errorbar",
    width = 0.35,
    size = 0.5,
    position = position_dodge(0.6)
  ) +
  scale_x_discrete(labels = geno_labels) +
  labs(
    title = "Home Microbiome by Genotype",
    x = "Duckweed Genotype",
    y = "Total Nitrogen (µg/L)",
    color = "Cyanobacteria \nTreatment"
  ) +
  TPTN_theme() +
  theme(legend.position = "none")
home_plot

# 2. Other microbiomes
othermicro_plot <- tn_dat %>%
  filter(micro != "H") %>%
  ggplot(aes(x = micro, y = actual_TN, color = cyano)) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 16,
    size = 2.5,
    position = position_dodge(0.6)
  ) +
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "errorbar",
    width = 0.35,
    size = 0.5,
    position = position_dodge(0.6)
  ) +
  scale_x_discrete(labels = micro_labels) +
  labs(
    title = "Other Microbiomes, Genotypes Combined",
    x = "Microbiome Treatment",
    y = "",
    color = "Cyanobacteria \nTreatment"
  ) +
  TPTN_theme() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
othermicro_plot

# --------------------------
# Variance explained
# --------------------------
tn_variance_data <- tibble(
  Factor = c("Cyanobacteria", "Genotype", "Microbiome"),
  Variance = c(
    ssbyvar(tn_dat$actual_TN, tn_dat$cyano),
    ssbyvar(tn_dat$actual_TN, tn_dat$geno),
    ssbyvar(tn_dat$actual_TN, tn_dat$micro)
  )
)

tn_variance_plot <- ggplot(
  tn_variance_data,
  aes(x = "", y = Variance, fill = Factor)
) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Variation Explained") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  TPvariance_theme()
tn_variance_plot

# --------------------------
# Combine plots
# --------------------------
TN_combined_plot <- plot_grid(
  home_plot,
  othermicro_plot,
  tn_variance_plot,
  nrow = 1,
  rel_widths = c(1.5, 1, 0.65)
)

TN_combined_plot

ggsave(
  file.path("plots", "Expt1_TN_combined_plot.jpg"),
  TN_combined_plot,
  width = 8,
  height = 3.25,
  dpi = 500
)
