source("config_paths.R")
source("globals.R")
source("thesis_theme.R")
source(file.path(path, "R/variance_explained.R"))

# load data
tndata <- read.csv(file.path(
  data,
  "total_nitrogen.csv"
)) |>
  mutate(actual_TN = TN_avg * 125) # account for dilution

# load metadata
metadata <- read.csv(file.path(
  data,
  "total_nitrogen_meta.csv"
)) |>
  filter(sample %in% tndata$sample)

tndata <- tndata |>
  left_join(metadata, by = "sample") |>
  separate(
    treatment,
    into = c("geno", "cyano", "micro"),
    sep = "_",
    remove = FALSE
  ) |>
  filter(
    !is.na(micro),
    actual_TN >= 0,
    sample != "36a2" # drop outlier
  ) |>
  mutate(
    geno = factor(geno, levels = c("DR", "LR", "M", "TF", "UM", "W")),
    micro = factor(micro, levels = c("N", "H", "KF", "ODR"))
  )

# checking normality
shapiro.test(tndata$actual_TN)

# ---------------------------------------------------------------------------------------------
# linear models

# first, see which effect contributes to most variance
mod <- MCMCglmm(
  actual_TN ~ micro + cyano + geno,
  data = tndata,
  verbose = FALSE,
  nitt = 100000,
  thin = 500,
  burnin = 5000
)
summary(mod)

# cyano and micro seem to have a greater affect than geno
# diving into micro effects more
microH <- filter(tndata, micro == "H")
mod_microH <- MCMCglmm(
  actual_TN ~ -1 + cyano:geno,
  data = microH, # removing intercept to show absolute means
  verbose = FALSE,
  nitt = 101000,
  thin = 10,
  burnin = 1000
)
summary(mod_microH)

# other microbiomes, cyano Y vs N
noHmicro <- tndata |>
  filter(micro != "H") |>
  mutate(micro = factor(micro, levels = c("N", "KF", "ODR")))
mod_noHmicro <- MCMCglmm(
  actual_TN ~ -1 + cyano:micro,
  data = noHmicro,
  verbose = FALSE,
  nitt = 101000,
  thin = 10,
  burnin = 1000
)
summary(mod_noHmicro)

# Posterior summary plot
# shows credible intervals for each interaction
# no CI overlap -> significantly different
post_summ <- summary(mod_microH)$solutions # replace with mod_microH/mod_noHmicro to compare
ci_df <- as.data.frame(post_summ)
ci_df$Effect <- rownames(ci_df)
names(ci_df)[1:3] <- c("PostMean", "Lower95CI", "Upper95CI")

ggplot(ci_df, aes(x = Effect, y = PostMean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower95CI, ymax = Upper95CI), width = 0.2) +
  labs(
    title = "Posterior Means and 95% Credible Intervals",
    y = "Posterior Mean ± 95% CI",
    x = "Effect"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---------------------------------------------------------------------------------------------
# Plots
microH$cyano <- factor(microH$cyano, levels = c("N", "Y"))
noHmicro$cyano <- factor(noHmicro$cyano, levels = c("N", "Y"))

# 1. Home microbiome alone
home_plot <- tndata |>
  filter(micro == "H") |>
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
  scale_color_manual(
    values = c("N" = "black", "Y" = "aquamarine4"),
  ) +
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
othermicro_plot <- tndata |>
  filter(micro != "H") |>
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
  scale_color_manual(
    values = c("N" = "black", "Y" = "aquamarine4"),
  ) +
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

# applying variance explained function
tn_variance_data <- tibble(
  Factor = c("Cyanobacteria", "Genotype", "Microbiome"),
  Variance = c(
    ssbyvar(tndata$actual_TN, tndata$cyano),
    ssbyvar(tndata$actual_TN, tndata$geno),
    ssbyvar(tndata$actual_TN, tndata$micro)
  )
)

# variance explained plot
tn_variance_plot <- ggplot(
  tn_variance_data,
  aes(x = "", y = Variance, fill = Factor)
) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Variation Explained") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  TPvariance_theme()
tn_variance_plot

# ---------------------------------------------------------------------------------------------
# Combine and save
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
