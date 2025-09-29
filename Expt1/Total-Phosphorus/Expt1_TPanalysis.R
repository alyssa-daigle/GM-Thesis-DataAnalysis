library(MCMCglmm)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(tibble)

# --------------------------
# themes and variance function
# --------------------------
source("thesis_theme.R")
source("variance_explained.R")
source("config_paths.R")

# --------------------------
# data
# --------------------------
TP_dat <- read.csv(file.path(
    data_folder,
    "Total-Phosphorus",
    "Expt1_TP.csv"
)) |>
    separate(
        treatment,
        into = c("geno", "cyano", "micro"),
        sep = "_",
        remove = FALSE
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
# normality check
# --------------------------
shapiro.test(TP_dat$ppb)

# --------------------------
# GLMMs
# --------------------------
# Home microbiome
tp_microH <- filter(TP_dat, micro == "H")
mod_microH <- MCMCglmm(
    ppb ~ -1 + cyano:geno,
    data = tp_microH,
    verbose = FALSE,
    nitt = 101000,
    thin = 10,
    burnin = 1000
)
summary(mod_microH)

# Other microbiomes
tpdat_noHmicro <- TP_dat |>
    filter(micro != "H") |>
    mutate(micro = factor(micro, levels = c("N", "KF", "ODR")))

mod_noHmicro <- MCMCglmm(
    ppb ~ -1 + cyano:micro,
    data = tpdat_noHmicro,
    verbose = FALSE,
    nitt = 101000,
    thin = 10,
    burnin = 1000
)
summary(mod_noHmicro)

# Posterior summary plot
post_summ <- summary(mod_noHmicro)$solutions
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
    TPTN_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --------------------------
# Plots
# --------------------------
# 1. Home Microbiome
TPplot_microH <- ggplot(tp_microH, aes(x = geno, y = ppb, color = cyano)) +
    stat_summary(
        fun = mean,
        geom = "point",
        shape = 16,
        size = 2.5,
        position = position_dodge(width = 0.6)
    ) +
    stat_summary(
        fun.data = mean_sdl,
        fun.args = list(mult = 1),
        geom = "errorbar",
        width = 0.35,
        size = 0.5,
        position = position_dodge(width = 0.6)
    ) +
    scale_x_discrete(labels = geno_labels) +
    ylim(0, 50000) +
    labs(
        title = "Home Microbiome by Genotype",
        x = "Duckweed Genotype",
        y = "Total Phosphorus (µg/L)"
    ) +
    TPTN_theme() +
    theme(legend.position = "none") +
    geom_hline(
        yintercept = 45570,
        linetype = "dashed",
        color = "red",
        size = 1
    ) +
    annotate(
        "text",
        x = 1,
        y = 45570 + 1500,
        label = "Initial TP (45,570 µg/L)",
        color = "red",
        size = 2,
        hjust = 0.2
    )

# 2. Other Microbiomes
TPplot_others <- ggplot(
    tpdat_noHmicro,
    aes(x = micro, y = ppb, color = cyano)
) +
    stat_summary(
        fun = mean,
        geom = "point",
        shape = 16,
        size = 2.5,
        position = position_dodge(width = 0.6)
    ) +
    stat_summary(
        fun.data = mean_sdl,
        fun.args = list(mult = 1),
        geom = "errorbar",
        width = 0.3,
        size = 0.5,
        position = position_dodge(width = 0.6)
    ) +
    scale_x_discrete(labels = micro_labels) +
    ylim(0, 50000) +
    labs(
        title = "Other Microbiomes, Genotypes Combined",
        x = "Microbiome",
        y = ""
    ) +
    TPTN_theme() +
    theme(
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
    ) +
    geom_hline(yintercept = 45570, linetype = "dashed", color = "red", size = 1)

# --------------------------
# Variance explained
# --------------------------
TPvariance_data <- tibble(
    Factor = c("Cyanobacteria", "Genotype", "Microbiome"),
    Variance = c(
        ssbyvar(TP_dat$ppb, TP_dat$cyano),
        ssbyvar(TP_dat$ppb, TP_dat$geno),
        ssbyvar(TP_dat$ppb, TP_dat$micro)
    )
)

TP_variance_plot <- ggplot(
    TPvariance_data,
    aes(x = "", y = Variance, fill = Factor)
) +
    geom_bar(stat = "identity", position = "stack") +
    labs(y = "Variation Explained") +
    scale_y_continuous(
        limits = c(0, 1),
        expand = expansion(mult = c(0, 0.05))
    ) +
    TPvariance_theme()

# --------------------------
# Combine and save
# --------------------------
TPmicro_combined_plot <- plot_grid(
    TPplot_microH,
    TPplot_others,
    TP_variance_plot,
    nrow = 1,
    rel_widths = c(1.5, 1, 0.65)
)

TPmicro_combined_plot

ggsave(
    file.path("plots", "Expt1_TPmicro_combined_plot.jpg"),
    TPmicro_combined_plot,
    width = 8,
    height = 3.25,
    dpi = 500
)
