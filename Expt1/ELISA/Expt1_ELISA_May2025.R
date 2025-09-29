library(MCMCglmm)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(ggtext)

# --------------------------
# Themes and variance function
# --------------------------
source("thesis_theme.R") # TPTN_theme() now includes cyanobacteria colors
source("variance_explained.R") # ssbyvar()
source("config_paths.R") # data_folder path
plots_folder <- "plots"

# --------------------------
# Load data
# --------------------------
ELISA_dat <- read.csv(file.path(
    data_folder,
    "ELISA",
    "Expt1_Phase2_ELISA.csv"
)) %>%
    select(-c(X, X.1, X.2, empty.tube.weight)) %>%
    separate(
        treatment,
        into = c("geno", "cyano", "micro"),
        sep = "_",
        remove = FALSE
    )

# Remove outliers
ELISA_dat_filtered <- ELISA_dat %>%
    filter(!is.na(cyano)) %>%
    filter(
        !(sample == "25_c2" & treatment == "M_N_N"),
        !(sample == "26_b3" & treatment == "W_Y_N")
    )

# --------------------------
# GLMM
# --------------------------
ELISA_glmm <- MCMCglmm(
    MC ~ -1 + cyano:micro,
    data = ELISA_dat_filtered,
    verbose = FALSE,
    nitt = 101000,
    thin = 10,
    burnin = 1000
)
summary(ELISA_glmm)

# Posterior means and 95% CI
post_summ <- summary(ELISA_glmm)$solutions
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
# Plot MC by microbiome source
# --------------------------
ELISA_dat_filtered <- ELISA_dat_filtered %>%
    mutate(micro = factor(micro, levels = c("N", "H", "KF", "ODR")))

MCplot <- ggplot(ELISA_dat_filtered, aes(x = micro, y = MC, color = cyano)) +
    stat_summary(
        fun = mean,
        geom = "point",
        size = 4,
        shape = 16,
        position = position_dodge(0.5)
    ) +
    stat_summary(
        fun.data = mean_sdl,
        fun.args = list(mult = 1),
        geom = "errorbar",
        width = 0.2,
        size = 1,
        position = position_dodge(0.5)
    ) +
    scale_x_discrete(
        labels = c(
            "N" = "Uninoculated",
            "H" = "Home",
            "KF" = "Kingman Farm",
            "ODR" = "Dairy Farm"
        )
    ) +
    labs(
        x = "Microbiome Source",
        y = "Microcystin (µg/g Duckweed)",
        color = "Cyanobacteria Treatment"
    ) +
    TPTN_theme() +
    theme(legend.position = "none")

MCplot

# --------------------------
# Variance explained
# --------------------------
MC_variance_data <- tibble(
    Factor = c("Cyanobacteria", "Microbiome", "Genotype"),
    Variance = c(
        ssbyvar(ELISA_dat_filtered$MC, ELISA_dat_filtered$cyano),
        ssbyvar(ELISA_dat_filtered$MC, ELISA_dat_filtered$micro),
        ssbyvar(ELISA_dat_filtered$MC, ELISA_dat_filtered$geno)
    )
)

MC_variance_plot <- ggplot(
    MC_variance_data,
    aes(x = "", y = Variance, fill = Factor)
) +
    geom_bar(stat = "identity", position = "stack") +
    labs(y = "Variation Explained") +
    scale_y_continuous(
        limits = c(0, 1),
        expand = expansion(mult = c(0, 0.05))
    ) +
    TPvariance_theme()

MC_variance_plot

# --------------------------
# Combined plot
# --------------------------
MC_combined_plot <- plot_grid(
    MCplot,
    MC_variance_plot,
    nrow = 1,
    rel_widths = c(5, 2)
)

MC_combined_plot

ggsave(
    file.path(plots_folder, "Expt1_MC_plot.jpg"),
    plot = MC_combined_plot,
    width = 6,
    height = 3.5,
    dpi = 500
)

# --------------------------
# Line plot by genotype
# --------------------------
MC_line_plot <- ELISA_dat %>%
    drop_na() %>%
    group_by(cyano, geno) %>%
    summarise(mean_conc = mean(MC), sd_conc = sd(MC), .groups = 'drop') %>%
    mutate(se_conc = sd_conc / sqrt(n())) %>%
    ggplot(aes(
        x = cyano,
        y = mean_conc,
        fill = geno,
        color = geno,
        group = geno
    )) +
    geom_path(linewidth = 1, lineend = "round") +
    geom_errorbar(
        aes(ymin = mean_conc - se_conc, ymax = mean_conc + se_conc),
        width = 0.05,
        size = 0.75
    ) +
    scale_x_discrete(expand = c(0, 0.1), labels = c("N" = "No", "Y" = "Yes")) +
    labs(
        x = "*M. aeruginosa* Spike",
        y = "Average Microcystin (µg/g Duckweed)",
        color = "Duckweed\nGenotype"
    ) +
    scale_color_viridis(discrete = TRUE) +
    theme(
        legend.position = "right",
        axis.title.x = ggtext::element_markdown(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)
    )

MC_line_plot

ggsave(
    file.path(plots_folder, "Expt1_MC_line_plot.jpg"),
    plot = MC_line_plot,
    width = 4,
    height = 3,
    dpi = 300
)

# --------------------------
# Split by microbiome
# --------------------------
micro_mod <- MCMCglmm(
    MC ~ -1 + micro:cyano,
    data = ELISA_dat_filtered,
    verbose = FALSE,
    nitt = 101000,
    thin = 10,
    burnin = 1000
)
summary(micro_mod)

micro_labels <- c(
    "H" = "Home",
    "ODR" = "Dairy \nFarm",
    "N" = "None",
    "KF" = "Kingman \nFarm"
)

MCplot_micro <- ggplot(
    ELISA_dat_filtered,
    aes(x = micro, y = MC, color = cyano)
) +
    stat_summary(
        fun = mean,
        geom = "point",
        size = 4,
        shape = 16,
        position = position_dodge(0.6)
    ) +
    stat_summary(
        fun.data = mean_sdl,
        fun.args = list(mult = 1),
        geom = "errorbar",
        width = 0.3,
        size = 1,
        position = position_dodge(0.6)
    ) +
    scale_color_manual(labels = c("N" = "No", "Y" = "Yes")) +
    scale_x_discrete(labels = micro_labels) +
    labs(
        x = "Microbiome Source",
        y = "Microcystin (µg/g Duckweed)",
        color = "*M. aeruginosa* <br> Spike"
    ) +
    TPTN_theme() +
    theme(
        legend.position = "right",
        legend.title = ggtext::element_markdown(size = 9),
        legend.text = element_text(size = 8)
    )

MCplot_micro

#ggsave(
#    file.path(plots_folder, "Expt1_MC_micro_plot.jpg"),
#    plot = MCplot_micro,
#    width = 5,
#    height = 3,
#   dpi = 500
#)
