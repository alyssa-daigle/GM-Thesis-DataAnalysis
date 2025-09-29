library(MCMCglmm)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(viridis)

setwd(
    "~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/Expt1/expt1_dataanalysis/ELISA"
)

ELISA_dat <- read.csv("/Users/alyssadaigle/Desktop/csvs/Expt1_Phase2_ELISA.csv")

ELISA_dat <- ELISA_dat |> select(-c(X, X.1, X.2, empty.tube.weight))

ELISA_dat <- ELISA_dat |>
    separate(
        treatment,
        into = c("geno", "cyano", "micro"),
        sep = "_",
        remove = FALSE
    )

#removing outliers
ELISA_dat_filtered <- ELISA_dat %>%
    filter(!is.na(cyano)) %>%
    filter(
        !(sample == "25_c2" & treatment == "M_N_N") &
            !(sample == "26_b3" & treatment == "W_Y_N")
    )

#linear modeling
ELISA_glmm <- MCMCglmm(
    MC ~ -1 + cyano:micro,
    data = ELISA_dat_filtered,
    verbose = F,
    nitt = 101000,
    thin = 10,
    burnin = 1000
)
summary(ELISA_glmm)


post_summ <- summary(ELISA_glmm)$solutions
ci_df <- as.data.frame(post_summ)
ci_df$Effect <- rownames(ci_df)
names(ci_df)[c(1, 2, 3)] <- c("PostMean", "Lower95CI", "Upper95CI")

ggplot(ci_df, aes(x = Effect, y = PostMean)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = Lower95CI, ymax = Upper95CI), width = 0.2) +
    theme_minimal() +
    labs(
        title = "Posterior Means and 95% Credible Intervals",
        y = "Posterior Mean ± 95% CI",
        x = "Effect"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


#plotting cyanoY vs cyanoN for MC concentration
ELISA_dat_filtered$micro <- factor(
    ELISA_dat_filtered$micro,
    levels = c("N", "H", "KF", "ODR")
)

MCplot <- ggplot(ELISA_dat_filtered, aes(x = micro, y = MC, color = cyano)) +
    stat_summary(
        fun = mean,
        geom = "point",
        shape = 16,
        size = 4,
        position = position_dodge(width = .5)
    ) +
    stat_summary(
        fun.data = mean_sdl,
        fun.args = list(mult = 1),
        geom = "errorbar",
        width = 0.2,
        size = 1,
        position = position_dodge(width = .5)
    ) +
    scale_color_manual(values = c("N" = "black", "Y" = "aquamarine4")) +
    scale_x_discrete(
        labels = c(
            "N" = "Uninoculated",
            "H" = "Home",
            "KF" = "Kingman Farm",
            "ODR" = "Dairy Farm"
        )
    ) +
    labs(x = "Microbiome Source", y = "Microcystin (µg/g Duckweed)") +
    theme_cowplot() +
    theme(
        legend.position = "none",
        axis.title.x = ggtext::element_markdown(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)
    )

MCplot

#sums of squares function
ssbyvar <- function(response, category.vec) {
    means <- tapply(response, category.vec, mean, na.rm = T)
    ssresid <- sum(sapply(sort(unique(category.vec)), function(z) {
        sum(
            (response[category.vec == z] - means[names(means) == z])^2,
            na.rm = T
        )
    }))
    sstot <- sum((response - mean(response, na.rm = T))^2, na.rm = T)
    sst <- (sstot - ssresid)
    return(sst / sstot)
}

variance_cyano <- ssbyvar(ELISA_dat_filtered$MC, ELISA_dat_filtered$cyano)
variance_micro <- ssbyvar(ELISA_dat_filtered$MC, ELISA_dat_filtered$micro)
variance_geno <- ssbyvar(ELISA_dat_filtered$MC, ELISA_dat_filtered$geno)

MC_variance_data <- data.frame(
    Factor = c("Cyanobacteria", "Microbiome", "Genotype"),
    Variance = c(variance_cyano, variance_micro, variance_geno)
)

#stacked bar plot
MC_variance_plot <- ggplot(
    MC_variance_data,
    aes(x = "", y = Variance, fill = Factor)
) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "", y = "Variation Explained", fill = "Effect") +
    scale_fill_manual(
        values = c(
            "Cyanobacteria" = "aquamarine4",
            "Genotype" = "darkblue",
            "Microbiome" = "wheat"
        )
    ) +
    guides(
        fill = guide_legend(
            title.position = "top",
            title.hjust = .5,
            override.aes = list(size = 1)
        )
    ) +
    scale_y_continuous(
        limits = c(0, 1), # Set the limits between 0 and 1
        expand = expansion(mult = c(0, 0.05))
    ) +
    theme_cowplot() +
    theme(
        legend.position = "right",
        legend.justification = "center",
        legend.box.just = "center",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, "cm"),
        legend.box = "vertical",
        legend.box.spacing = unit(0, "cm"),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.ticks.x = element_blank(),
    )

MC_variance_plot

MC_combined_plot <- plot_grid(
    MCplot,
    MC_variance_plot,
    nrow = 1,
    rel_widths = c(5, 2)
)

MC_combined_plot

ggsave2(
    "Expt1_MC_plot.jpg",
    plot = MC_combined_plot,
    width = 6,
    height = 3.5,
    dpi = 500
)
ggsave(
    "~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/WRITING/plots/Expt1_MC_plot.jpg",
    MC_combined_plot,
    width = 6,
    height = 3.5,
    dpi = 500
)


#Line plot, MC by genotype
MC_line_plot <- ELISA_dat |>
    drop_na() |>
    group_by(cyano, geno) |>
    summarise(
        mean_conc = mean(MC),
        sd_conc = sd(MC),
        .groups = 'drop'
    ) |>
    mutate(
        se_conc = sd_conc / sqrt(n()) # Calculate standard error
    ) |>
    ggplot(aes(
        x = cyano,
        y = mean_conc,
        fill = geno,
        color = geno,
        group = geno
    )) +
    scale_x_discrete(expand = c(0, 0.1), labels = c("N" = "No", "Y" = "Yes")) +
    geom_path(linewidth = 1, lineend = "round") +
    geom_errorbar(
        aes(ymin = mean_conc - se_conc, ymax = mean_conc + se_conc),
        width = 0.05,
        size = .75
    ) +
    labs(
        x = "*M. aeruginosa* Spike",
        y = "Average Microcystin (µg/g Duckweed)",
        color = "Duckweed\nGenotype"
    ) +
    guides(color = guide_legend(ncol = 1)) +
    theme_cowplot() +
    theme(
        legend.position = "right",
        axis.title.x = ggtext::element_markdown(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)
    ) +
    scale_color_viridis(discrete = TRUE)

# Display the line plot
MC_line_plot

ggsave(
    "Expt1_MC_line_plot.jpg",
    plot = MC_line_plot,
    width = 4,
    height = 3,
    dpi = 300
)


##splitting it up by micro

micro_mod <- MCMCglmm(
    MC ~ -1 + micro:cyano,
    data = ELISA_dat_filtered,
    verbose = F,
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

ELISA_dat_filtered <- ELISA_dat_filtered |>
    mutate(micro = factor(micro, levels = c("N", "H", "KF", "ODR")))

MCplot_micro <- ggplot(
    ELISA_dat_filtered,
    aes(x = micro, y = MC, color = cyano)
) +
    #geom_jitter(shape = 1, width = 0.2, size = 2) +
    stat_summary(
        fun = mean,
        geom = "point",
        shape = 16,
        size = 4,
        position = position_dodge(width = 0.6)
    ) +
    stat_summary(
        fun.data = mean_sdl,
        fun.args = list(mult = 1),
        geom = "errorbar",
        width = 0.3,
        size = 1,
        position = position_dodge(width = 0.6)
    ) +
    scale_color_manual(
        values = c("N" = "black", "Y" = "aquamarine4"),
        labels = c("N" = "No", "Y" = "Yes")
    ) +
    labs(
        x = "Microbiome Source",
        color = "*M. aeruginosa* <br> Spike",
        y = "Microcystin (µg/g Duckweed)"
    ) +
    scale_x_discrete(labels = micro_labels) +
    theme_cowplot() +
    theme(
        legend.position = "right",
        legend.title = ggtext::element_markdown(size = 9),
        legend.text = element_text(size = 8),
        axis.title.y = element_text(size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 9)
    )


MCplot_micro

ggsave(
    "Expt1_MC_micro_plot.jpg",
    plot = MCplot_micro,
    width = 5,
    height = 3,
    dpi = 500
)
