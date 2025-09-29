library(MCMCglmm)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsignif)
library(cowplot)
library(tibble)
library(ggtext)

TP_dat <- read.csv("Expt1_TP.csv")

TP_dat <- TP_dat |>
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

#is the data normally distributed
shapiro.test(TP_dat$ppb)

#filter to include only home microbiome
tp_microH <- TP_dat |>
    filter(micro == "H")

mod_microH <- MCMCglmm(
    ppb ~ -1 + cyano:geno,
    data = tp_microH,
    verbose = F,
    nitt = 101000,
    thin = 10,
    burnin = 1000
)
summary(mod_microH)

#filter to include only other microbiomes
tpdat_noHmicro <- TP_dat |>
    filter(micro != "H") |>
    mutate(micro = factor(micro, levels = c("N", "KF", "ODR")))

mod_noHmicro <- MCMCglmm(
    ppb ~ -1 + cyano:micro,
    data = tpdat_noHmicro,
    verbose = F,
    nitt = 101000,
    thin = 10,
    burnin = 1000
)
summary(mod_noHmicro)

post_summ <- summary(mod_noHmicro)$solutions
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

# First plot (Home Microbiome by Genotype)
TPplot_microH <- tp_microH |>
    ggplot(aes(x = geno, y = ppb, color = cyano)) +
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
        width = .35,
        size = .5,
        position = position_dodge(width = 0.6)
    ) +
    scale_color_manual(values = c("N" = "black", "Y" = "aquamarine4")) +
    scale_x_discrete(labels = geno_labels) +
    ylim(0, 50000) +
    labs(
        title = "Home Microbiome by Genotype",
        x = "Duckweed Genotype",
        y = "Total Phosphorus (µg/L)"
    ) +
    theme_classic() +
    theme(
        legend.position = "none",
        axis.title.x = element_text(size = 7, margin = margin(t = 8)),
        axis.title.y = element_text(size = 9, margin = margin(r = 10)),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 7)
    ) +
    geom_hline(
        yintercept = 45570,
        linetype = "dashed",
        color = "red",
        size = 1,
        show.legend = TRUE
    ) +
    annotate(
        "text",
        x = 1,
        y = 45570 + 1500,
        label = "Initial TP (45,570 µg/L)",
        color = "red",
        size = 2,
        hjust = 0.2,
        vjust = 0
    )
TPplot_microH

# Second plot (Other Microbiomes, Genotypes Combined)
TPplot_others <- tpdat_noHmicro |>
    ggplot(aes(x = micro, y = ppb, color = cyano)) +
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
        size = .5,
        position = position_dodge(width = 0.6)
    ) +
    scale_color_manual(
        values = c("N" = "black", "Y" = "aquamarine4"),
        name = "Cyanobacteria"
    ) +
    scale_x_discrete(labels = micro_labels) +
    ylim(0, 50000) +
    labs(
        title = "Other Microbiomes, Genotypes Combined",
        x = "Microbiome",
        y = ""
    ) +
    theme_classic() +
    theme(
        legend.position = "none",
        axis.title.x = element_text(size = 7, margin = margin(t = 8)),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_blank(),
        axis.line.y = element_line(
            colour = "black",
            linetype = "dashed",
            size = .5
        ),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 7)
    ) +
    geom_hline(
        yintercept = 45570,
        linetype = "dashed",
        color = "red",
        size = 1,
        show.legend = TRUE
    )


ssbyvar <- function(response, category.vec) {
    #sums of squares function
    means <- tapply(response, category.vec, mean, na.rm = T) #take the means by category
    ssresid <- sum(sapply(sort(unique(category.vec)), function(z) {
        sum(
            (response[category.vec == z] - means[names(means) == z])^2,
            na.rm = T
        )
    })) #square of difference of each datapoint from its associated treatment mean (residual variation)
    sstot <- sum((response - mean(response, na.rm = T))^2, na.rm = T) #square of difference of each datapoint from the grand mean (total variation)
    sst <- (sstot - ssresid) # total variation - residual variation = treatment variation
    return(sst / sstot) # treatment variance as a fraction of total variation
}

TPvariance_cyano <- ssbyvar(TP_dat$ppb, TP_dat$cyano)
TPvariance_geno <- ssbyvar(TP_dat$ppb, TP_dat$geno)
TPvariance_micro <- ssbyvar(TP_dat$ppb, TP_dat$micro)

# Create a data frame with variance and error
TPvariance_data <- data.frame(
    Factor = c("Cyanobacteria", "Genotype", "Microbiome"),
    Variance = c(TPvariance_cyano, TPvariance_geno, TPvariance_micro)
)

TP_variance_plot <- ggplot(
    TPvariance_data,
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


TP_variance_plot

# Combine the two main plots into one grid
TPmicro_combined_plot <- plot_grid(
    TPplot_microH,
    TPplot_others,
    TP_variance_plot,
    nrow = 1,
    rel_widths = c(1.5, 1, .65)
)

TPmicro_combined_plot

ggsave(
    "Expt1_TPmicro_combined_plot.jpg",
    TPmicro_combined_plot,
    width = 8,
    height = 3.25,
    dpi = 500
)
