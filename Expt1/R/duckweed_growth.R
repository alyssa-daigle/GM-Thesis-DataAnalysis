library(MCMCglmm)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)

setwd(
    "~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/Expt1/expt1_dataanalysis/Growth"
)

dat <- read.csv("~/Desktop/csvs/Expt1_Batch2_frondarea.csv")

dat <- dat %>%
    mutate(RGR = as.numeric(RGR))

dat <- dat |>
    separate(
        treatment,
        into = c("geno", "cyano", "micro"),
        sep = "_",
        remove = FALSE
    )

#both normally distributed
shapiro.test(dat$RGR)

mod <- MCMCglmm(
    RGR ~ -1 + cyano:micro,
    data = dat,
    verbose = F,
    nitt = 11000,
    thin = 10,
    burnin = 100
)
summary(mod)

post_summ <- summary(mod)$solutions
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

#order micro by N, H, KF, then ODR
dat$micro <- factor(dat$micro, levels = c("N", "H", "KF", "ODR"))

# Create plot
growth_plot <- ggplot(data = dat, aes(x = micro, y = RGR, color = cyano)) +
    #geom_jitter(shape = 1, width = .3, size = 2,) +
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
    labs(x = "Microbiome Source", y = "Duckweed Relative Growth Rate") +
    theme_cowplot() +
    theme(
        legend.position = "none",
        axis.title.x = ggtext::element_markdown(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)
    )

growth_plot

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

growth_variance_cyano <- ssbyvar(dat$RGR, dat$cyano)
growth_variance_geno <- ssbyvar(dat$RGR, dat$geno)
growth_variance_micro <- ssbyvar(dat$RGR, dat$micro)

# Create a data frame with variance and error
growth_variance_data <- data.frame(
    Factor = c("Cyanobacteria", "Genotype", "Microbiome"),
    Variance = c(
        growth_variance_cyano,
        growth_variance_geno,
        growth_variance_micro
    )
)

# variation explained plot
growth_variance_plot <- ggplot(
    growth_variance_data,
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

growth_variance_plot

growth_combined_plot <- plot_grid(
    growth_plot,
    growth_variance_plot,
    nrow = 1,
    rel_widths = c(5, 2)
)

growth_combined_plot

ggsave(
    "Expt1_growth_plot.jpg",
    growth_combined_plot,
    width = 6,
    height = 3.5,
    dpi = 500
)
ggsave(
    "~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/WRITING/plots/Expt1_growth_plot.jpg",
    growth_combined_plot,
    width = 6,
    height = 3.5,
    dpi = 500
)
