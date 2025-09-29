library(tidyverse)
library(cowplot)
library(ggpubr)
library(MCMCglmm)

setwd(
    "~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/Expt1/expt1_dataanalysis/lettuce"
)

lettuce_score <- read.csv("/Users/alyssadaigle/Desktop/csvs/lettuce_scores.csv")

lettuce_score <- lettuce_score |>
    select(!c(X12.Jul.24, X15.Jul.24, X15JulQualitativenotes, X))

lettuce_score <- lettuce_score %>%
    mutate(
        jul15NumericGerm = factor(
            jul15NumericGerm,
            levels = c("0", "1", "2", "3")
        ),
        Jul15HealthScore = factor(
            Jul15HealthScore,
            levels = c("0", "1", "2", "3", "4", "5", "6")
        )
    )

lettuce_score$J15numerichealth <- as.numeric(as.character(
    lettuce_score$Jul15HealthScore
))

lettuce_score_NC <- lettuce_score %>%
    filter(Genotype != "control") %>%
    drop_na(jul15NumericGerm, Jul15HealthScore)

shapiro.test(lettuce_score_NC$J15numerichealth)
mod <- MCMCglmm(
    J15numerichealth ~ -1 + Microbe:Cyano,
    data = lettuce_score_NC,
    verbose = F,
    nitt = 101000,
    thin = 10,
    burnin = 1000
)
summary(mod)

label_positions <- lettuce_score_NC %>%
    group_by(Cyano) %>%
    summarise(max_score = max(J15numerichealth, na.rm = TRUE)) %>%
    mutate(label = c("a", "a"))

lettuce_score_NC$Microbe <- factor(
    lettuce_score_NC$Microbe,
    levels = c("N", "H", "KF", "ODR")
)

micro_labels <- c(
    "H" = "Home",
    "ODR" = "Dairy Farm",
    "N" = "None",
    "KF" = "Kingman Farm"
)

lettuce_plot <- lettuce_score_NC |>
    ggplot(aes(x = Microbe, y = J15numerichealth, color = Cyano)) +
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
        width = 0.2,
        size = 1,
        position = position_dodge(width = 0.6)
    ) +
    scale_color_manual(values = c("N" = "black", "Y" = "aquamarine4")) +
    scale_x_discrete(labels = micro_labels) +
    scale_y_continuous(limits = c(0, 7), breaks = seq(0, 6, by = 1)) +
    labs(x = "Microbiome Source", y = "Lettuce Germinant Health Score") +
    theme_cowplot() +
    theme(
        legend.position = "none",
        axis.title.x = ggtext::element_markdown(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)
    )

lettuce_plot

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

lettuce_variance_cyano <- ssbyvar(
    lettuce_score$J15numerichealth,
    lettuce_score$Cyano
)
lettuce_variance_geno <- ssbyvar(
    lettuce_score$J15numerichealth,
    lettuce_score$Genotype
)
lettuce_variance_micro <- ssbyvar(
    lettuce_score$J15numerichealth,
    lettuce_score$Microbe
)

# Create a data frame with variance and error
lettuce_variance_data <- data.frame(
    Factor = c("Cyanobacteria", "Genotype", "Microbiome"),
    Variance = c(
        lettuce_variance_cyano,
        lettuce_variance_geno,
        lettuce_variance_micro
    )
)

#variation explained plot
lettuce_variance_plot <- ggplot(
    lettuce_variance_data,
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


lettuce_variance_plot

lettuce_combined_plot <- plot_grid(
    lettuce_plot,
    lettuce_variance_plot,
    nrow = 1,
    rel_widths = c(5, 2)
)

lettuce_combined_plot

ggsave2(
    "Expt1_lettuce_combined_plot.jpg",
    plot = lettuce_combined_plot,
    width = 6,
    height = 3.5,
    dpi = 500
)
ggsave(
    "~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/WRITING/plots/Expt1_lettuce_combined_plot.jpg",
    lettuce_combined_plot,
    width = 6,
    height = 3.5,
    dpi = 500
)
