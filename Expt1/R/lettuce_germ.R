source("config_paths.R")
source("globals.R")
source("thesis_theme.R")
source(file.path(path, "R/variance_explained.R"))

# load and reformat
lettuce_score <- read.csv(file.path(
    data,
    "lettuce_scores.csv"
)) |>
    select(!c(X12.Jul.24, X15.Jul.24, X15JulQualitativenotes, X)) |>
    mutate(
        # setting the germination to be a score, 0-3
        jul15NumericGerm = factor(
            jul15NumericGerm,
            levels = c("0", "1", "2", "3")
        ),

        # setting the health to be a score, 0-6
        Jul15HealthScore = factor(
            Jul15HealthScore,
            levels = c("0", "1", "2", "3", "4", "5", "6")
        )
    )

# setting scores as numeric
lettuce_score$J15numerichealth <- as.numeric(as.character(
    lettuce_score$Jul15HealthScore
))

# filtering out the control
lettuce_score_NC <- lettuce_score |>
    filter(Genotype != "control") |>
    drop_na(jul15NumericGerm, Jul15HealthScore)


# check normality
# W = 0.9212, p-value = 2.578e-13
shapiro.test(lettuce_score_NC$J15numerichealth)

# linear modeling
mod <- MCMCglmm(
    J15numerichealth ~ -1 + Microbe:Cyano, # removing intercept to show absolute means
    data = lettuce_score_NC,
    verbose = F,
    nitt = 101000,
    thin = 10,
    burnin = 1000
)
summary(mod)

label_positions <- lettuce_score_NC |>
    group_by(Cyano) |>
    summarise(max_score = max(J15numerichealth, na.rm = TRUE)) |>
    mutate(label = c("a", "a"))

lettuce_score_NC$Microbe <- factor(
    lettuce_score_NC$Microbe,
    levels = c("N", "H", "KF", "ODR")
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
    expt1_theme()

# using variance function
lettuce_variance_data <- tibble(
    Factor = c("Cyanobacteria", "Genotype", "Microbiome"),
    Variance = c(
        ssbyvar(lettuce_score$J15numerichealth, lettuce_score$Cyano),
        ssbyvar(lettuce_score$J15numerichealth, lettuce_score$Genotype),
        ssbyvar(lettuce_score$J15numerichealth, lettuce_score$Microbe)
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
    variance_theme()

lettuce_combined_plot <- plot_grid(
    lettuce_plot,
    lettuce_variance_plot,
    nrow = 1,
    rel_widths = c(5, 2)
)

ggsave(
    file.path("plots", "Expt1_lettuce_germ.jpg"),
    lettuce_combined_plot,
    width = 8,
    height = 3.25,
    dpi = 500
)
