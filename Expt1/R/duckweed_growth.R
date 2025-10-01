source("config_paths.R")
source("globals.R")
source("thesis_theme.R")
source(file.path(path, "R/variance_explained.R"))

# load data
growthdata <- read.csv(file.path(
    data,
    "frondarea.csv"
)) |>
    mutate(RGR = as.numeric(RGR)) |>
    separate(
        treatment,
        into = c("geno", "cyano", "micro"),
        sep = "_",
        remove = FALSE
    )

# checking normality
shapiro.test(growthdata$RGR)

# ---------------------------------------------------------------------------------------------
# linear models
mod <- MCMCglmm(
    RGR ~ -1 + cyano:micro,
    data = growthdata,
    verbose = F,
    nitt = 11000,
    thin = 10,
    burnin = 100
)
summary(mod)

# Posterior summary plot
# shows credible intervals for each interaction
# no CI overlap -> significantly different
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
growthdata$micro <- factor(growthdata$micro, levels = c("N", "H", "KF", "ODR"))

# Create plot
growth_plot <- ggplot(
    data = growthdata,
    aes(x = micro, y = RGR, color = cyano)
) +
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
    scale_x_discrete(labels = micro_labels) +
    labs(x = "Microbiome Source", y = "Duckweed Relative Growth Rate") +
    theme_cowplot() +
    expt1_theme() +
    theme(legend.position = "none")
growth_plot

# applying variance explained function
growth_variance_data <- tibble(
    Factor = c("Cyanobacteria", "Genotype", "Microbiome"),
    Variance = c(
        ssbyvar(growthdata$RGR, growthdata$cyano),
        ssbyvar(growthdata$RGR, growthdata$cyano),
        ssbyvar(growthdata$RGR, growthdata$cyano)
    )
)


# variance explained plot
growth_variance_plot <- ggplot(
    growth_variance_data,
    aes(x = "", y = Variance, fill = Factor)
) +
    geom_bar(stat = "identity", position = "stack") +
    labs(y = "Variation Explained") +
    scale_y_continuous(
        limits = c(0, 1),
        expand = expansion(mult = c(0, 0.05))
    ) +
    variance_theme()
growth_variance_plot

growth_combined_plot <- plot_grid(
    growth_plot,
    growth_variance_plot,
    nrow = 1,
    rel_widths = c(5, 2)
)

growth_combined_plot

ggsave(
    file.path("plots", "Expt1_growth_plot.jpg"),
    growth_combined_plot,
    width = 8,
    height = 3.25,
    dpi = 500
)
