# theme_duckweed.R
library(ggplot2)
library(grid) # for unit()

TPTN_theme <- function(base_size = 7, base_family = "") {
  list(
    # Color scale for cyanobacteria treatment
    scale_color_manual(values = c("N" = "black", "Y" = "aquamarine4")),

    # Theme
    theme_classic(base_size = base_size, base_family = base_family) %+replace%
      theme(
        # Titles
        plot.title = element_text(hjust = 0.5, face = "bold", size = base_size),

        # Axis
        axis.title.x = element_text(size = base_size, margin = margin(t = 8)),
        axis.title.y = element_text(
          size = base_size,
          margin = margin(r = 10),
          angle = 90
        ),
        axis.text.x = element_text(size = base_size - 1),
        axis.text.y = element_text(size = base_size - 1),
        axis.ticks.x = element_line(size = 0.4),
        axis.ticks.y = element_line(size = 0.4),

        # Legend
        legend.position = "right",
        legend.justification = "center",
        legend.title = element_text(size = base_size),
        legend.text = element_text(size = base_size - 1),
        legend.key.size = unit(0.25, "cm"),
        legend.box = "vertical",
        legend.box.just = "center",

        # Panel/Background
        panel.grid = element_blank(),
        panel.border = element_blank()
      )
  )
}

TPvariance_theme <- function(base_size = 7, base_family = "") {
  list(
    # Fill scale
    scale_fill_manual(
      values = c(
        "Cyanobacteria" = "aquamarine4",
        "Genotype" = "darkblue",
        "Microbiome" = "wheat"
      )
    ),

    # Legend
    guides(
      fill = guide_legend(
        title = "Effect",
        title.position = "top",
        title.hjust = .5,
        override.aes = list(size = 1)
      )
    ),

    # Theme
    theme_classic(base_size = base_size, base_family = base_family) %+replace%
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_text(
          size = base_size,
          margin = margin(r = 8),
          angle = 90
        ),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = base_size - 1),
        axis.ticks.x = element_blank(),

        legend.position = "right",
        legend.justification = "center",
        legend.box = "vertical",
        legend.box.just = "center",
        legend.title = element_text(size = base_size),
        legend.text = element_text(size = base_size - 1),
        legend.key.size = unit(0.25, "cm"),

        plot.title = element_text(hjust = 0.5, face = "bold", size = base_size),

        panel.grid = element_blank(),
        panel.border = element_blank()
      )
  )
}
