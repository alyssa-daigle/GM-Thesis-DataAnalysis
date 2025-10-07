library(ggplot2)
library(grid)

#pond mapping and levels
pond_name_mapping <- c(
  "MP-1" = "Mill Pond",
  "ODR-2" = "Dairy Farm\nPond 1",
  "ODR-3" = "Dairy Farm\nPond 2",
  "TF-1" = "Thompson \nFarm Pond 1",
  "TF-2" = "Thompson \nFarm Pond 2",
  "UM-1" = "Upper \nMill Pond"
)

pond_levels <- c(
  "Mill Pond",
  "Upper \nMill Pond",
  "Dairy Farm\nPond 1",
  "Dairy Farm\nPond 2",
  "TF-1" = "Thompson \nFarm Pond 1",
  "TF-2" = "Thompson \nFarm Pond 2"
)

# theme.R

library(ggplot2)
library(cowplot)
library(ggtext)

# Define a consistent fill scale for duckweed species
duckweed_fill_scale <- function() {
  scale_fill_manual(
    values = c(
      "Lemna" = "#9ACD32",
      "Spirodela" = "#20B2AA",
      "Wolffia" = "#9370DB"
    ),
    labels = c(
      "Lemna" = "*Lemna*",
      "Spirodela" = "*Spirodela*",
      "Wolffia" = "*Wolffia*"
    )
  )
}

# Define your standard plot theme
thesis_theme <- function() {
  theme_cowplot() +
    theme(
      axis.text.x = element_text(size = 9),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_markdown(size = 10)
    )
}
