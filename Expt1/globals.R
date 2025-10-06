invisible(
  c(
    "MCMCglmm", #
    "readxl", #
    "dplyr", #
    "tidyr",
    "purrr", #
    "stringr",
    "ggplot2",
    "cowplot",
    "forcats",
    "ggpubr",
    "tibble",
    "viridis"
  ) |>
    lapply(function(x) {
      if (suppressMessages(!require(x, character.only = TRUE))) {
        install.packages(x)
        library(x, character.only = TRUE)
      }
    })
)
