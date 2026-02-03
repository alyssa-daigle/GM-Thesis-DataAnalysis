source("config_paths.R")
source("globals.R")
source("thesis_theme.R")

# "n" is for neuston in this case, so fresh duckweed
elisa_n <- read.csv(file.path(
  datapath,
  "microcystin_neuston.csv"
)) |>
  separate(
    sample,
    into = c("loc", "pond", "rep", "type"),
    sep = "-",
    fill = "right"
  ) |>
  mutate(ppb = as.numeric(ppb)) |>
  mutate(loc_pond = paste(loc, pond, sep = "_")) |>
  drop_na()

loc_labels <- c(
  "fert_1" = "Conventional \nFertilizer",
  "cntr_1" = "Control",
  "ODR_2" = "Dairy Farm 1",
  "ODR_3" = "Dairy Farm 2",
  "TF_1" = "Thompson \nFarm 1",
  "TF_2" = "Thompson \nFarm 2",
  "MP_1" = "Mill Pond",
  "UM_1" = "Upper \nMill Pond"
)

elisa_n$loc_pond <- factor(
  elisa_n$loc_pond,
  levels = c(
    "fert_1",
    "cntr_1",
    "ODR_2",
    "ODR_3",
    "TF_1",
    "TF_2",
    "MP_1",
    "UM_1"
  )
)

elisa_n |>
  ggplot(aes(x = loc_pond, y = ppb)) +
  geom_boxplot() +
  scale_x_discrete(labels = loc_labels) +
  labs(x = "", y = "Microcystin - Neuston (ppb)") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(size = 9))


# duckweed when added to pot, where "pp" means "pre-pot"
elisa_pp <- read.csv(file.path(
  datapath,
  "microcystin_prepot.csv"
)) |>
  separate(
    sample,
    into = c("loc", "pond", "rep", "type"),
    sep = "-",
    fill = "right"
  ) |>
  mutate(ppb = as.numeric(ppb)) |>
  mutate(loc_pond = paste(loc, pond, sep = "_")) |>
  drop_na()

elisa_pp$loc_pond <- factor(
  elisa_pp$loc_pond,
  levels = c(
    "fert_1",
    "cntr_1",
    "ODR_2",
    "ODR_3",
    "TF_1",
    "TF_2",
    "MP_1",
    "UM_1"
  )
)

elisa_pp |>
  ggplot(aes(x = loc_pond, y = ppb)) +
  geom_boxplot() +
  scale_x_discrete(labels = loc_labels) +
  labs(x = "", y = "Microcystin - Prepot (ppb)") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(size = 9))

combined_data <- bind_rows(elisa_n, elisa_pp)

combined_data$loc_pond <- factor(
  combined_data$loc_pond,
  levels = c(
    "fert_1",
    "cntr_1",
    "ODR_2",
    "ODR_3",
    "TF_1",
    "TF_2",
    "MP_1",
    "UM_1"
  )
)
combined_data$condition <- factor(
  combined_data$condition,
  levels = c(
    "Duckweed at Collection (9/20/24)",
    "Duckweed at Potting (9/24/24)"
  )
)

prepot_neuston_plot <- combined_data |>
  ggplot(aes(x = loc_pond, y = ppb, color = condition)) +
  scale_x_discrete(labels = loc_labels) +
  geom_boxplot() +
  scale_color_manual(
    values = c(
      "Duckweed at Collection (9/20/24)" = "aquamarine4",
      "Duckweed at Potting (9/24/24)" = "darkblue"
    )
  ) +
  labs(x = "Location", y = "Microcystin (ppb)", color = "") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(size = 9)) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "red", size = 1)

prepot_neuston_plot

# ggsave(
#   filename = file.path(plots_folder, "pp_n_plot.png"),
#   width = 7,
#   height = 4
# )

stats_results <- combined_data |>
  group_by(loc_pond) |>
  summarise(t_test = list(t.test(ppb ~ condition)), .groups = 'drop') |>
  mutate(
    p_value = sapply(t_test, function(x) x$p.value),
    significant = ifelse(p_value < 0.05, "Yes", "No")
  )

stats_results

mod <- MCMCglmm(
  ppb ~ loc_pond:condition,
  data = combined_data,
  verbose = F,
  nitt = 100000,
  thin = 500,
  burnin = 5000
)
summary(mod)

#pooled at time of potting

combined_data |>
  filter(condition == "Duckweed at Potting (9/24/24)") |>
  ggplot(aes(x = condition, y = ppb)) +
  geom_boxplot(aes(fill = condition)) +
  scale_fill_manual(
    values = c("Duckweed at Potting (9/24/24)" = "aquamarine4")
  ) +
  scale_y_continuous(
    breaks = seq(0, max(combined_data$ppb, na.rm = TRUE), by = 0.1)
  ) +
  labs(
    x = "Duckweed Green Manure at Potting",
    y = "Microcystin (µg/L)",
    color = ""
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(), # Removes the x-axis tick label
    axis.ticks.x = element_blank()
  ) + # Removes x-axis ticks
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "red", size = 1)

# ggsave(
#   filename = file.path(plots_folder, "pooled_pp_n_plot.png"),
#   width = 7,
#   height = 4
# )

mod2 <- MCMCglmm(
  ppb ~ loc_pond:condition,
  data = combined_data,
  verbose = F,
  nitt = 100000,
  thin = 500,
  burnin = 5000
)
summary(mod2)

# plotting SW on the same plot

# strainer water from the duckweed when it first collected

elisa_sw <- read.csv(file.path(
  datapath,
  "microcystin_strainerwater.csv"
)) |>
  separate(
    sample,
    into = c("loc", "pond", "rep", "type"),
    sep = "-",
    fill = "right"
  ) |>
  mutate(ppb = as.numeric(mc_ppb)) |>
  mutate(loc_pond = paste(loc, pond, sep = "_")) |>
  drop_na()

elisa_sw <- elisa_sw |> mutate(condition = "SW") # Change condition to "SW"
sw_pp_combined_data <- bind_rows(elisa_sw, elisa_pp) # Combine the dataframes

# Set factor levels so "SW" appears first
sw_pp_combined_data <- sw_pp_combined_data |>
  mutate(
    condition = factor(
      condition,
      levels = c("SW", unique(condition[condition != "SW"]))
    )
  )

sw_pp_prepot_neuston_plot <- sw_pp_combined_data |>
  ggplot(aes(x = condition, y = ppb)) +
  geom_boxplot(aes(fill = condition)) +
  scale_fill_manual(
    values = c(
      "SW" = "#1f78b4",
      "Duckweed at Potting (9/24/24)" = "aquamarine4"
    )
  ) +
  scale_x_discrete(
    labels = c(
      "SW" = "Strainer Water",
      "Duckweed at Potting (9/24/24)" = "Duckweed Green Manure"
    )
  ) +
  scale_y_continuous(
    breaks = seq(0, max(combined_data$ppb, na.rm = TRUE), by = 0.1)
  ) +
  labs(y = "Microcystin (µg/L)", color = "") +
  theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "red", size = 1)

sw_pp_prepot_neuston_plot

# ggsave(
#   filename = file.path(plots_folder, "sw_pp_n_plot.png"),
#   width = 7,
#   height = 4
# )

# Combine all three datasets
sw_pp_n_combined <- bind_rows(elisa_sw, elisa_n, elisa_pp)

# Ensure correct factor order for x-axis
sw_pp_n_combined <- sw_pp_n_combined |>
  mutate(
    condition = factor(
      condition,
      levels = c(
        "SW",
        "Duckweed at Collection (9/20/24)",
        "Duckweed at Potting (9/24/24)"
      )
    )
  )

mod3 <- MCMCglmm(
  ppb ~ condition,
  data = sw_pp_n_combined,
  verbose = F,
  nitt = 11000,
  thin = 10,
  burnin = 1000
)
summary(mod3)

summary_df <- summary(mod3)$solutions %>%
  as.data.frame() %>%
  rownames_to_column("Effect") |>
  mutate(
    Effect = recode(
      Effect,
      "(Intercept)" = "Strainer Water (9/20) (Intercept)",
      "conditionDuckweed at Collection (9/20/24)" = "Duckweed at Collection (9/20)",
      "conditionDuckweed at Potting (9/24/24)" = "Duckweed at Potting (9/24)"
    ),
    Significance = case_when(
      pMCMC < 0.001 ~ "***",
      pMCMC < 0.01 ~ "**",
      pMCMC < 0.05 ~ "*",
      pMCMC < 0.1 ~ ".",
      TRUE ~ ""
    ),

    # Format pMCMC with stars
    pMCMC = paste0(sprintf("%.4f", pMCMC), Significance)
  ) %>%
  dplyr::select(-Significance)

# Create the table
kable(
  summary_df,
  digits = 4,
  caption = "Location effects: MC ~ Sample Type"
) %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE)

# perform ANOVA and Tukey tests
anova_results <- aov(ppb ~ condition, data = sw_pp_n_combined)
tukey_results <- TukeyHSD(anova_results)
tukey_groups <- multcompLetters(
  tukey_results$condition[, "p adj"] < 0.05
)$Letters

summary_df <- sw_pp_n_combined |>
  group_by(condition) |>
  summarise(
    mean_ppb = mean(ppb, na.rm = TRUE),
    sd_ppb = sd(ppb, na.rm = TRUE),
    max_ppb = max(ppb, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(group = tukey_groups[match(condition, names(tukey_groups))])

sw_n_pp_plot <- sw_pp_n_combined |>
  ggplot(aes(x = condition, y = ppb)) +
  geom_boxplot(aes(fill = condition)) +
  scale_fill_manual(
    values = c(
      "SW" = "#1f78b4",
      "Duckweed at Collection (9/20/24)" = "aquamarine",
      "Duckweed at Potting (9/24/24)" = "aquamarine4"
    )
  ) +
  scale_x_discrete(
    labels = c(
      "SW" = "Strainer Water \n(9/20/24)",
      "Duckweed at Collection (9/20/24)" = "Duckweed at \nCollection (9/20/24)",
      "Duckweed at Potting (9/24/24)" = "Duckweed at \nPotting (9/24/24)"
    )
  ) +
  labs(x = "", y = "Microcystin (µg/L)", color = "") +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "red", size = 1) +
  geom_text(
    data = summary_df,
    aes(x = condition, y = max_ppb + 0.1, label = group),
    inherit.aes = FALSE,
    vjust = -0.5,
    size = 4
  ) +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(size = 10))

sw_n_pp_plot

# ggsave(
#   filename = file.path(plots_folder, "Expt2_sw_n_pp_plot.jpg"),
#   width = 5,
#   height = 4
# )
