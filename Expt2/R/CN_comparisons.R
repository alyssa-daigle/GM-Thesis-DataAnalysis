source("config_paths.R")
source("globals.R")
source("thesis_theme.R")

data <- read_excel(
  file.path(datapath, "CN_ratio_data.xlsx"),
  sheet = "CN_applied"
)

pond_name_mapping <- c(
  "MP_1" = "Mill Pond",
  "ODR_2" = "Dairy Farm\nPond 1",
  "ODR_3" = "Dairy Farm\nPond 2",
  "TF_1" = "Thompson \nFarm Pond 1",
  "TF_2" = "Thompson \nFarm Pond 2",
  "UM_1" = "Upper \nMill Pond"
)

pond_levels <- c(
  "Mill Pond",
  "Upper \nMill Pond",
  "Dairy Farm\nPond 1",
  "Dairy Farm\nPond 2",
  "Thompson \nFarm Pond 1",
  "Thompson \nFarm Pond 2"
)

# Clean and separate
df_separated <- data |>
  dplyr::select(sample, CN) |>
  separate(
    sample,
    into = c("location_pond", "replicate"),
    sep = "-(?=[^-.]+$)",
    extra = "merge"
  ) |>
  mutate(
    location_pond = gsub("-", "_", location_pond),
    pond_label = pond_name_mapping[location_pond],
    pond_label = factor(pond_label, levels = pond_levels)
  )

# Custom colors
custom_colors <- c(
  "#44AA99",
  "#CC6677",
  "#88CCEE",
  "#117733",
  "#DDCC77",
  "#332288"
)

# Calculate summary statistics (mean and standard deviation) for CN by pond_label
grouped_summary_df <- df_separated |>
  group_by(pond_label) |>
  summarise(
    mean_CN = mean(CN, na.rm = TRUE),
    sd_CN = sd(CN, na.rm = TRUE),
    .groups = "drop"
  )

# Join the summary statistics back into df_separated
df_separated <- df_separated |>
  left_join(grouped_summary_df, by = "pond_label")

# Plotting
GM_CN_plot <- ggplot(data = df_separated, aes(x = pond_label, y = CN)) +
  geom_jitter(size = 2, width = 0.2, shape = 1) +
  geom_point(aes(x = pond_label, y = mean_CN), size = 3, shape = 19) + # Changed shape to a filled circle for better visibility
  geom_errorbar(
    aes(x = pond_label, ymin = mean_CN - sd_CN, ymax = mean_CN + sd_CN),
    width = 0.15
  ) +
  labs(x = "Green Manure Source", y = "Green Manure C:N Applied") + # General label change for clarity
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 8), legend.position = "none")

# Display plot
print(GM_CN_plot)

# Save plot
ggsave(
  file.path("plots", "Expt2_GM_CN_plot.jpg"),
  GM_CN_plot,
  width = 6,
  height = 4,
  dpi = 500
)

# Running MCMC model
mod <- MCMCglmm(
  CN ~ -1 + location_pond,
  data = df_separated,
  verbose = F,
  nitt = 101000,
  thin = 10,
  burnin = 1000
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
