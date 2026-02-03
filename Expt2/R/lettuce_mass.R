source("config_paths.R")
source("globals.R")
source("thesis_theme.R")

mass_dat <- read.csv(file.path(
  datapath,
  "prod_traits.csv"
)) |>
  mutate(
    height_cm = as.numeric(height_cm),
    width_cm = as.numeric(width_cm),
    mass_g = as.numeric(mass_g)
  )

mass_dat$group <- sub("-.*", "", mass_dat$sample)

mass_dat <- na.omit(mass_dat)


mass_dat$date <- as.factor(mass_dat$date)
mass_dat$group <- as.factor(mass_dat$group)

mass_dat <- mass_dat %>%
  mutate(
    group = case_when(
      grepl("ODR-2-", sample) ~ "ODR-2",
      grepl("ODR-3-", sample) ~ "ODR-3",
      grepl("TF-1-", sample) ~ "TF-1",
      grepl("TF-2-", sample) ~ "TF-2",
      grepl("MP-1-", sample) ~ "MP-1",
      grepl("UM-1-", sample) ~ "UM-1",
      TRUE ~ group
    )
  )

mass_dat <- mass_dat %>%
  mutate(
    group = factor(
      group,
      levels = c(
        "FERT",
        "CNTR",
        "ODR-2",
        "ODR-3",
        "TF-1",
        "TF-2",
        "MP-1",
        "UM-1"
      )
    )
  )

loc_labels <- c(
  "FERT" = "Conventional \nFertilizer",
  "CNTR" = "Control",
  "ODR-2" = "Dairy Farm \nPond 1",
  "ODR-3" = "Dairy Farm \nPond 2",
  "TF-1" = "Thompson \nFarm Pond 1",
  "TF-2" = "Thompson \nFarm Pond 2",
  "MP-1" = "Mill Pond",
  "UM-1" = "Upper \nMill Pond"
)

mass_dat <- mass_dat |>
  mutate(
    pooled_group = case_when(
      group %in% c("FERT", "CNTR") ~ group,
      TRUE ~ "Pooled"
    ),
    pooled_group = factor(pooled_group, levels = c("FERT", "CNTR", "Pooled")) # Set the order
  )

mass_dat <- mass_dat %>%
  mutate(
    group = factor(
      group,
      levels = c(
        "FERT",
        "CNTR",
        "MP-1",
        "UM-1",
        "ODR-2",
        "ODR-3",
        "TF-1",
        "TF-2"
      )
    )
  )

# Summarize mean and SD for each group
mass_summary <- mass_dat %>%
  group_by(group) %>%
  summarise(
    mean_mass = mean(mass_g, na.rm = TRUE),
    sd_mass = sd(mass_g, na.rm = TRUE)
  )

# Plot
plot <- ggplot(mass_dat, aes(x = group, y = mass_g, fill = group)) +
  geom_jitter(size = 2, width = 0.2, shape = 1, aes(color = group)) +
  geom_point(
    data = mass_summary,
    aes(x = group, y = mean_mass, color = group),
    size = 3,
    shape = 19,
    inherit.aes = TRUE
  ) +
  geom_errorbar(
    data = mass_summary,
    aes(
      x = group,
      ymin = mean_mass - sd_mass,
      ymax = mean_mass + sd_mass,
      color = group
    ),
    width = 0.2,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "FERT" = "#1f78b4",
      "CNTR" = "firebrick",
      "MP-1" = "aquamarine4",
      "UM-1" = "aquamarine4",
      "ODR-2" = "aquamarine4",
      "ODR-3" = "aquamarine4",
      "TF-1" = "aquamarine4",
      "TF-2" = "aquamarine4"
    )
  ) +
  scale_color_manual(
    values = c(
      "FERT" = "#1f78b4",
      "CNTR" = "firebrick",
      "MP-1" = "aquamarine4",
      "UM-1" = "aquamarine4",
      "ODR-2" = "aquamarine4",
      "ODR-3" = "aquamarine4",
      "TF-1" = "aquamarine4",
      "TF-2" = "aquamarine4"
    )
  ) +
  guides(fill = "none", color = "none") +
  scale_x_discrete(labels = loc_labels) +
  labs(
    x = "Treatment Group",
    y = "Mass of Lettuce Head (g)",
    shape = "Green Manure \nSource"
  ) +
  theme_cowplot() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(hjust = .5),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.title.x.top = element_blank()
  ) +
  scale_y_break(c(36, 140))

plot

# ggsave(
#   file.path("plots", "Expt2_mass_pooled.jpg"),
#   plot,
#   width = 7,
#   height = 4,
#   dpi = 500
# )

mod <- MCMCglmm(
  data = mass_dat,
  mass_g ~ -1 + group,
  nitt = 101000,
  thin = 10,
  burnin = 1000,
  verbose = FALSE
)
summary(mod)

plot <- ggplot(mass_dat, aes(x = group, y = mass_g, fill = group)) +
  geom_jitter(size = 2, width = 0.2, shape = 1, aes(color = group)) +
  geom_point(
    data = mass_summary,
    aes(x = group, y = mean_mass, color = group),
    size = 3,
    shape = 19,
    inherit.aes = TRUE
  ) +
  geom_errorbar(
    data = mass_summary,
    aes(
      x = group,
      ymin = mean_mass - sd_mass,
      ymax = mean_mass + sd_mass,
      color = group
    ),
    width = 0.2,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "FERT" = "#1f78b4",
      "CNTR" = "firebrick",
      "MP-1" = "aquamarine4",
      "UM-1" = "aquamarine4",
      "ODR-2" = "aquamarine4",
      "ODR-3" = "aquamarine4",
      "TF-1" = "aquamarine4",
      "TF-2" = "aquamarine4"
    )
  ) +
  scale_color_manual(
    values = c(
      "FERT" = "#1f78b4",
      "CNTR" = "firebrick",
      "MP-1" = "aquamarine4",
      "UM-1" = "aquamarine4",
      "ODR-2" = "aquamarine4",
      "ODR-3" = "aquamarine4",
      "TF-1" = "aquamarine4",
      "TF-2" = "aquamarine4"
    )
  ) +
  guides(fill = "none", color = "none") +
  scale_x_discrete(labels = loc_labels) +
  labs(
    x = "Treatment Group",
    y = "Mass of Lettuce Head (g)",
    shape = "Green Manure \nSource"
  ) +
  theme_cowplot() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(hjust = .5),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.title.x.top = element_blank()
  ) +
  ylim(0, 39)

plot

# ggsave(
#   file.path("plots", "Expt2_mass_pooled_noCF.jpg"),
#   plot,
#   width = 7,
#   height = 4,
#   dpi = 500
# )
