source("config_paths.R")
source("globals.R")
source("thesis_theme.R")

# load data
elisa_dat <- read.csv(file.path(
  datapath,
  "microcystin_lettuce.csv"
)) |>
  separate(
    sample,
    into = c("loc", "pond", "rep", "type"),
    sep = "-",
    fill = "right"
  ) |>
  mutate(ng_g_let = as.numeric(ng_g_let)) |>
  mutate(loc_pond = paste(loc, pond, sep = "_")) |>
  drop_na()

elisa_dat$loc_pond <- factor(
  elisa_dat$loc_pond,
  levels = c(
    "fert_1",
    "cntr_1",
    "MP_1",
    "UM_1",
    "ODR_2",
    "ODR_3",
    "TF_1",
    "TF_2"
  )
)

elisa_dat <- elisa_dat |>
  mutate(
    pooled_group = case_when(
      loc_pond %in% c("fert_1", "cntr_1") ~ loc_pond,
      TRUE ~ "Pooled"
    ),
    pooled_group = factor(
      pooled_group,
      levels = c("fert_1", "cntr_1", "Pooled")
    )
  )

loc_labels <- c(
  "fert_1" = "Conventional \nFertilizer",
  "cntr_1" = "Control",
  "ODR_2" = "Dairy Farm \nPond 1",
  "ODR_3" = "Dairy Farm \nPond 2",
  "TF_1" = "Thompson \nFarm Pond 1",
  "TF_2" = "Thompson \nFarm Pond 2",
  "MP_1" = "Mill Pond",
  "UM_1" = "Upper \nMill Pond"
)


MC_summary <- elisa_dat %>%
  group_by(loc_pond, pooled_group) %>%
  summarise(
    mean_MC = mean(ng_g_let, na.rm = TRUE),
    sd_MC = sd(ng_g_let, na.rm = TRUE)
  )


lettuce_pooled <- elisa_dat |>
  ggplot(aes(x = loc_pond, y = ng_g_let, color = loc_pond)) +
  geom_jitter(size = 2, width = 0.2, shape = 1, aes(color = loc_pond)) +
  geom_point(
    data = MC_summary,
    aes(x = loc_pond, y = mean_MC),
    size = 3,
    shape = 19,
    inherit.aes = TRUE
  ) +
  geom_errorbar(
    data = MC_summary,
    aes(
      x = loc_pond,
      ymin = mean_MC - sd_MC,
      ymax = mean_MC + sd_MC,
      color = loc_pond
    ),
    width = 0.2,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "fert_1" = "#1f78b4",
      "cntr_1" = "firebrick",
      "MP_1" = "aquamarine4",
      "UM_1" = "aquamarine4",
      "ODR_2" = "aquamarine4",
      "ODR_3" = "aquamarine4",
      "TF_1" = "aquamarine4",
      "TF_2" = "aquamarine4"
    )
  ) +
  scale_color_manual(
    values = c(
      "fert_1" = "#1f78b4",
      "cntr_1" = "firebrick",
      "MP_1" = "aquamarine4",
      "UM_1" = "aquamarine4",
      "ODR_2" = "aquamarine4",
      "ODR_3" = "aquamarine4",
      "TF_1" = "aquamarine4",
      "TF_2" = "aquamarine4"
    )
  ) +
  guides(fill = "none") +
  scale_x_discrete(labels = loc_labels) +
  labs(
    x = "Treatment Group",
    y = "Microcystin (µg/g Lettuce)",
    shape = "Green Manure \nSource"
  ) +
  theme_cowplot() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 9),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8)
  )
# geom_hline(yintercept = 0.1, linetype = "dashed", color = "red", size = 1) +
# annotate(
#     "text",
#     x = 1.2,
#     y = 0.12,
#     label = "Limit of Detection",
#     color = "red",
#     size = 3
# )

lettuce_pooled

# ggsave(
#   filename = file.path(plots_folder, "Expt2_lettuce_elisa.jpg"),
#   width = 7,
#   height = 4
# )

mod <- MCMCglmm(
  data = elisa_dat,
  ng_g_let ~ -1 + loc_pond,
  nitt = 101000,
  thin = 10,
  burnin = 1000,
  verbose = FALSE
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
