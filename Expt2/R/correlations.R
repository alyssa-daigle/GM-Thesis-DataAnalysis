source("config_paths.R")
source("globals.R")
source("thesis_theme.R")

data <- read_excel(
    file.path(datapath, "CN_ratio_data.xlsx"),
    sheet = "CN_applied"
)

pond_name_mapping <- c(
    "MP_1" = "Mill Pond",
    "ODR_2" = "Dairy Farm Pond 1",
    "ODR_3" = "Dairy Farm Pond 2",
    "TF_1" = "Thompson Farm Pond 1",
    "TF_2" = "Thompson Farm Pond 2",
    "UM_1" = "Upper Mill Pond"
)

pond_levels <- c(
    "MP_1" = "Mill Pond",
    "ODR_2" = "Dairy Farm Pond 1",
    "ODR_3" = "Dairy Farm Pond 2",
    "TF_1" = "Thompson Farm Pond 1",
    "TF_2" = "Thompson Farm Pond 2",
    "UM_1" = "Upper Mill Pond"
)

# Clean and separate CN data
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

# Load mass data

mass_dat <- read.csv(file.path(
    datapath,
    "prod_traits.csv"
))

mass_separated <- mass_dat |>
    dplyr::select(sample, mass_g) |>
    separate(
        sample,
        into = c("location_pond", "replicate"),
        sep = "-(?=[^-.]+$)",
        extra = "merge"
    ) |>
    mutate(location_pond = gsub("-", "_", location_pond))

# Combine mass and CN data
combined_massCN_data <- merge(
    mass_separated,
    df_separated,
    by = c("location_pond", "replicate")
)

# Load microcystin data
mc_dat <- read.csv(file.path(
    datapath,
    "microcystin_lettuce.csv"
))
mc_separated <- mc_dat |>
    dplyr::select(sample, ng_g_let) |>
    extract(
        sample,
        into = c("location_pond", "replicate", "sample_type"),
        regex = "(.*-\\d+)-(\\d+)-([A-Z]+)",
        remove = FALSE
    ) |>
    mutate(location_pond = gsub("-", "_", location_pond))

# Combine mass and microcystin data
combined_massMC_data <- merge(
    mass_separated,
    mc_separated,
    by = c("location_pond", "replicate")
)

# Add the 'pond_label' to combined_massMC_data
combined_massMC_data <- combined_massMC_data |>
    mutate(pond_label = pond_name_mapping[location_pond])

# View the first few rows of the combined data
head(combined_massMC_data)

combined_MCCN_data <- merge(
    mc_separated,
    df_separated,
    by = c("location_pond", "replicate")
)

#drop NAs from combined data
combined_massCN_data <- combined_massCN_data |>
    filter(!is.na(mass_g) & !is.na(CN))
combined_MCCN_data <- combined_MCCN_data |>
    filter(!is.na(ng_g_let) & !is.na(CN))
combined_massMC_data <- combined_massMC_data |>
    filter(!is.na(mass_g) & !is.na(ng_g_let))


# Create a custom mapping for pond_label to specific shapes (e.g., 16, 17, etc.)
pond_labels <- c(
    "Mill Pond" = 0,
    "Dairy Farm Pond 1" = 1,
    "Dairy Farm Pond 2" = 2,
    "Thompson Farm Pond 1" = 3,
    "Thompson Farm Pond 2" = 4,
    "Upper Mill Pond" = 5
)

# Prepare JPEG output
jpeg("Expt2_correlation_MCMCglmm.jpeg", width = 1100, height = 1200, res = 150)

# Set layout: 3 rows, 1 column with space for legend on the right
layout(
    matrix(c(1, 2, 3), nrow = 3, byrow = TRUE),
    widths = c(4, 4, 4),
    heights = c(1, 1, 1)
)
par(
    mar = c(5, 5, 2, 15),
    oma = c(0, 0, 0, 5.7),
    cex.axis = 1.3,
    cex.lab = 1.3,
    cex.main = 1.3,
    xpd = TRUE
)

## === 1. Mass vs CN ===
mass_CN_mod <- MCMCglmm(
    mass_g ~ CN,
    data = combined_massCN_data,
    verbose = FALSE,
    nitt = 10100,
    thin = 10,
    burnin = 100
)
cn_seq <- seq(
    from = min(combined_massCN_data$CN),
    to = max(combined_massCN_data$CN),
    length.out = 1000
)
mod_sol <- mass_CN_mod$Sol
mean_pred <- sapply(cn_seq, function(x) mean(mod_sol[, 1] + mod_sol[, 2] * x))
hpdi_pred <- sapply(cn_seq, function(x) {
    HPDinterval(mod_sol[, 1] + mod_sol[, 2] * x, prob = 0.95)
})
ymin <- min(c(hpdi_pred, combined_massCN_data$mass_g))
ymax <- max(c(hpdi_pred, combined_massCN_data$mass_g))
plot(
    cn_seq,
    mean_pred,
    type = "l",
    lwd = 2,
    col = "blue",
    ylim = c(ymin, ymax),
    xlab = "Green Manure C:N",
    ylab = "Lettuce Head Mass (g)"
)
polygon(
    c(cn_seq, rev(cn_seq)),
    c(hpdi_pred[1, ], rev(hpdi_pred[2, ])),
    col = rgb(0, 0, 1, 0.2),
    border = NA
)
points(
    combined_massCN_data$CN,
    combined_massCN_data$mass_g,
    pch = pond_labels[combined_massCN_data$pond_label],
    col = "black"
)
pMCMC_mass_CN <- 2 * min(mean(mod_sol[, 2] > 0), mean(mod_sol[, 2] < 0))
text(
    x = 16,
    y = 34,
    sprintf("pMCMC = %.3f (n.s.)", pMCMC_mass_CN),
    cex = 1.2,
    col = "blue"
)


## === 2. MC vs CN ===
mc_CN_mod <- MCMCglmm(
    ng_g_let ~ CN,
    data = combined_MCCN_data,
    verbose = FALSE,
    nitt = 10100,
    thin = 10,
    burnin = 100
)
mod_sol <- mc_CN_mod$Sol
mean_pred <- sapply(cn_seq, function(x) mean(mod_sol[, 1] + mod_sol[, 2] * x))
hpdi_pred <- sapply(cn_seq, function(x) {
    HPDinterval(mod_sol[, 1] + mod_sol[, 2] * x, prob = 0.95)
})
ymin <- min(c(hpdi_pred, combined_MCCN_data$ng_g_let))
ymax <- max(c(hpdi_pred, combined_MCCN_data$ng_g_let))
plot(
    cn_seq,
    mean_pred,
    type = "l",
    lwd = 2,
    col = "blue",
    ylim = c(ymin, ymax),
    xlab = "Green Manure C:N",
    ylab = "Microcystin (µg/g lettuce)"
)
polygon(
    c(cn_seq, rev(cn_seq)),
    c(hpdi_pred[1, ], rev(hpdi_pred[2, ])),
    col = rgb(0, 0, 1, 0.2),
    border = NA
)
points(
    combined_MCCN_data$CN,
    combined_MCCN_data$ng_g_let,
    pch = pond_labels[combined_MCCN_data$pond_label],
    col = "black"
)
pMCMC_MC_CN <- 2 * min(mean(mod_sol[, 2] > 0), mean(mod_sol[, 2] < 0))
text(
    x = 16,
    y = .55,
    sprintf("pMCMC = %.3f (n.s.)", pMCMC_MC_CN),
    cex = 1.2,
    col = "blue"
)

# Add the shared legend outside the second plot
legend(
    "right",
    inset = c(-0.5, 0),
    legend = names(pond_labels),
    pch = pond_labels,
    title = "Green Manure Source",
    cex = 1.2,
    xpd = TRUE
)

## === 3. Mass vs MC ===
mass_MC_mod <- MCMCglmm(
    mass_g ~ ng_g_let,
    data = combined_massMC_data,
    verbose = FALSE,
    nitt = 10100,
    thin = 10,
    burnin = 100
)
mc_seq <- seq(
    from = min(combined_massMC_data$ng_g_let),
    to = max(combined_massMC_data$ng_g_let),
    length.out = 1000
)
mod_sol <- mass_MC_mod$Sol
mean_pred <- sapply(mc_seq, function(x) mean(mod_sol[, 1] + mod_sol[, 2] * x))
hpdi_pred <- sapply(mc_seq, function(x) {
    HPDinterval(mod_sol[, 1] + mod_sol[, 2] * x, prob = 0.95)
})
ymin <- min(c(hpdi_pred, combined_massMC_data$mass_g))
ymax <- max(c(hpdi_pred, combined_massMC_data$mass_g))
plot(
    mc_seq,
    mean_pred,
    type = "l",
    lwd = 2,
    col = "blue",
    ylim = c(ymin, ymax),
    xlab = "Microcystin (µg/g lettuce)",
    ylab = "Lettuce Head Mass (g)"
)
polygon(
    c(mc_seq, rev(mc_seq)),
    c(hpdi_pred[1, ], rev(hpdi_pred[2, ])),
    col = rgb(0, 0, 1, 0.2),
    border = NA
)
points(
    combined_massMC_data$ng_g_let,
    combined_massMC_data$mass_g,
    pch = pond_labels[combined_massMC_data$pond_label],
    col = "black"
)
pMCMC_mass_MC <- 2 * min(mean(mod_sol[, 2] > 0), mean(mod_sol[, 2] < 0))
text(
    x = 0.26,
    y = 11,
    sprintf("pMCMC = %.3f (n.s.)", pMCMC_mass_MC),
    cex = 1.2,
    col = "blue"
)

# Finish JPEG
dev.off()
