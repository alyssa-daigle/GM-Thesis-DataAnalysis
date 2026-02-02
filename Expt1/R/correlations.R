# using source() is basically "calling" those other scripts, which loads anything within those scripts

source("config_paths.R") # loads the necessary file paths
source("globals.R") # load the librarys needed for this script
source("thesis_theme.R") # loads the plot theming
source(file.path(path, "R/variance_explained.R")) # loads the ssbyvar variance function

# Read and summarize data

# microcystin data
elisa <- read.csv(file.path(
    data,
    "microcystin.csv"
)) |>
    select(treatment, MC) |>
    group_by(treatment) |>
    summarise(mean_MC = mean(MC, na.rm = TRUE))

# total phosphorus data
tp <- read.csv(file.path(
    data,
    "total_phosphorus.csv"
)) |>
    select(treatment, ppb) |>
    group_by(treatment) |>
    summarise(mean_tp = mean(ppb, na.rm = TRUE))

# relative (duckweed) growth rate data
rgr <- read.csv(file.path(
    data,
    "frondarea.csv"
)) |>
    select(treatment, RGR) |>
    group_by(treatment) |>
    summarise(mean_rgr = mean(RGR, na.rm = TRUE))

# Merge all
combined_dat <- elisa |>
    left_join(tp, by = "treatment") |>
    left_join(rgr, by = "treatment")


# Save to JPEG file
# using base R to plot these
jpeg(
    file.path(plots_folder, "Expt1_correlations.jpeg"),
    width = 800,
    height = 1200,
    res = 150
)

# Set up plotting parameters: 3 rows, 1 column; adjust margins; text size
par(
    mfrow = c(3, 1),
    mar = c(5, 5, 2, 2),
    cex.axis = 1.3,
    cex.lab = 1.3,
    cex.main = 1.3
)

## === Plot 1: RGR vs TP ===
rgr_tp_mod <- MCMCglmm(
    mean_rgr ~ mean_tp,
    data = combined_dat,
    verbose = F,
    nitt = 10100,
    thin = 10,
    burnin = 100
)
tp_seq <- seq(
    from = min(combined_dat$mean_tp),
    to = max(combined_dat$mean_tp),
    length.out = 1000
)
mod_sol <- rgr_tp_mod$Sol
mean_pred <- sapply(tp_seq, function(x) mean(mod_sol[, 1] + mod_sol[, 2] * x))
hpdi_pred <- sapply(tp_seq, function(x) {
    HPDinterval(mod_sol[, 1] + mod_sol[, 2] * x, prob = 0.95)
})
ymin <- min(c(hpdi_pred, combined_dat$mean_rgr))
ymax <- max(c(hpdi_pred, combined_dat$mean_rgr))
plot(
    tp_seq,
    mean_pred,
    type = "l",
    lwd = 2,
    col = "blue",
    ylim = c(ymin, ymax),
    xlab = "Mean Total Phosphorus (µg/L)",
    ylab = "Mean Relative Growth Rate"
)
polygon(
    c(tp_seq, rev(tp_seq)),
    c(hpdi_pred[1, ], rev(hpdi_pred[2, ])),
    col = rgb(0, 0, 1, 0.2),
    border = NA
)
points(combined_dat$mean_tp, combined_dat$mean_rgr, pch = 16, col = "black")

# Extract pMCMC for RGR vs TP model
pMCMC_rgr_tp <- 2 * min(mean(mod_sol[, 2] > 0), mean(mod_sol[, 2] < 0))
text(
    x = 39000,
    y = 6.45,
    sprintf("pMCMC = %.3f", pMCMC_rgr_tp),
    cex = 1.2,
    col = "blue"
)


## === Plot 2: MC vs TP ===
mc_tp_mod <- MCMCglmm(
    mean_MC ~ mean_tp,
    data = combined_dat,
    verbose = F,
    nitt = 10100,
    thin = 10,
    burnin = 100
)
mod_sol <- mc_tp_mod$Sol
mean_pred <- sapply(tp_seq, function(x) mean(mod_sol[, 1] + mod_sol[, 2] * x))
hpdi_pred <- sapply(tp_seq, function(x) {
    HPDinterval(mod_sol[, 1] + mod_sol[, 2] * x, prob = 0.95)
})
ymin <- min(c(hpdi_pred, combined_dat$mean_MC))
ymax <- max(c(hpdi_pred, combined_dat$mean_MC))
plot(
    tp_seq,
    mean_pred,
    type = "l",
    lwd = 2,
    col = "blue",
    ylim = c(ymin, ymax),
    xlab = "Mean Total Phosphorus (µg/L)",
    ylab = "Mean Microcystin (µg/g Duckweed)"
)
polygon(
    c(tp_seq, rev(tp_seq)),
    c(hpdi_pred[1, ], rev(hpdi_pred[2, ])),
    col = rgb(0, 0, 1, 0.2),
    border = NA
)
points(combined_dat$mean_tp, combined_dat$mean_MC, pch = 16, col = "black")

pMCMC_mc_tp <- 2 * min(mean(mod_sol[, 2] > 0), mean(mod_sol[, 2] < 0))

# Extract pMCMC for MC vs TP model
text_label <- if (pMCMC_mc_tp >= 0.05) {
    sprintf("pMCMC = %.3f (n.s.)", pMCMC_mc_tp)
} else {
    sprintf("pMCMC = %.3f", pMCMC_mc_tp)
}
text(x = 38500, y = 9, text_label, cex = 1.2, col = "blue")


## === Plot 3: RGR vs MC ===
rgr_mc_mod <- MCMCglmm(
    mean_rgr ~ mean_MC,
    data = combined_dat,
    verbose = F,
    nitt = 10100,
    thin = 10,
    burnin = 100
)
mc_seq <- seq(
    from = min(combined_dat$mean_MC),
    to = max(combined_dat$mean_MC),
    length.out = 1000
)
mod_sol <- rgr_mc_mod$Sol
mean_pred <- sapply(mc_seq, function(x) mean(mod_sol[, 1] + mod_sol[, 2] * x))
hpdi_pred <- sapply(mc_seq, function(x) {
    HPDinterval(mod_sol[, 1] + mod_sol[, 2] * x, prob = 0.95)
})
ymin <- min(c(hpdi_pred, combined_dat$mean_rgr))
ymax <- max(c(hpdi_pred, combined_dat$mean_rgr))
plot(
    mc_seq,
    mean_pred,
    type = "l",
    lwd = 2,
    col = "blue",
    ylim = c(ymin, ymax),
    xlab = "Mean Microcystin (µg/g Duckweed)",
    ylab = "Mean Relative Growth Rate"
)
polygon(
    c(mc_seq, rev(mc_seq)),
    c(hpdi_pred[1, ], rev(hpdi_pred[2, ])),
    col = rgb(0, 0, 1, 0.2),
    border = NA
)
points(combined_dat$mean_MC, combined_dat$mean_rgr, pch = 16, col = "black")

# Extract pMCMC for RGR vs MC model
pMCMC_rgr_mc <- 2 * min(mean(mod_sol[, 2] > 0), mean(mod_sol[, 2] < 0))
text(
    x = 9,
    y = 5.9,
    sprintf("pMCMC = %.3f", pMCMC_rgr_mc),
    cex = 1.2,
    col = "blue"
)

# Finish JPEG
dev.off()
