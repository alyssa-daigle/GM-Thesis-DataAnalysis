source("config_paths.R")
source("globals.R")
source("thesis_theme.R")

# load data
species_comp <- read.csv(file.path(
    datapath,
    "species_comp.csv"
)) |>
    mutate(pond = sub("^([^-]+-[^-]+).*", "\\1", sample)) |>
    select(-c(mass_addded, sample))

#average species composition by pond
species_comp_avg <- species_comp |>
    group_by(pond) |>
    summarise(
        wolf_avg = mean(wolf_percent, na.rm = TRUE),
        lem_avg = mean(lem_percent, na.rm = TRUE),
        spiro_avg = mean(spiro_percent, na.rm = TRUE)
    )

#transform data
species_comp_long <- species_comp_avg |>
    pivot_longer(
        cols = c(wolf_avg, lem_avg, spiro_avg),
        names_to = "species",
        values_to = "avg_percent"
    ) |>
    mutate(
        species = case_when(
            species == "wolf_avg" ~ "Wolffia",
            species == "lem_avg" ~ "Lemna",
            species == "spiro_avg" ~ "Spirodela"
        ),
        pond_label = pond_name_mapping[pond],
        pond_label = factor(pond_label, levels = pond_levels)
    )

#plot data
comp_plot <- species_comp_long |>
    ggplot(aes(x = pond_label, y = avg_percent, fill = species)) +
    geom_bar(stat = "identity") +
    duckweed_fill_scale() +
    labs(
        x = "Green Manure Source",
        y = "Average Percent Abundance",
        fill = "Duckweed \nSpecies"
    ) +
    thesis_theme()

comp_plot

ggsave(
    file.path("plots", "Expt2_species_comp.jpg"),
    comp_plot,
    width = 8,
    height = 5,
    dpi = 500
)

# potential other way to represent data
# pie chart?
comp_plot <- species_comp_long |>
    ggplot(aes(x = "", y = avg_percent, fill = species)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    facet_wrap(~pond_label) +
    duckweed_fill_scale() +
    labels = c(
    "Lemna" = "*Lemna*",
    "Spirodela" = "*Spirodela*",
    "Wolffia" = "*Wolffia*"
) +
    labs(x = NULL, y = NULL, fill = "Duckweed \nSpecies") +
    theme_bw() +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 10),
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_markdown(size = 10)
    )

comp_plot
