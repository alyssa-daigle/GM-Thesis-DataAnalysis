library(ggplot2)
library(tidyverse)
library(lubridate)
library(cowplot)
library(leaflet)
library(mapview)
library(webshot)

setwd(
    "~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/Expt1"
)

#making custom df
locations <- c(
    "Durham Reservoir (DR)",
    "LaRoche Pond (LR)",
    "Mill Pond (M)",
    "Thompson Farm (TF)",
    "Upper Mill Pond (UM)",
    "Woodman Rd (W)",
    "Kingman Farm",
    "Organic Dairy Research Farm",
    "UNH Durham Campus"
)

longitude <- c(
    -70.94412859656731,
    -70.94562421030663,
    -70.92152623,
    -70.94553034136233,
    -70.92529325317224,
    -70.91917424282973,
    -70.9331772897215,
    -70.99231110992503,
    -70.93700374516798
)

latitude <- c(
    43.147559724237034,
    43.120993650017304,
    43.130809896234574,
    43.10826058413447,
    43.12248009213886,
    43.1373552934115,
    43.16780084737212,
    43.097277151100016,
    43.13876007949401
)

sampling_locations <- data.frame(
    Location = locations,
    long = longitude,
    lat = latitude
)

color_list <- c(
    "#EE6677", # coral
    "#228833", # green
    "#4477AA", # blue
    "#CCBB44", # olive
    "#66CCEE", # light blue
    "#AA3377"
) # purple


# Split the data
gm_locations <- sampling_locations %>%
    filter(Location %in% c("Kingman Farm", "Organic Dairy Research Farm"))

unh_location <- sampling_locations %>%
    filter(Location == "UNH Durham Campus")

regular_locations <- sampling_locations %>%
    filter(
        !(Location %in%
            c(
                "Kingman Farm",
                "Organic Dairy Research Farm",
                "UNH Durham Campus"
            ))
    )


# Define the custom star icon
star_icon <- makeIcon(
    iconUrl = "star.png", # path to your downloaded image
    iconWidth = 30, # adjust size as needed
    iconHeight = 30,
    iconAnchorX = 15, # center anchor
    iconAnchorY = 15
)

farm_icon <- makeIcon(
    iconUrl = "farm.png", # path to your downloaded image
    iconWidth = 30, # adjust size as needed
    iconHeight = 30,
    iconAnchorX = 15, # center anchor
    iconAnchorY = 15
)

# Now make your map
sampling_map <- leaflet() |>
    addProviderTiles(providers$CartoDB.Positron) |>

    # Add regular circle markers for the others
    addCircleMarkers(
        data = regular_locations,
        lng = ~long,
        lat = ~lat,
        color = color_list[1:6],
        opacity = 1
    ) |>

    # Add the UNH location using the custom star image
    addMarkers(
        data = unh_location,
        lng = ~long,
        lat = ~lat,
        icon = star_icon,
        label = ~Location,
        labelOptions = labelOptions(
            noHide = TRUE,
            direction = "top",
            textOnly = TRUE
        )
    ) |>

    addMarkers(
        data = gm_locations,
        lng = ~long,
        lat = ~lat,
        icon = farm_icon,
        label = ~Location,
        labelOptions = labelOptions(
            noHide = TRUE,
            direction = "top",
            textOnly = TRUE
        )
    ) |>

    # Add scale bar
    addScaleBar(position = "topright") |>

    # Add legend
    addLegend(
        position = "bottomright",
        colors = color_list,
        labels = regular_locations$Location,
        opacity = 1,
        title = "Sampling Locations"
    )

sampling_map
# Save the map as an HTML file
mapshot(sampling_map, file = "sampling_map.jpg")
