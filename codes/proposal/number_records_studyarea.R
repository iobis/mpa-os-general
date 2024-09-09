# Plot records over Europe

library(h3jsr)
library(sf)
library(terra)
library(ggplot2)
library(dplyr)
library(DBI)
library(duckdb)
library(ggplot2)
library(patchwork)
sf::sf_use_s2(FALSE)
outfolder <- "figures"
fs::dir_create(outfolder)

europe <- st_read("~/Research/mpa_europe/mpaeu_sdm/data/shapefiles/mpa_europe_starea_v2.shp")

europe_bbox <- st_as_sfc(st_bbox(europe))
europe_bbox <- st_intersection(europe_bbox, st_make_grid(europe_bbox, n = c(2,2)))

records <- list()

# Do in batches
for (i in 1:4) {
    europe_h3 <- polygon_to_cells(europe_bbox[i,], res = 7)
    europe_h3 <- unlist(europe_h3)
    cells <- data.frame(cell = europe_h3)

    # Set up duckdb connection and register cells table
    con <- dbConnect(duckdb())
    dbSendQuery(con, "install httpfs; load httpfs;")
    duckdb_register(con, "cells", cells)

    # Join cells list and gridded species dataset
    species <- dbGetQuery(con, "
  select *
  from cells
  inner join read_parquet('s3://obis-products/speciesgrids/h3_7/*') h3 on cells.cell = h3.h3_07
")

    species$h3_5 <- h3jsr::get_parent(species$h3_07, res = 5)

    records[[i]] <- species
}

dbDisconnect(con)

records <- bind_rows(records)

# Because we run in batches, just check if there are no duplicated cells
records_undup <- records %>%
    group_by(species) %>%
    distinct(cell, .keep_all = T)

europe_h3_shp <- st_as_sf(cell_to_polygon(unique(records_undup$h3_5)))

europe_h3_shp$h3 <- unique(records_undup$h3_5)
names(europe_h3_shp)[1] <- "geometry"
st_geometry(europe_h3_shp) <- "geometry"

records_summary <- records_undup %>%
    group_by(h3_5) %>%
    summarise(records = sum(records))

europe_h3_shp <- left_join(europe_h3_shp, records_summary %>% rename(h3 = h3_5))

base_map <- rnaturalearth::ne_countries(returnclass = "sf", scale = 50)

sf_use_s2(FALSE)
base_map <- st_crop(base_map, europe)

to_remove <- st_intersects(europe_h3_shp, base_map, sparse = F)
to_remove <- apply(to_remove, 1, any)

# Check if onland
to_remove <- europe_h3_shp[to_remove, "h3"]
to_remove_crd <- h3jsr::cell_to_point(to_remove$h3)
to_remove_crd <- st_coordinates(to_remove_crd)
colnames(to_remove_crd) <- c("decimalLongitude", "decimalLatitude")
to_remove_checked <- obistools::check_onland(as.data.frame(to_remove_crd), report = T, buffer = 100)

to_remove <- to_remove[to_remove_checked$row,]

europe_h3_shp <- europe_h3_shp %>%
    filter(!h3 %in% to_remove$h3)

records_undup_selarea <- records %>%
    group_by(species) %>%
    distinct(cell, .keep_all = T) %>%
    filter(h3_5 %in% europe_h3_shp$h3)

# Get stats ----
sum(europe_h3_shp$records)
length(unique(records_undup_selarea$species))
min(records_undup_selarea$min_year, na.rm = T)
max(records_undup_selarea$max_year, na.rm = T)


# Plot ----
# Prepare base settings
sca <- function(limits, guide_title, ...) {
  scale_fill_gradientn(
    colours = rev(c("#7d1500", "#da4325", "#eca24e", "#e7e2bc", "#5cc3af", "#0a6265")),
    limits = limits,
    guide = guide_colorbar(title = guide_title,
                           show.limits = TRUE,
                           barheight = unit(0.12, "in"),
                           barwidth = unit(3.5, "in"),
                           ticks = F,
                           ticks.colour = "grey20",
                           frame.colour = "grey20",
                           title.position = "top"),
    ...#,
    #labels = c("Low", "", "", "", "High")#,
    #na.value = "#2b0700"
  )
}

wlt <- theme_classic()+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(colour = "grey60", size = 11),
        axis.title.x.top = element_text(vjust = 2, size = 16),
        axis.title.y.left = element_text(size = 16),
        panel.background = element_blank(),
        #panel.border = element_rect(fill = NA, color = "grey60"),
        panel.border = element_blank(),
        panel.grid.major = element_line(linetype = 'dashed', 
                                        colour = "grey70",
                                        linewidth = .1),
        legend.position="bottom",
        legend.title.align=0.5,
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.background = element_rect(fill = "white")
        
  )

ggplot() +
  geom_sf(data = europe_h3_shp, aes(fill = records, color = after_scale(fill)), linewidth = 0.1) +
#   sca(limits = NULL, guide_title = "log10(Number of records)",
#       trans = "log10",
#       na.value = "grey60") +
  scale_fill_viridis_c(limits = NULL, name = "log10(Number of records)",
      trans = "log10", option = "C") + 
  geom_sf(data = base_map, fill = "white", color = "grey90") +
  wlt +
  #scale_fill_viridis_c(option = "plasma", trans = "log") +
  coord_sf(crs = "EPSG:3035")

ggsave(file.path(outfolder, "number_records.jpg"),
       width = 20, height = 18, units = "cm", quality = 100)


