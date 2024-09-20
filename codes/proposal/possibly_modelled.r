# identify species that can be modeled

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

iho <- st_read("~/Research/mpa_europe/mpaeu_sandbox/data/Intersect_EEZ_IHO_v4_2020/Intersect_EEZ_IHO_v4_2020.shp")
europe <- st_read("~/Research/mpa_europe/mpaeu_sdm/data/shapefiles/mpa_europe_starea_v2.shp")

europe_bbox <- st_bbox(europe)
europe_bbox[1] <- -100
europe_bbox <- st_as_sfc(europe_bbox)

to_cut <- iho$MARREGION[grepl("Canadian|United States|Mexico|Bahamas|Bermu|Turks and C|Dominican|Saint-Pierre", iho$MARREGION)]
iho_sel <- iho[iho$MARREGION %in% to_cut,]

wrld <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf", continent = "North America")

europe_bbox <- st_intersection(europe_bbox, st_make_grid(europe_bbox, n = c(4,4)))

records <- list()

# Do in batches
for (i in 1:16) {
    cat("Processing", i, "\n")
    europe_h3 <- polygon_to_cells(europe_bbox[i,], res = 7)
    europe_h3 <- unlist(europe_h3)
    cells <- data.frame(cell = europe_h3)

    # Set up duckdb connection and register cells table
    con <- dbConnect(duckdb())
    dbSendQuery(con, "install httpfs; load httpfs;")
    duckdb_register(con, "cells", cells)

    # Join cells list and gridded species dataset
    # Use this one to load from the S3 bucket
#     species <- dbGetQuery(con, "
#   select *
#   from cells
#   inner join read_parquet('s3://obis-products/speciesgrids/h3_7/*') h3 on cells.cell = h3.h3_07
# ")
    species <- dbGetQuery(con, "
      select *
      from cells
      inner join read_parquet('../../mpa_europe/mpaeu_shared/grided/h3_7/*') h3 on cells.cell = h3.h3_07
    ")

    species$h3_5 <- h3jsr::get_parent(species$h3_07, res = 5)

    records[[i]] <- species

    dbDisconnect(con)
}

records <- bind_rows(records)


# Because we run in batches, just check if there are no duplicated cells
records_undup <- records %>%
    group_by(species) %>%
    distinct(cell, .keep_all = T)

# Remove those in the iho part
iho_sel_h5 <- polygon_to_cells(iho_sel, res = 5)

for (i in 1:length(iho_sel_h5)) {
  cat("Processing", i, "\n")
  records_undup <- records_undup %>%
    filter(!h3_5 %in% unlist(iho_sel_h5[[i]]))
}

# Remove those in land americas
wrld_sel_h5 <- polygon_to_cells(st_as_sf(buffer(vect(wrld), 4000)), res = 5)

for (i in 1:length(wrld_sel_h5)) {
  cat("Processing", i, "\n")
  records_undup <- records_undup %>%
    filter(!h3_5 %in% unlist(wrld_sel_h5[[i]]))
}

records_summary <- records_undup %>%
    group_by(h3_5) %>%
    summarise(records = sum(records))

europe_h3_shp <- st_as_sf(cell_to_polygon(unique(records_undup$h3_5)))

europe_h3_shp$h3 <- unique(records_undup$h3_5)
names(europe_h3_shp)[1] <- "geometry"
st_geometry(europe_h3_shp) <- "geometry"

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

# Fix those still on Americas
wrld <- st_as_sf(buffer(vect(wrld), 15000))

to_remove <- st_intersects(europe_h3_shp, wrld, sparse = F)
to_remove <- apply(to_remove, 1, any)
to_remove <- europe_h3_shp[to_remove, "h3"]

europe_h3_shp <- europe_h3_shp %>%
    filter(!h3 %in% to_remove$h3)

# Still a problem with one point
one_rm <- st_as_sfc(st_bbox(ext(-80, -75, 25, 28)))
st_crs(one_rm) <- "EPSG:4326"

to_remove <- st_intersects(europe_h3_shp, one_rm, sparse = F)
to_remove <- apply(to_remove, 1, any)
to_remove <- europe_h3_shp[to_remove, "h3"]

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

pts_per_sp <- records_undup_selarea %>%
    group_by(species, AphiaID) %>%
    summarise(available_records = n())

new_sp <- unique(records_undup_selarea$AphiaID)

mpaeu_list <- read.csv("~/Research/mpa_europe/mpaeu_sdm/data/all_splist_20240724.csv")

new_sp <- new_sp[!new_sp %in% mpaeu_list$AphiaID]
length(new_sp)

new_sp_df <- data.frame(cell = new_sp)

con <- dbConnect(duckdb())
    dbSendQuery(con, "install httpfs; load httpfs;")
    duckdb_register(con, "cells", new_sp_df)

    speciesnew <- dbGetQuery(con, "
      select species, AphiaID, h3_07
      from cells
      inner join read_parquet('../../mpa_europe/mpaeu_shared/grided/h3_7/*') h3 on cells.cell = h3.AphiaID
    ")

    dbDisconnect(con)

head(speciesnew)


pts_per_sp <- speciesnew %>%
    group_by(species, AphiaID) %>%
    summarise(available_records = n())

pts_per_sp %>%
filter(available_records > 30) %>%
nrow()
