library(tidyverse)
library(sf)
library(rnaturalearth)

summarise_location <- function(dat, ...) {
  group_by(dat, ...) %>% summarise(across(c(lat, long), mean))
}

load_america_map <- function() {
  america_ct <- c("North America", "Central America", "South America")
  rnaturalearth::ne_countries(continent = america_ct, returnclass = "sf")
}

load_br_map <- function() {
  rnaturalearth::ne_states(country = "Brazil", returnclass = "sf")
}

base_america_map <- function(color = "#cccccc", fill = "#fefefe", map = NULL) {
  if (is.null(map))
    america_sf <- load_america_map()
  else
    america_sf <- map
  
  list(
    add_map(map = america_sf, fill = fill, color = color),
    scale_shape_manual("Age", values = c(16, 15)),
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(fill = NA),
      panel.background = element_rect(fill = "lightblue1"),
      legend.key.height = unit(0.7,"line"),
      legend.key.width = unit(1, "line"),
      legend.position = c(0.01, 0.01),
      legend.justification = c(0, 0),
      legend.box = "vertical",
      legend.direction = "horizontal",
      legend.background = element_blank()
    ),
    coord_sf(ylim = c(-55, 45), xlim = c(-120, -40))
  )
}

add_map <- function(map, color = "#cccccc", fill = "#fefefe") {
  geom_sf(data = map, fill = fill, color = color)
}

idw_america <- function(dat) {
  suppressMessages({
    require(sp)
    require(gstat)
    require(raster)
  })
  america_sf <- load_america_map()
  america_union <- sf::st_buffer(america_sf, 0.1) %>%
    sf::st_union() %>%
    as("Spatial")
  america_grid  <- makegrid(america_union)
  colnames(america_grid) <- c("long", "lat")
  coordinates(america_grid) <- ~ long + lat

  dat %>%
    map(function(x) { sp::coordinates(x) <- ~ long + lat; return(x) }) %>%
    ## Interpolate and rasterize dstat
    map(gstat::idw, formula = zv ~ 1, america_grid) %>%
    map_dfr(function(x) {
      rasterFromXYZ(as.data.frame(x)[,1:3]) %>%
        mask(america_union) %>%
        rasterToPoints() %>%
        as_tibble()
    }, .id = "Desc")
}