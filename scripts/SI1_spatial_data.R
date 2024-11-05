#####################################################
##                                                 ##
##         Lythrum Optimal Flowering Time          ##
##                                                 ##
##               Spatial data map                  ##
##                                                 ##
##               JJ- October 2024                  ##
##                                                 ##
#####################################################
rm(list = ls())
options(width = 100)

library(tidyverse)
library(sf)
library(rnaturalearth)

#_______________________________________________________________________________
#### 1. load data and wrangle ####

sweden_map <- ne_countries(scale = 10, returnclass = "sf", type = "map_units") %>% 
  filter(name == "Sweden")

load("data/lythrum_population_coordinates.RData")
load("data/solidago_population_coordinates.RData")

lythrum_sf <- st_as_sf(lythrum_coords, coords = c("long_dec", "lat_dec")) %>% 
  st_set_crs(st_crs(sweden_map)) %>% 
  rename(growing_season_length = capped_season_length)
solidago_sf <- st_as_sf(solidago_coordinates, coords = c("lon_dec", "lat_dec")) %>% 
  st_set_crs(st_crs(sweden_map)) %>% 
  rename(growing_season_length = capped_season_length)

umea_uppsala <- tibble(site = c("Uppsala", "Umeå"),
                       lon_dec = c(17.6292, 20.3022),
                       lat_dec = c(59.8525, 63.8190))
  # st_as_sf(coords = c("lon_dec", "lat_dec")) %>% 
  # st_set_crs(st_crs(sweden_map))

#_______________________________________________________________________________
#### 2. Map ####

mapplot <- ggplot() +
  geom_sf(data = sweden_map) +
  geom_sf(data = lythrum_sf, alpha = 0.8, size = 2.8,
          aes(colour = growing_season_length, shape = "Lythrum")) +
  geom_sf(data = solidago_sf, alpha = 0.8, size = 2.8,
          aes(colour = growing_season_length, shape = "Solidago")) +
  geom_point(data = umea_uppsala, aes(x = lon_dec, y = lat_dec), shape = 3) +
  geom_text(data = umea_uppsala, aes(x = lon_dec + c(3.3,2.5), y = lat_dec, label = site)) +
  scale_colour_viridis_c(option = "D",begin = 0.1, end = 0.9) +
  labs(shape = "Species", colour = "Growing\nSeason\nLength", 
       x = NULL, y = NULL) +
  scale_y_continuous(breaks = seq(54, 70, by = 2)) +
  scale_x_continuous(breaks = seq(12, 24, by = 4)) +
  guides(shape = guide_legend(order = 1), colour = guide_colourbar(order = 2)) +
  theme_bw() + theme(panel.grid = element_blank())

ggsave(mapplot, filename = "output/population_map.jpeg", 
       width = 10, height = 14, units = "cm", dpi = 800)  
