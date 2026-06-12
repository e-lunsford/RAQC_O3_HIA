#################################################################################
# title: Setting area of interest 
# author: Beth Lunsford 
# date: 2026-06-01 
#
# This code is to set the area of interest for further spatial coding and
# map making.

# Last Run: 06/03/2026 and code was in working order using R 4.5.3
#################################################################################


# ------------------------------------------------------------------------------
# Set working directory to external hard drive folder.
# ------------------------------------------------------------------------------

setwd("M:/ALA_O3_HIA/RAQC_ALA_O3_HIA")

# ------------------------------------------------------------------------------
# Load Libraries
# ------------------------------------------------------------------------------
library(broom)
library(colorspace)
library(ggmap)
library(ggspatial)
library(ggthemes)
library(ggplot2)
library(keyring)
library(knitr)
library(leaflet)
library(lubridate)
library(purrr)
library(RAQSAPI)
library(raster)
library(sf)
library(stringr)
library(tidycensus)
library(tidyverse)
library(tigris)
library(units)



# ------------------------------------------------------------------------------
# Set parameter variables.
# ------------------------------------------------------------------------------

## Years of interest
years <- c(2020:2025)

# Colorado Boundary
colorado_sf_nad83 <- states(year = 2025, cb = T) %>%
  filter(GEOID == "08")

# Change CRS to WGS 84
colorado_sf_wgs84 <- colorado_sf_nad83 %>%
  sf::st_transform(crs = 4326)

save(file = "output/colorado_sf_wgs84.RData", x = colorado_sf_wgs84)

# Colorado bbox
colorado_bbox <- st_bbox(colorado_sf_nad83)

save(file = "ouput/colorado_bbox_nad83.RData", colorado_bbox)

# ------------------------------------------------------------------------------
# Get counties of interest
# ------------------------------------------------------------------------------

# Colorado Counties
colorado_counties_nad83 <- counties(state = "08",
                                    year = 2025,
                                    cb = T)

save(file = "output/colorado_counties_nad83.RData", x = colorado_counties_nad83)

# Transform to WGS 84
colorado_counties_wgs84 <- colorado_counties_nad83 %>%
  sf::st_transform(crs = 4326)

save(file = "output/colorado_counties_wgs84.RData", colorado_counties_wgs84)

# The Denver Metropolitan/North Front Range nonattainment region for ozone includes 
# nine counties: Adams, Arapahoe, Boulder, Broomfield, Denver, Douglas, Jefferson, 
# Weld and a portion of Larimer Counties. 
# Please note only a portion of Weld is included in the 2008 standard while all 
# of Weld is included in the 2015 standard. These counties are included because 
# they have higher ozone levels than the rest of Colorado.

counties_of_interest <- c("001", # Adams
                          "005", # Arapahoe
                          "013", # Boulder
                          "014", # Broomfield
                          "031", # Denver
                          "035", # Douglas
                          "059", # Jefferson
                          "123", # Weld
                          "069") # Larimer

DMNFR_counties_wgs84 <- colorado_counties_wgs84 %>%
  filter(COUNTYFP  %in%  counties_of_interest)

save(file = "output/DMNFR_counties_wgs84.RData", DMNFR_counties_wgs84)

# ------------------------------------------------------------------------------
# Census tract boundaries
# ------------------------------------------------------------------------------

# Colorado census tracts
colorado_ct_nad83 <- tracts(state = "08", year = 2025, cb = T)

save(file = "output/colorado_ct_nad83.RData", x = colorado_ct_nad83)

# Set CRS to WGS 84
colorado_ct_wgs_84 <- colorado_ct_nad83 %>%
  sf::st_transform(crs = 4326)

save(file = "output/colorado_ct_wgs_84.RData", colorado_ct_wgs_84)

# Denver metropolitan census tracts
dmnfr_tract_wgs84 <- colorado_ct_wgs_84 %>%
  filter(COUNTYFP  %in%  counties_of_interest)

save(file = "output/dmnfr_tract_wgs84.RData", dmnfr_tract_wgs84)



