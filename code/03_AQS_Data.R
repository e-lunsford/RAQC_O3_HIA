#################################################################################
# title: Obtain EPA AQS Data
# author: Beth Lunsford 
# date: 2026-06-05 
#
# This code is to obtain annual and daily ozone monitoring station data.
#
# Last Run: 06/05/2026 and code was in working order using R 4.5.3
#
   
#################################################################################

# ------------------------------------------------------------------------------
# Set working directory to external hard drive folder.
# ------------------------------------------------------------------------------

setwd("M:/ALA_O3_HIA/RAQC_ALA_O3_HIA")

# ------------------------------------------------------------------------------
# Load Libraries
# ------------------------------------------------------------------------------
library(broom)
library(dplyr)
library(ggmap)
library(ggplot2)
library(ggspatial)
library(gstat)
library(ggthemes)
library(keyring)
library(knitr)
library(leaflet)
library(lubridate)
library(purrr)
library(RAQSAPI)
library(raster)
library(sf)
library(skimr)
library(sp)
library(spatialEco)
library(stars)
library(stringr)
library(terra)
library(tidycensus)
library(tidyverse)
library(tigris)
library(units)


# ------------------------------------------------------------------------------
# Load key-ring information to set AQS API credentials.
# ------------------------------------------------------------------------------



aqs_credentials(username = datamartAPI_user, 
                key = key_get(service = server,
                              username = datamartAPI_user))

# ------------------------------------------------------------------------------
# Load parameters 01_setting_area file.
# ------------------------------------------------------------------------------

# Years of interest
years <- c(2020:2024)

# Colorado Boundary
load(file = "output/colorado_sf_wgs84.RData")

# Colorado bbox
co_bbox <- st_bbox(colorado_sf_wgs84)

# Counties
load(file = "output/DMNFR_counties_wgs84.RData")

# Census tracts
load(file = "output/dmnfr_tract_wgs84.RData")

# DMNFR bbox
dmnfr_bbox <- st_bbox(DMNFR_counties_wgs84) + c(-0.005,-0.005,0.005,0.005)


# ------------------------------------------------------------------------------
# Obtaining EPA AQS for O3.
# ------------------------------------------------------------------------------

# Set air parameters of interest.
## Ozone - 44201


# ------------------------------------------------------------------------------
# Get the annual data.
## Obtain Ozone monitor data across Colorado for 2020-2025.
## Use the aqs annual summary function under the RAQSAPI package.
# ------------------------------------------------------------------------------

ozone_aqs_annual <- aqs_annualsummary_by_state(parameter = "44201", #O3,
                                        bdate = as.Date("2020-01-01"), #Begin Date
                                        edate = as.Date("2024-12-31"), # End Date
                                        stateFIPS = "08") # Colorado


# Save data frame as a csv file. 
write_csv(file = "data/Colorado_O3_2020_2025_annual.csv", x = ozone_aqs_annual)

ozone_aqs_annual <- read_csv(file = "data/Colorado_O3_2020_2025_annual.csv")

# ------------------------------------------------------------------------------
# Data Manipulation.
## Manipulate data into usable format. 
## - Convert data frame into an sf object.
## - Set the CRS to WGS84.
## - Mutate date into "date" format, and pull year using the lubridate package.
## - Remove some of the columns for easier future use.
# ------------------------------------------------------------------------------

ozone_aqs2_annual <- ozone_aqs_annual %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs("EPSG:4326") %>%
#  mutate(year = lubridate::year(year)) %>%
  rename("o3_mean_ppm"="arithmetic_mean",
         "o3_max_ppm" = "first_max_value") %>%
  dplyr::select(-parameter_code, -poc, -datum, -parameter,
                -sample_duration_code, -sample_duration,
                -metric_used, -method, 
                -units_of_measure, -event_type,
                -observation_count, -observation_percent,
                -validity_indicator, -required_day_count,
                -exceptional_data_count, -certification_indicator)

# Save as a csv file. 
write_csv(file = "data/ozone_aqs2_annual.csv", x = ozone_aqs2_annual)

# Save as .RData file.
save(file = "data/ozone_aqs2_annual.RData", ozone_aqs2_annual)

load(file = "data/ozone_aqs2_annual.RData")

# ------------------------------------------------------------------------------
# Use the glimpse function from dplyr to view data.
# ------------------------------------------------------------------------------

dplyr::glimpse(ozone_aqs2_annual)



# ------------------------------------------------------------------------------
# Get the daily data.
## Obtain Ozone monitor data across Colorado for 2020-2025
# #Use the aqs daily summary function under the RAQSAPI package.
# ------------------------------------------------------------------------------

ozone_aqs_daily <- aqs_dailysummary_by_state(parameter="44201",#O3
                                       bdate=as.Date("2020-01-01"),
                                       edate=as.Date("2024-12-31"),
                                       stateFIPS = "08")

# Save data frame as a csv file. 
write_csv(file = "data/Colorado_O3_2020_2025_daily.csv", x = ozone_aqs_daily)

ozone_aqs_daily <- read_csv(file = "data/Colorado_O3_2020_2025_daily.csv")

# ------------------------------------------------------------------------------
# Data Manipulation.  
## Manipulate data into usable format.
## -Convert data frame into an sf object.
## -Set the CRS to WGS84.
## -Mutate date into "date" format, and pull year using the lubridate package.
## -Remove some of the columns for easier future use.
# ------------------------------------------------------------------------------

ozone_aqs2_daily <- ozone_aqs_daily %>%
  st_as_sf(coords = c("longitude", "latitude")) %>%
  st_set_crs("EPSG:4326") %>%
  mutate(sample_date = as.Date(date_local),
         #pollutant = "ozone",
         year = lubridate::year(date_local)) %>%
  rename("o3_mean_ppm"="arithmetic_mean",
         "o3_max_ppm" = "first_max_value") %>%
  dplyr::select(-parameter_code, -poc, -datum, -parameter,
                -sample_duration_code, -sample_duration, 
                -units_of_measure, -event_type,
                -observation_count, -observation_percent,
                -validity_indicator, -first_max_hour,
                -method_code, -method)
# geometry, state_code, county_code, site_number, 
# Date, Pollutant, Concentration)

# Save as a csv file. 
write_csv(file = "data/ozone_aqs2_daily.csv", x = ozone_aqs2_daily)

# Save as .RData file.
save(file = "data/ozone_aqs2_daily.RData", ozone_aqs2_daily)

load(file= "data/ozone_aqs2_daily.RData")

# ------------------------------------------------------------------------------
# Use the glimpse function from dplyr to view data.
# ------------------------------------------------------------------------------

dplyr::glimpse(ozone_aqs2_daily)


# ------------------------------------------------------------------------------
# Calculate summary statistics by county, site number, and year.
# ------------------------------------------------------------------------------

ozone_sum_stats_annual <- ozone_aqs2_annual %>%
  st_drop_geometry() %>%
  group_by(site_number, year, site_address, local_site_name) %>%
  summarize(count = n(),
            min_av = min(o3_mean_ppm),
            max_av = max(o3_mean_ppm),
            average_av = mean(o3_mean_ppm),
            med_av = median(o3_mean_ppm),
            min_mx = min(o3_max_ppm),
            max_mx = max(o3_max_ppm),
            average_mx = mean(o3_max_ppm),
            med_mx = median(o3_max_ppm))

knitr::kable(ozone_sum_stats_annual, digits = 2)

write_csv(file = "output/o3_locations.csv", x = ozone_sum_stats)

ozone_sum_stats_daily <- ozone_aqs2_daily %>%
  st_drop_geometry() %>%
  group_by(site_number, year, site_address, local_site_name) %>%
  summarize(count = n(),
            min_av = min(o3_mean_ppm),
            max_av = max(o3_mean_ppm),
            average_av = mean(o3_mean_ppm),
            med_av = median(o3_mean_ppm),
            min_mx = min(o3_max_ppm),
            max_mx = max(o3_max_ppm),
            average_mx = mean(o3_max_ppm),
            med_mx = median(o3_max_ppm))

knitr::kable(ozone_sum_stats_daily, digits = 2)

# ------------------------------------------------------------------------------
## Check monitor locations.
# ------------------------------------------------------------------------------

ggplot() +
  geom_sf(data = dmnfr_tract_wgs84,
          inherit.aes = F,
          fill = NA,
          colour = "blue",
          linewidth = 1) +
  geom_sf(data = ozone_aqs2_annual)

# ------------------------------------------------------------------------------
# Grids, buffers, and points.
# Filter locations to those within 50 km buffer area of DMNFR
# ------------------------------------------------------------------------------

# Transform to PCS for metric units
dmnfr_proj <- st_transform(x = DMNFR_counties_wgs84,
                           crs = 102003) # USA_Contiguous_Albers_Equal_Area_Conic

# Set buffer area of dmnfr counties
dmnfr_buffer <- dmnfr_proj %>%
  st_union() %>%
  st_buffer(dist = 50000) # 50 km (meters)

# Transform back to WGSG84 to match O3 data
dmnfr_buffer_wgs <- st_transform(x = dmnfr_buffer,
                                 crs = 4326)

# Check with plot
plot(st_geometry(DMNFR_counties_wgs84))
plot(st_geometry(dmnfr_buffer_wgs), add = T, border = "green")


# Filter monitors to only include those within buffer area
ozone_dmnfr_annual <- ozone_aqs2_annual %>%
  st_filter(y = dmnfr_buffer_wgs,
            .predicate = st_intersects)

ozone_dmnfr_daily <- ozone_aqs2_daily %>%
  st_filter(y = dmnfr_buffer_wgs,
            .predicate = st_within)

# Check with plot
ggplot() +
  geom_sf(data = dmnfr_buffer_wgs,
          fill = NA,
          color = "green",
          linewidth = 1) +
  geom_sf(data = ozone_dmnfr_annual,
          color = "red",
          size = 2) +
  geom_sf(data = DMNFR_counties_wgs84,
          fill = NA,
          color = "black")



# ------------------------------------------------------------------------------
# Separate out 2008 standard and 2015 standard
# ------------------------------------------------------------------------------
years_of_interest

dmnfr_o3_daily_2008stnd <- ozone_dmnfr_daily %>%
  filter(stringr::str_detect(pollutant_standard, "2008")) %>%
  filter(year %in% years_of_interest )


dmnfr_o3_daily_2015stnd <- ozone_dmnfr_daily %>%
  filter(stringr::str_detect(pollutant_standard, "2015"))


#-------------------------------------------------------------------------------
# Calculate scenarios
# ------------------------------------------------------------------------------
dmnfr_o3_hif <- dmnfr_o3_daily_2008stnd %>%
  mutate(
    
    # Scenario 1 - reduce all exceedance to 70 ppb
    delta_x_70 = pmax(o3_max_ppm - 0.070, 0),
    
    # Scenario 2 - reduce all exceedance to 70 ppb
    delta_x_65 = pmax(o3_max_ppm - 0.065, 0),
    
    # Scenario 3 - all ozone worsens by 10 ppb everywhere
    delta_x_plus10 = 0.010
)

# ------------------------------------------------------------------------------
# Exceedance days summary information
# ------------------------------------------------------------------------------

# For HIF Calculations - summary by site days

hif_exceedance_summary <- dmnfr_o3_hif %>%
  sf::st_drop_geometry() %>%
  dplyr::select(delta_x_70, delta_x_65) %>%
  pivot_longer(
    everything(),
    names_to = "scenario",
    values_to = "exceedance"
  ) %>%
  group_by(scenario) %>%
  summarize(
    total_site_days = n(),
    exceedance_site_days = sum(exceedance > 0, na.rm = TRUE),
    site_day_ratio = exceedance_site_days / total_site_days,
    
    delta_x_min = min(exceedance, na.rm = TRUE),
    delta_x_median = median(exceedance[exceedance > 0], na.rm = TRUE),
    delta_x_q1 = quantile(exceedance[exceedance > 0], 0.25, na.rm = TRUE),
    delta_x_q3 = quantile(exceedance[exceedance > 0], 0.75, na.rm = TRUE),
    delta_x_iqr = IQR(exceedance[exceedance > 0], na.rm = TRUE),
    delta_x_max = max(exceedance, na.rm = TRUE),
    .groups = "drop"
  )

plus10_hif <- dmnfr_o3_hif %>%
  sf::st_drop_geometry() %>%
  summarize(
    scenario = "delta_x_plus10",
    total_site_days = n(),
    exceedance_site_days = sum(delta_x_plus10 > 0),
    site_day_ratio = exceedance_site_days / total_site_days,
    
    delta_x_min = min(delta_x_plus10, na.rm = TRUE),
    delta_x_median = median(delta_x_plus10),
    delta_x_q1 = quantile(delta_x_plus10, 0.25),
    delta_x_q3 = quantile(delta_x_plus10, 0.75),
    delta_x_iqr = IQR(delta_x_plus10),
    delta_x_max = max(delta_x_plus10)
  )

hif_exceedance_summary <- bind_rows(
  hif_exceedance_summary,
  plus10_hif
)


# For report of exceedances by calendar days

report_exceedance_summary <- dmnfr_o3_hif %>%
  sf::st_drop_geometry() %>%
  dplyr::select(date_local, delta_x_70, delta_x_65) %>%
  pivot_longer(
    cols = starts_with("delta_x_"),
    names_to = "scenario",
    values_to = "exceedance"
  ) %>%
  group_by(scenario, date_local) %>%
  summarize(
    any_exceedance = any(exceedance > 0),
    .groups = "drop"
  ) %>%
  group_by(scenario) %>%
  summarize(
    total_days = n_distinct(date_local),
    exceedance_days = sum(any_exceedance),
    day_ratio = exceedance_days / total_days,
    .groups = "drop"
  )

plus10_report <- dmnfr_o3_hif %>%
  sf::st_drop_geometry() %>%
  summarize(
    scenario = "delta_x_plus10",
    total_days = n_distinct(date_local),
    exceedance_days = n_distinct(date_local),
    day_ratio = 1
  )

report_exceedance_summary <- bind_rows(
  report_exceedance_summary,
  plus10_report
)



# Join as 1
o3_hif <- hif_exceedance_summary %>%
  left_join(y = report_exceedance_summary)
            
            
            
save(file = "output/o3_hif.RData", o3_hif)

write_csv(file = "output/o3_hif.csv", o3_hif)



### OLD CODE ###
# ------------------------------------------------------------------------------
# Calculate 5-year averages for ozone daily summary
# Need
# - IQR
# - average
# - Median
# ------------------------------------------------------------------------------
dmnfr_o3_hif2 <- dmnfr_o3_hif %>%
 # st_drop_geometry() %>%
  group_by(site_number, pollutant_standard, year, 
           site_address, local_site_name, county) %>%
  summarize(
    
    count = n(),
    
    # means
    min_av = min(o3_mean_ppm, na.rm = TRUE),
    max_av = max(o3_mean_ppm, na.rm = TRUE),
    average_av = mean(o3_mean_ppm, na.rm = TRUE),
    med_av = median(o3_mean_ppm, na.rm = TRUE),
    q1_av = quantile(o3_mean_ppm, 0.25, na.rm = TRUE),
    q3_av = quantile(o3_mean_ppm, 0.75, na.rm = TRUE),
    iqr_av = IQR(o3_mean_ppm, na.rm = TRUE),
    
    # MDA
    min_mx = min(o3_max_ppm, na.rm = TRUE),
    max_mx = max(o3_max_ppm, na.rm = TRUE),
    average_mx = mean(o3_max_ppm, na.rm = TRUE),
    med_mx = median(o3_max_ppm, na.rm = TRUE),
    q1_mx = quantile(o3_max_ppm, 0.25, na.rm = TRUE),
    q3_mx = quantile(o3_max_ppm, 0.75, na.rm = TRUE),
    iqr_mx = IQR(o3_max_ppm, na.rm = TRUE)
    )

knitr::kable(dmnfr_o3_hif2, digits = 2)

save(file = "output/dmnfr_o3_hif2.RData", dmnfr_o3_hif2)

write_csv(file = "output/dmnfr_o3_hif.csv", x = dmnfr_o3_hif)


dmnfr_o3_hif3 <- dmnfr_o3_hif %>%
  # st_drop_geometry() %>%
  group_by(county) %>%
  summarize(
    
    count = n(),
    
    # means
    min_av = min(o3_mean_ppm, na.rm = TRUE),
    max_av = max(o3_mean_ppm, na.rm = TRUE),
    average_av = mean(o3_mean_ppm, na.rm = TRUE),
    med_av = median(o3_mean_ppm, na.rm = TRUE),
    q1_av = quantile(o3_mean_ppm, 0.25, na.rm = TRUE),
    q3_av = quantile(o3_mean_ppm, 0.75, na.rm = TRUE),
    iqr_av = IQR(o3_mean_ppm, na.rm = TRUE),
    
    # MDA
    min_mx = min(o3_max_ppm, na.rm = TRUE),
    max_mx = max(o3_max_ppm, na.rm = TRUE),
    average_mx = mean(o3_max_ppm, na.rm = TRUE),
    med_mx = median(o3_max_ppm, na.rm = TRUE),
    q1_mx = quantile(o3_max_ppm, 0.25, na.rm = TRUE),
    q3_mx = quantile(o3_max_ppm, 0.75, na.rm = TRUE),
    iqr_mx = IQR(o3_max_ppm, na.rm = TRUE)
  )

knitr::kable(dmnfr_o3_hif3, digits = 2)


write_csv(x = dmnfr_o3_hif3,
          file = "dmnfr_o3_hif3.csv")



dmnfr_o3_hif4 <- dmnfr_o3_hif %>%
  # st_drop_geometry() %>%
  # group_by(county) %>%
  summarize(
    
    count = n(),
    
    # means
    min_av = min(o3_mean_ppm, na.rm = TRUE),
    max_av = max(o3_mean_ppm, na.rm = TRUE),
    average_av = mean(o3_mean_ppm, na.rm = TRUE),
    med_av = median(o3_mean_ppm, na.rm = TRUE),
    q1_av = quantile(o3_mean_ppm, 0.25, na.rm = TRUE),
    q3_av = quantile(o3_mean_ppm, 0.75, na.rm = TRUE),
    iqr_av = IQR(o3_mean_ppm, na.rm = TRUE),
    
    # MDA
    min_mx = min(o3_max_ppm, na.rm = TRUE),
    max_mx = max(o3_max_ppm, na.rm = TRUE),
    average_mx = mean(o3_max_ppm, na.rm = TRUE),
    med_mx = median(o3_max_ppm, na.rm = TRUE),
    q1_mx = quantile(o3_max_ppm, 0.25, na.rm = TRUE),
    q3_mx = quantile(o3_max_ppm, 0.75, na.rm = TRUE),
    iqr_mx = IQR(o3_max_ppm, na.rm = TRUE)
  )

knitr::kable(dmnfr_o3_hif4, digits = 2)


write_csv(x = dmnfr_o3_hif4,
          file = "dmnfr_o3_hif4.csv")













