#################################################################################
# title: Health data cleaning
# author: Beth Lunsford 
# date: 2026-06-05 
#
# This code is to obtain and format health (y0) data.
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
library(knitr)
library(lubridate)
library(purrr)
library(readxl)
library(sf)
library(sp)
library(stars)
library(stringr)
library(tibble)
library(tidycensus)
library(tidyverse)
library(units)

# ------------------------------------------------------------------------------
# Set Parameters
# ------------------------------------------------------------------------------

years_of_interest <- c(2020, 2021, 2022, 2023, 2024)

county_of_interest <- c("Adams",
                        "Arapahoe",
                        "Boulder",
                        "Broomfield",
                        "Denver",
                        "Douglas",
                        "Jefferson",
                        "Weld",
                        "Larimer")

# ------------------------------------------------------------------------------
# Load CDPHE Asthma estimates
# ------------------------------------------------------------------------------

# ED AA Asthma 2011-2024 Annual
ed_aa_asthma_annual <- read_excel(
  path = "data/EPHT_REF_COEPHT Asthma Data_2024_EN.xlsx",
  sheet = 1)

ed_aa_asthma_annual_2020_2024 <- ed_aa_asthma_annual %>%
  filter(YEAR %in% years_of_interest) %>%
  filter(stringr::str_detect(GENDER, "Both"))

statewide_ed_aa_asthma_annual_2020_2024 <- ed_aa_asthma_annual_2020_2024 %>%
  filter(stringr::str_detect(string = COUNTY,
                             pattern = "Statewide"))

dmnfr_ed_aa_asthma_annual_2020_2024 <- ed_aa_asthma_annual_2020_2024 %>%
  filter(COUNTY %in% county_of_interest)

head(dmnfr_ed_aa_asthma_annual_2020_2024)

# ------------------------------------------------------------------------------
# Overall DMNFR average regional summary
# ------------------------------------------------------------------------------

# summarized distribution of county-year rates with median, quartiles, and IQR.
dmnfr_ed_aa_asthma_summary <- dmnfr_ed_aa_asthma_annual_2020_2024 %>%
  summarize(
    weighted_rate = weighted.mean(RATE, VISITS, na.rm = TRUE),
    
    median_rate = median(RATE, na.rm = TRUE),
    q1_rate = quantile(RATE, 0.25, na.rm = TRUE),
    q3_rate = quantile(RATE, 0.75, na.rm = TRUE),
    iqr_rate = IQR(RATE, na.rm = TRUE),
    
    median_l95cl = median(L95CL, na.rm = TRUE),
    median_u95cl = median(U95CL, na.rm = TRUE),
    
    total_visits = sum(VISITS, na.rm = TRUE),
    start_year = min(YEAR, na.rm = TRUE),
    end_year = max(YEAR, na.rm = TRUE),
    
    measure = first(MEASURE)
  )

# ED AS Asthma 2011-2024 Annual
ed_as_asthma_annual <- read_excel(path = "data/EPHT_REF_COEPHT Asthma Data_2024_EN.xlsx",
                                  sheet = 3)

dmnfr_ed_as_asthma_annual_2020_2024 <- ed_as_asthma_annual %>%
  
  # Filter years 2020-2024
  filter(YEAR %in% years_of_interest) %>%
  
  # Filter both genders
  filter(stringr::str_detect(GENDER, "Both")) %>%
  
  # Filter to region
  filter(COUNTY %in% county_of_interest) %>%
  
  # Prep age column for future merging
  mutate(
    AGE = dplyr::case_when(
      str_detect(AGE, "0-4") ~ "0_4",
      str_detect(AGE, "5-14") ~ "5_14",
      str_detect(AGE, "15-34") ~ "15_34",
      str_detect(AGE, "35-64") ~ "35_64",
      str_detect(AGE, "65") ~ "65_99",
      TRUE ~ AGE
    )
  )

# HIF-ready age-specific rates
dmnfr_asthma_rates <- dmnfr_ed_as_asthma_annual_2020_2024 %>%
  filter(AGE != "All ages") %>%
  group_by(AGE) %>%
  summarize(
    y0_rate = mean(RATE, na.rm = TRUE),
    y0_mean_l95ci = mean(L95CL, na.rm = TRUE),
    y0_mean_u95ci = mean(U95CL, na.rm = TRUE),
    
    y0_median_rate = median(RATE, na.rm = TRUE),
    y0_q1_rate = quantile(RATE, 0.25, na.rm = TRUE),
    y0_q3_rate = quantile(RATE, 0.75, na.rm = TRUE),
    y0_iqr_rate = IQR(RATE, na.rm = TRUE),
    
    y0_median_l95ci = median(L95CL, na.rm = TRUE),
    y0_median_u95ci = median(U95CL, na.rm = TRUE),
    
    y0_total_visits = sum(VISITS, na.rm = TRUE),
    start_year = min(YEAR, na.rm = TRUE),
    end_year = max(YEAR, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  rename(age = AGE)


save(file = "output/dmnfr_asthma_rates.RData", dmnfr_asthma_rates)


write_csv(file = "output/dmnfr_asthma_rates.csv", x = dmnfr_asthma_rates)



# ------------------------------------------------------------------------------
# BenMAP health data
# CO is part of the West region
# ------------------------------------------------------------------------------


# Hospitalizations - respiratory from BenMap page 308
hosp_resp <- read_excel(path = "data/BenMap_Health_Rates.xlsx",
                        sheet = 4)

head(hosp_resp)

hosp_resp2 <- hosp_resp %>%
  rename(age = `Hospitalization Category`,
         y0_rate = `All Respiratory`) %>%
  mutate(health_outcome = "HOSP_resp") %>%
  dplyr::select(age, y0_rate, health_outcome)


# MRAD is 7.8 cases per person-year (BenMAP page 312)
mrad_health <- read_excel(path = "data/BenMap_Health_Rates.xlsx",
                          sheet = 6) 

mrad_health2 <- mrad_health %>%
  rename(
    health_outcome = Endpoint,
    age = Age,
    y0_rate = Rate
  ) %>%
  dplyr::select(-Parameter,
                -Source)


# School loss days is 2.2 students per year for respiratory days (BenMap 311)
# BenMAP estimates year as 180 days
sld_health <- tibble::tibble(
  health_outcome = "School_loss_day",
  age = "5_18",
  y0_rate = 2.2
)


# Non accidental morality from BenMAP page 295
non_accident_mortality <- read_excel(
  path = "data/BenMap_Health_Rates.xlsx",
  sheet = 1)

non_accident_mortality_wide <- non_accident_mortality %>%
  slice(-1) %>%
  pivot_longer(
    cols = -`Mortality Category`,
    names_to = "mort_type",
    values_to = "y0_rate"
  ) %>%
  rename(age = `Mortality Category`) %>%
  # Keep only non accidental mortality
  filter(stringr::str_detect(string = mort_type,
                             pattern = "Mortality, Non-Accidental")) %>% 
  # Relabel non-accidental to match CRF df
  mutate(health_outcome = "Non_accidental_mortality") %>%
  # Set y0 as numeric
  dplyr::mutate(y0_rate = as.numeric(y0_rate)) %>%
  # drop mort_type
  dplyr::select(-mort_type) %>%
  dplyr::mutate(age = recode(age, 
                             "Infant*" = "0_1"))
  
head(non_accident_mortality_wide)

# Combine into 1 DF
benmap_health_outcomes <- rbind(non_accident_mortality_wide,
                                hosp_resp2,
                                mrad_health2,
                                sld_health)

# Save
save(file = "output/benmap_health_outcomes.RData", benmap_health_outcomes)



# ------------------------------------------------------------------------------
# BANKED CODE from CDPHE Data
# Retain for potential future use
# ------------------------------------------------------------------------------

# Statewide ED AS Asthma
statewide_ed_as_asthma_annual_2020_2024 <- ed_as_asthma_annual_2020_2024 %>%
  filter(stringr::str_detect(string = COUNTY,
                             pattern = "Statewide"))

# Monthly ED AA Asthma 2011-2024
ed_aa_asthma_monthly <- read_excel(
  path = "data/EPHT_REF_COEPHT Asthma Data_2024_EN.xlsx",
  sheet = 2)

# ED AS Asthma 2011-2024 Annual
ed_as_asthma_annual <- read_excel(
  path = "data/EPHT_REF_COEPHT Asthma Data_2024_EN.xlsx",
  sheet = 3)

# HOSP AA Asthma 2004-2024 Annual
hosp_aa_asthma_annual <- read_excel(
  path = "data/EPHT_REF_COEPHT Asthma Data_2024_EN.xlsx",
  sheet = 4)

# HOSP AA Asthma 2004-2024 Monthly
hosp_aa_asthma_monthly <- read_excel(
  path = "data/EPHT_REF_COEPHT Asthma Data_2024_EN.xlsx",
  sheet = 5)

# HOSP AS Asthma 2004-2024 Annual
hosp_as_asthma_annual <- read_excel(
  path = "data/EPHT_REF_COEPHT Asthma Data_2024_EN.xlsx",
  sheet = 6)


