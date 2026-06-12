#################################################################################
# title: Finalize Data
# author: Beth Lunsford 
# date: 2026-06-05
#
#
# Last Run: 06/05/2026 and code was in working order using R 4.5.3
#
# This code is to combine all data points into single datasets to run HIF
#                                              
################################################################################

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
library(purrr)
library(sf)
library(skimr)
library(sp)
library(stars)
library(stringr)
library(tidyverse)
library(units)


# ------------------------------------------------------------------------------
# Load parameters
# ------------------------------------------------------------------------------

# DMNFR
county_of_interest <- c("Adams",
                        "Arapahoe",
                        "Boulder",
                        "Broomfield",
                        "Denver",
                        "Douglas",
                        "Jefferson",
                        "Weld",
                        "Larimer")
# Load CRFs
load(file = "crf_tidy.RData")

# Load o3 summary
load(file = "output/o3_hif.RData")

# Population at risk
load(file = "output/dmnfr_age_pop_at_risk.RData")

# CDPHE Asthma Health
load(file = "output/dmnfr_asthma_rates.RData")

# BenMap Health
load(file = "output/benmap_health_outcomes.RData")

# ------------------------------------------------------------------------------
# Asthma ED
# ------------------------------------------------------------------------------
asthma_ed_hif_input <- dmnfr_age_pop_at_risk %>%
  # Filter for asthma ed
  filter(health_outcome == "Asthma_ED") %>%
  # Add in baseline rates (y0)
  left_join(y = dmnfr_asthma_rates,
            by = "age") %>%
  # Add in CRF
  left_join(y = CRF_tidy,
            by = "health_outcome")

save(file = "output/asthma_ed_hif_input.RData", asthma_ed_hif_input)

# load(file = output/asthma_ed_hif_input.RData)

# ------------------------------------------------------------------------------
# BenMAP health outcomes
# Respiratory hospitalizations, non-accidental mortality,
# Minor restricted activity days, school loss days
# ------------------------------------------------------------------------------

benmap_hif_input <- dmnfr_age_pop_at_risk %>%
  # Remove CDPHE data to focus on benmap health outcomes
  filter(health_outcome != "Asthma_ED") %>%
  filter(health_outcome != "Asthma_symptom_day") %>%
 
  # Join health outcome rates and CRFs
  left_join(y = benmap_health_outcomes,
            by = c("age", "health_outcome")) %>%
  left_join(y = CRF_tidy,
            by = "health_outcome") %>%
  
  # Add denominator used by each health outcome rate
  mutate(y0_per = case_when(
    health_outcome == "HOSP_resp" ~ 100,
    health_outcome == "School_loss_day" ~ 1,
    health_outcome == "Non_accidental_mortality" ~ 100,
    health_outcome == "MRAD" ~ 100000,
    TRUE ~ NA_real_
  ))

glimpse(benmap_hif_input)

save(file = "output/benmap_hif_input.RData", benmap_hif_input)





