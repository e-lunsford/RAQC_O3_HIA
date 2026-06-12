#################################################################################
# title: Obtain BenMAP CRFs Data
# author: Beth Lunsford 
# date: 2026-06-05 
#
# This code is to establish the concentration-response functions, obtained
# from BenMAP and literature
#
# Last Run: 06/05/2026 and code was in working order using R 4.5.3
#
#################################################################################

# ------------------------------------------------------------------------------
# Set working directory to external hard drive CEBA folder.
# ------------------------------------------------------------------------------

setwd("M:/ALA_O3_HIA/RAQC_ALA_O3_HIA")

# ------------------------------------------------------------------------------
# Load Libraries
# ------------------------------------------------------------------------------
library(dplyr)
library(knitr)
library(purrr)
library(readxl)
library(stringr)
library(tidyverse)
library(units)

# ------------------------------------------------------------------------------
# CRFs from:
# Martenies SE, Akherati A, Jathar S, Magzamen S. Health and 
# Environmental Justice Implications of Retiring Two Coal-Fired Power Plants 
# in the Southern Front Range Region of Colorado. GeoHealth. 2019;3(9):266-283.
# doi:10.1029/2019GH000206
# ------------------------------------------------------------------------------

# Load 2019 paper's CRFs for O3
CRF_coal <- 
  # read file
  read_excel(path = "data/drafted_crf.xlsx", 
             sheet = 1)

glimpse(CRF_coal)

# Tidy CRFs for data merging
CRF_tidy <- CRF_coal %>%
  rename(
    health_outcome = `Outcome Name`,
    pooled_beta = `Pooled CR`,
    se_pooled_beta = `SE of Pooled CR`
  )

save(file = "output/crf_tidy.RData", CRF_tidy)



# ------------------------------------------------------------------------------
# Load in BenMAP obtained CRFs
# Load in O3 CRF raw .xlsx file and clean output to be wider table
# ------------------------------------------------------------------------------
cr_functions <- 
  
  # read file
  read_excel(path = "data/o3_CRF.xlsx",
             col_names = FALSE) %>%
  
  # Convert col. 1 into var-names
      rename(variable = ...1) %>%
  
  # Pivot Longer
    pivot_longer(
      cols = -variable,
      names_to = "study_id",
      values_to = "value"
  ) %>%
  
  # Pivot wider
    pivot_wider(
       names_from = variable,
       values_from = value
  ) %>%
  
  # convert fields to numeric
  mutate(
    Year = as.numeric(Year),
    Beta = as.numeric(Beta),
    `Std Err` = as.numeric(`Std Err`)
  ) %>%
  
  # format for future use
  rename_with(tolower) %>%
  rename(stnd_err = `std err`)


# Check results
glimpse(cr_functions)

