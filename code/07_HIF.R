#################################################################################
# title: Run HIF 
# author: Beth Lunsford 
# date: 2026-06-08 
#
#
# Last Run: 06/09/2026 and code was in working order using R 4.5.3
#
# This code is to run the HIF
#
# Formula: delta-Y=Pop * Y0 * Total_Days * (1−exp)^(−beta*delta-x)
#
# Where population represents the population at risk
# Y0 represents baseline incidence
# Total days are the duration of the study (5-years; 2020-2024)
# beta are crfs obtained from benmap and literature
# delta-x represents the three different scenarios
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

load(file = "output/asthma_ed_hif_input.RData")

load(file = "output/benmap_hif_input.RData")

################################################################################
# Asthma age-specific ED (CDPHE-data) (by total study day)
################################################################################

# Formula: delta-Y=Pop * Y0 * Total_Days * (1−exp)^(−beta*delta-x)

cdphe_monte_carlo_hif <- function(asthma_ed_hif_input, # age, pop, y0, beta
                                  o3_hif, # delta-x, total days
                                  n_iter = 10000) {
  
  # -----------------------------
  # STEP 1: expand scenarios
  # -----------------------------
  o3_long <- o3_hif 
  
  # -----------------------------
  # STEP 2: create full dataset
  # (age × scenario)
  # -----------------------------
  data <- asthma_ed_hif_input %>%
    tidyr::crossing(o3_long)
  
  # -----------------------------
  # STEP 3: run MC row-wise
  # -----------------------------
  results <- lapply(seq_len(nrow(data)), function(i) {
    
    row <- data[i, ]
    
    # -------------------------
    # uncertainty: y0
    # -------------------------
    se_y0 <- (row$y0_mean_u95cl - row$y0_mean_l95cl) / (2 * 1.96)
    y0_samples <- rnorm(n = n_iter,
                        mean = row$y0_rate,
                        sd = se_y0)
    
    # -------------------------
    # uncertainty: pop
    # -------------------------
    pop_samples <- rnorm(n = n_iter, 
                         mean = row$pop, 
                         sd = row$pop_se)
    
    # -------------------------
    # uncertainty: beta
    # -------------------------
    beta_samples <- rnorm(n = n_iter, 
                          mean = row$pooled_cr, 
                          sd = row$se_pooled_cr)
    
    # -------------------------
    # scenario exposure (sample from Q1-Q3)
    # -------------------------
    delta_x_samples <- runif(n = n_iter,
                             min = row$delta_x_min,
                             max = row$delta_x_max)
    
    # -------------------------
    # Keep biologically/statistically valid values
    # -------------------------
    pop_samples <- pmax(pop_samples, 0)
    y0_samples <- pmax(y0_samples, 0)
    delta_x_samples <- pmax(delta_x_samples, 0)
    
    # -------------------------
    # Set duration
    # -------------------------
    total_days <- row$total_days
    
    # -------------------------
    # HIF equation
    # -------------------------
    delta_y <- 
      # P
      pop_samples * 
      # Y0 
      (y0_samples / 10000 / 365.25) * # Asthma ED is per 10,000 annual
      # D
      total_days * 
      # Beta and Delta X
      (1 - exp(-beta_samples * delta_x_samples))
    
    # -------------------------
    # summary output
    # -------------------------
    annual_factor <- 365.25 / total_days
    
    mean_study <- mean(delta_y)
    lower_study <- quantile(delta_y, 0.025)
    upper_study <- quantile(delta_y, 0.975)
    
    tibble::tibble(
      health_outcome = row$health_outcome,
      age = row$age,
      scenario = row$scenario,
      
      #iter = seq_len(n_iter),
      #delta_y = delta_y,
      
      mean_delta_y_study = mean_study,
      lower_95 = lower_study,
      upper_95 = upper_study,
      
      mean_delta_y_annual = mean_study * annual_factor,
      annual_lower_95 = lower_study * annual_factor,
      annual_upper_95 = upper_study * annual_factor
    )
  })
  
  dplyr::bind_rows(results)
}
# ------------------------------------------------------------------------------
# Run the model
# ------------------------------------------------------------------------------
asthma_ed_hia_results3 <- cdphe_monte_carlo_hif(
  asthma_ed_hif_input,
  o3_hif
)

write_csv(file = "output/asthma_ed_hia_results3.csv",
          x = asthma_ed_hia_results)



################################################################################
# BenMAP health outcomes
# Respiratory hospitalizations, non-accidental mortality,
# Minor restricted activity days, school loss days
# By total study days
################################################################################

# ------------------------------------------------------------------------------
# Formula: delta-Y = Pop * Y0 * Total_Days * (1−exp) ^ (−beta * delta-x)
# ------------------------------------------------------------------------------

benmap_monte_carlo_hif <- function(benmap_hif_input, # age, pop, y0, beta
                                   o3_hif, # delta-x, total days
                                   n_iter = 10000) {
  
  # -----------------------------
  # STEP 1: expand scenarios
  # -----------------------------
  o3_long <- o3_hif %>%
    dplyr::select(scenario, delta_x_median, total_days,
                  delta_x_min, delta_x_max)
  
  # -----------------------------
  # STEP 2: create full dataset
  # (age × scenario)
  # -----------------------------
  data <- benmap_hif_input %>%
    tidyr::crossing(o3_long)
  
  # -----------------------------
  # STEP 3: run MC row-wise
  # -----------------------------
  results <- lapply(seq_len(nrow(data)), function(i) {
    
    row <- data[i, ]
    
    # -------------------------
    # uncertainty: y0
    # -------------------------
    #se_y0 <- (row$y0_median_u95cl - row$y0_median_l95cl) / (2 * 1.96)
    y0_samples <- rnorm(n = n_iter, 
                        mean = row$y0_rate)
    
    # -------------------------
    # uncertainty: pop
    # -------------------------
    pop_samples <- rnorm(n = n_iter, 
                         mean = row$pop, 
                         sd = row$pop_se)
    
    # -------------------------
    # uncertainty: beta
    # -------------------------
    beta_samples <- rnorm(n = n_iter,
                          mean = row$pooled_cr,
                          sd = row$se_pooled_cr)
    
    # -------------------------
    # scenario exposure
    # -------------------------
    delta_x_samples <- runif(n = n_iter,
                             min = row$delta_x_min,
                             max = row$delta_x_max)
    
    # -------------------------
    # Keep biologically/statistically valid values
    # -------------------------
    pop_samples <- pmax(pop_samples, 0)
    y0_samples <- pmax(y0_samples, 0)
    delta_x_samples <- pmax(delta_x_samples, 0)
    
    # -------------------------
    # Set duration
    # -------------------------
    total_days <- row$total_days
    y0_per <- row$y0_per
    
    # -------------------------
    # HIF equation
    # -------------------------
    delta_y <- pop_samples *
      (y0_samples / y0_per / 365.25) * 
      total_days *
      (1 - exp(-beta_samples * delta_x_samples))
    
    # -------------------------
    # summary output
    # -------------------------
    annual_factor <- 365.25 / total_days
    
    mean_study <- mean(delta_y)
    lower_study <- quantile(delta_y, 0.025)
    upper_study <- quantile(delta_y, 0.975)
    
    tibble::tibble(
      health_outcome = row$health_outcome,
      age = row$age,
      scenario = row$scenario,
      
      #iter = seq_len(n_iter),
     # delta_y = delta_y,
      
      mean_delta_y_study = mean_study,
      lower_95 = lower_study,
      upper_95 = upper_study,
      
      mean_delta_y_annual = mean_study * annual_factor,
      annual_lower_95 = lower_study * annual_factor,
      annual_upper_95 = upper_study * annual_factor
    )
  })
  
  dplyr::bind_rows(results)
}
# -------------------------
# Run the model
# -------------------------

benmap_hia_results3 <- benmap_monte_carlo_hif(
  benmap_hif_input = benmap_hif_input,
  o3_hif = o3_hif,
  n_iter = 10000
)

write_csv(file = "output/benmap_hia_results.csv",
          x = benmap_hia_results)






################################################################################
# Asthma age-specific ED (CDPHE-data) (by total exceedance days)
################################################################################

# Formula: delta-Y=Pop * Y0 * Total_Days * (1−exp)^(−beta*delta-x)

cdphe_monte_carlo_hif2 <- function(asthma_ed_hif_input, # age, pop, y0, beta
                                  o3_hif, # delta-x, total days
                                  n_iter = 10000) {
  
  # -----------------------------
  # STEP 1: expand scenarios
  # -----------------------------
  o3_long <- o3_hif
  
  # -----------------------------
  # STEP 2: create full dataset
  # (age × scenario)
  # -----------------------------
  data <- asthma_ed_hif_input %>%
    tidyr::crossing(o3_long)
  
  # -----------------------------
  # STEP 3: run MC row-wise
  # -----------------------------
  results <- lapply(seq_len(nrow(data)), function(i) {
    
    row <- data[i, ]
    
    # -------------------------
    # uncertainty: y0
    # -------------------------
    se_y0 <- (row$y0_mean_u95cl - row$y0_mean_l95cl) / (2 * 1.96)
    y0_samples <- rnorm(n = n_iter,
                        mean = row$y0_rate,
                        sd = se_y0)
    
    # -------------------------
    # uncertainty: pop
    # -------------------------
    pop_samples <- rnorm(n = n_iter, 
                         mean = row$pop, 
                         sd = row$pop_se)
    
    # -------------------------
    # uncertainty: beta
    # -------------------------
    beta_samples <- rnorm(n = n_iter, 
                          mean = row$pooled_cr, 
                          sd = row$se_pooled_cr)
    
    # -------------------------
    # scenario exposure (sample from Q1-Q3)
    # -------------------------
    delta_x_samples <- runif(n = n_iter,
                             min = row$delta_x_min,
                             max = row$delta_x_max)
    
    # -------------------------
    # Keep biologically/statistically valid values
    # -------------------------
    pop_samples <- pmax(pop_samples, 0)
    y0_samples <- pmax(y0_samples, 0)
    delta_x_samples <- pmax(delta_x_samples, 0)
    
    # -------------------------
    # Set duration
    # -------------------------
    total_days <- row$total_days
    duration_days <- ifelse(
      row$scenario %in% c("delta_x_65", "delta_x_70"),
      row$exceedance_days,
      row$total_days
    )
    
    # -------------------------
    # HIF equation
    # -------------------------
    delta_y <- 
      # P
      pop_samples * 
      # Y0 
      (y0_samples / 10000 / 365.25) * # Asthma ED is per 10,000 annual
      # D
      duration_days * 
      # Beta and Delta X
      (1 - exp(-beta_samples * delta_x_samples))
    
    # -------------------------
    # summary output
    # -------------------------
    annual_factor <- 365.25 / duration_days
    
    mean_study <- mean(delta_y)
    lower_study <- quantile(delta_y, 0.025)
    upper_study <- quantile(delta_y, 0.975)
    
    tibble::tibble(
      health_outcome = row$health_outcome,
      age = row$age,
      scenario = row$scenario,
      duration_days = duration_days,
      
      mean_delta_y_study = mean_study,
      lower_95 = lower_study,
      upper_95 = upper_study,
      
      mean_delta_y_annual = mean_study * annual_factor,
      annual_lower_95 = lower_study * annual_factor,
      annual_upper_95 = upper_study * annual_factor
    )
  })
  
  dplyr::bind_rows(results)
}
# ------------------------------------------------------------------------------
# Run the model
# ------------------------------------------------------------------------------
asthma_ed_hia_results2 <- cdphe_monte_carlo_hif2(
  asthma_ed_hif_input,
  o3_hif
)

write_csv(file = "output/asthma_ed_hia_results2.csv",
          x = asthma_ed_hia_results2)



################################################################################
# BenMAP health outcomes
# Respiratory hospitalizations, non-accidental mortality,
# Minor restricted activity days, school loss days
# By total study days
################################################################################

# ------------------------------------------------------------------------------
# Formula: delta-Y = Pop * Y0 * Total_Days * (1−exp) ^ (−beta * delta-x)
# ------------------------------------------------------------------------------

benmap_monte_carlo_hif2 <- function(benmap_hif_input, # age, pop, y0, beta
                                   o3_hif, # delta-x, total days
                                   n_iter = 10000) {
  
  # -----------------------------
  # STEP 1: expand scenarios
  # -----------------------------
  o3_long <- o3_hif 
  
  # -----------------------------
  # STEP 2: create full dataset
  # (age × scenario)
  # -----------------------------
  data <- pop_at_risk2 %>%
    tidyr::crossing(o3_long)
  
  # -----------------------------
  # STEP 3: run MC row-wise
  # -----------------------------
  results <- lapply(seq_len(nrow(data)), function(i) {
    
    row <- data[i, ]
    
    # -------------------------
    # uncertainty: y0
    # -------------------------
    #se_y0 <- (row$y0_median_u95cl - row$y0_median_l95cl) / (2 * 1.96)
    y0_samples <- rnorm(n = n_iter, 
                        mean = row$y0_rate)
    
    # -------------------------
    # uncertainty: pop
    # -------------------------
    pop_samples <- rnorm(n = n_iter, 
                         mean = row$pop, 
                         sd = row$pop_se)
    
    # -------------------------
    # uncertainty: beta
    # -------------------------
    beta_samples <- rnorm(n = n_iter,
                          mean = row$pooled_cr,
                          sd = row$se_pooled_cr)
    
    # -------------------------
    # scenario exposure
    # -------------------------
    delta_x_samples <- runif(n = n_iter,
                             min = row$delta_x_min,
                             max = row$delta_x_max)
    
    # -------------------------
    # Keep biologically/statistically valid values
    # -------------------------
    pop_samples <- pmax(pop_samples, 0)
    y0_samples <- pmax(y0_samples, 0)
    delta_x_samples <- pmax(delta_x_samples, 0)
    
    # -------------------------
    # Set duration
    # -------------------------
    total_days <- row$total_days
    duration_days <- ifelse(
      row$scenario %in% c("delta_x_65", "delta_x_70"),
      row$exceedance_days,
      row$total_days
    )
    y0_per <- row$y0_per
    
    # -------------------------
    # HIF equation
    # -------------------------
    delta_y <- pop_samples *
      (y0_samples / y0_per / 365.25) * 
      duration_days *
      (1 - exp(-beta_samples * delta_x_samples))
    
    # -------------------------
    # summary output
    # -------------------------
    annual_factor <- 365.25 / total_days
    
    mean_study <- mean(delta_y)
    lower_study <- quantile(delta_y, 0.025)
    upper_study <- quantile(delta_y, 0.975)
    
    tibble::tibble(
      health_outcome = row$health_outcome,
      age = row$age,
      scenario = row$scenario,
      duration_days = duration_days,
      
      
      mean_delta_y_study = mean_study,
      lower_95 = lower_study,
      upper_95 = upper_study,
      
      mean_delta_y_annual = mean_study * annual_factor,
      annual_lower_95 = lower_study * annual_factor,
      annual_upper_95 = upper_study * annual_factor
    )
  })
  
  dplyr::bind_rows(results)
}
# -------------------------
# Run the model
# -------------------------

benmap_hia_results2 <- benmap_monte_carlo_hif2(
  benmap_hif_input = benmap_hif_input,
  o3_hif = o3_hif,
  n_iter = 10000
)

write_csv(file = "output/benmap_hia_results2.csv",
          x = benmap_hia_results2)














