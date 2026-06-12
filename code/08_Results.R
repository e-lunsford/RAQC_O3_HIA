#################################################################################
# title: Summary of Results
# author: Beth Lunsford 
# date: 2026-06-06 
#
#
# Last Run: 06/09/2026 and code was in working order using R 4.5.3
#
# This code is to calculate the total burden and avoided cases
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

# Results
asthma_ed_hia_results <- read_csv(file = "output/asthma_ed_hia_results.csv")

benmap_hia_results <- read_csv(file = "output/benmap_hia_results.csv")

# Combine results
hia_results <- rbind(asthma_ed_hia_results,
                     benmap_hia_results)

# Results 2
asthma_ed_hia_results2 <- read_csv(file = "output/asthma_ed_hia_results2.csv")

benmap_hia_results2 <- read_csv(file = "output/benmap_hia_results2.csv")

hia_results2 <- rbind(asthma_ed_hia_results2,
                      benmap_hia_results2)

# Results 3
asthma_ed_hia_results3

benmap_hia_results3

hia_results3 <- rbind(asthma_ed_hia_results3,
                      benmap_hia_results3)

################################################################################
# Results by total study days
################################################################################


# Total Burden
total_burden3 <- hia_results3 %>%
  group_by(scenario, health_outcome) %>%
  summarize(
    total_delta_y = sum(mean_delta_y_study, na.rm = TRUE),
    lower_total = sum(lower_95, na.rm = TRUE),
    upper_total = sum(upper_95, na.rm = TRUE),
    total_delta_y_annual = sum(mean_delta_y_annual, na.rm = TRUE),
    lower_annual = sum(annual_lower_95),
    upper_annual = sum(annual_upper_95),
    .groups = "drop"
  )

write_csv(file = "output/total_burden3.csv",
          x = total_burden3)

# Avoided Cases
avoided_cases3 <- total_burden3 %>%
  pivot_wider(
    names_from = scenario,
    values_from = c(
      total_delta_y,
      lower_total,
      upper_total,
      total_delta_y_annual,
      lower_annual,
      upper_annual
    )
  )

avoided_cases3 <- avoided_cases3 %>%
  mutate(
    
    # Study-period SEs
    se_65 =
      (upper_total_delta_x_65 -
         lower_total_delta_x_65) / (2 * 1.96),
    
    se_70 =
      (upper_total_delta_x_70 -
         lower_total_delta_x_70) / (2 * 1.96),
    
    se_plus10 =
      (upper_total_delta_x_plus10 -
         lower_total_delta_x_plus10) / (2 * 1.96),


    additional_burden_65_vs_70 =
      total_delta_y_delta_x_65 -
      total_delta_y_delta_x_70,
    
    additional_burden_plus10_vs_70 =
      total_delta_y_delta_x_plus10 -
      total_delta_y_delta_x_70,
    
    additional_burden_plus10_vs_65 =
      total_delta_y_delta_x_plus10 -
      total_delta_y_delta_x_65,
    
    se_65_vs_70 =
      sqrt(se_65^2 + se_70^2),
    
    se_plus10_vs_70 =
      sqrt(se_plus10^2 + se_70^2),
    
    se_plus10_vs_65 =
      sqrt(se_plus10^2 + se_65^2),
    
    lower_65_vs_70 =
      additional_burden_65_vs_70 -
      1.96 * se_65_vs_70,
    
    upper_65_vs_70 =
      additional_burden_65_vs_70 +
      1.96 * se_65_vs_70,
    
    lower_plus10_vs_70 =
      additional_burden_plus10_vs_70 -
      1.96 * se_plus10_vs_70,
    
    upper_plus10_vs_70 =
      additional_burden_plus10_vs_70 +
      1.96 * se_plus10_vs_70,
    
    lower_plus10_vs_65 =
      additional_burden_plus10_vs_65 -
      1.96 * se_plus10_vs_65,
    
    upper_plus10_vs_65 =
      additional_burden_plus10_vs_65 +
      1.96 * se_plus10_vs_65,
    
    se_65_annual =
      (upper_annual_delta_x_65 -
         lower_annual_delta_x_65) / (2 * 1.96),
    
    se_70_annual =
      (upper_annual_delta_x_70 -
         lower_annual_delta_x_70) / (2 * 1.96),
    
    se_plus10_annual =
      (upper_annual_delta_x_plus10 -
         lower_annual_delta_x_plus10) / (2 * 1.96),
    
    additional_burden_65_vs_70_annual =
      total_delta_y_annual_delta_x_65 -
      total_delta_y_annual_delta_x_70,
    
    additional_burden_plus10_vs_70_annual =
      total_delta_y_annual_delta_x_plus10 -
      total_delta_y_annual_delta_x_70,
    
    additional_burden_plus10_vs_65_annual =
      total_delta_y_annual_delta_x_plus10 -
      total_delta_y_annual_delta_x_65,
    
    se_65_vs_70_annual =
      sqrt(se_65_annual^2 + se_70_annual^2),
    
    se_plus10_vs_70_annual =
      sqrt(se_plus10_annual^2 + se_70_annual^2),
    
    se_plus10_vs_65_annual =
      sqrt(se_plus10_annual^2 + se_65_annual^2),
    
    lower_65_vs_70_annual =
      additional_burden_65_vs_70_annual -
      1.96 * se_65_vs_70_annual,
    
    upper_65_vs_70_annual =
      additional_burden_65_vs_70_annual +
      1.96 * se_65_vs_70_annual,
    
    lower_plus10_vs_70_annual =
      additional_burden_plus10_vs_70_annual -
      1.96 * se_plus10_vs_70_annual,
    
    upper_plus10_vs_70_annual =
      additional_burden_plus10_vs_70_annual +
      1.96 * se_plus10_vs_70_annual,
    
    lower_plus10_vs_65_annual =
      additional_burden_plus10_vs_65_annual -
      1.96 * se_plus10_vs_65_annual,
    
    upper_plus10_vs_65_annual =
      additional_burden_plus10_vs_65_annual +
      1.96 * se_plus10_vs_65_annual
  )


avoided_cases3_final <- avoided_cases3 %>%
  dplyr::select(
    health_outcome,
    
    additional_burden_65_vs_70,
    lower_65_vs_70,
    upper_65_vs_70,
    
    additional_burden_plus10_vs_70,
    lower_plus10_vs_70,
    upper_plus10_vs_70,
    
    additional_burden_plus10_vs_65,
    lower_plus10_vs_65,
    upper_plus10_vs_65,
    
    additional_burden_65_vs_70_annual,
    lower_65_vs_70_annual,
    upper_65_vs_70_annual,
    
    additional_burden_plus10_vs_70_annual,
    lower_plus10_vs_70_annual,
    upper_plus10_vs_70_annual,
    
    additional_burden_plus10_vs_65_annual,
    lower_plus10_vs_65_annual,
    upper_plus10_vs_65_annual
  )

write_csv(
  avoided_cases3_final,
  "output/avoided_cases3.csv"
)




# ------------------------------------------------------------------------------
# Total burden by scenario
# ------------------------------------------------------------------------------
total_burden <- hia_results %>%
  group_by(scenario, health_outcome) %>%
  summarise(
    total_delta_y = sum(mean_delta_y_study, na.rm = TRUE),
    lower_total = sum(lower_95, na.rm = TRUE),
    upper_total = sum(upper_95, na.rm = TRUE),
    total_delta_y_annual = sum(mean_delta_y_annual, na.rm = TRUE),
    lower_annual = sum(annual_lower_95),
    upper_annual = sum(annual_upper_95),
    .groups = "drop"
  )

write_csv(x = total_burden,
          file = "output/total_burden.csv")


total_burden_age <- hia_results %>%
  group_by(scenario, health_outcome, age) %>%
  summarise(
    total_delta_y = sum(mean_delta_y_study, na.rm = TRUE),
    lower_total = sum(lower_95, na.rm = TRUE),
    upper_total = sum(upper_95, na.rm = TRUE),
    total_delta_y_annual = sum(mean_delta_y_annual, na.rm = TRUE),
    lower_annual = sum(annual_lower_95),
    upper_annual = sum(annual_upper_95),
    .groups = "drop"
  )

write_csv(x = total_burden_age,
          file = "output/total_burden_age.csv")



# ------------------------------------------------------------------------------
# Avoided Cases
# Note: delta_x_65 - delta_x_70 gives the additional burden 
# captured/removed under the 65 ppb scenario compared with 70 ppb.
# ------------------------------------------------------------------------------

avoided_cases <- total_burden %>%
  dplyr::select(
    health_outcome,
    scenario,
    total_delta_y,
    lower_total,
    upper_total,
    total_delta_y_annual,
    lower_annual,
    upper_annual
  ) %>%
  pivot_wider(
    names_from = scenario,
    values_from = c(
      total_delta_y,
      lower_total,
      upper_total,
      total_delta_y_annual,
      lower_annual,
      upper_annual
    )
  ) %>%
  mutate(
    # Study period avoided/additional cases
    avoided_65_vs_70_study =
      total_delta_y_delta_x_70- total_delta_y_delta_x_65,
    
    avoided_plus10_vs_70_study =
      total_delta_y_delta_x_plus10 - total_delta_y_delta_x_70,
    
    avoided_plus10_vs_65_study =
      total_delta_y_delta_x_plus10 - total_delta_y_delta_x_65,
    
    # Annual avoided/additional cases
    avoided_70_vs_65_annual =
      total_delta_y_annual_delta_x_65 - total_delta_y_annual_delta_x_70,
    
    avoided_plus10_vs_70_annual =
      total_delta_y_annual_delta_x_plus10 - total_delta_y_annual_delta_x_70,
    
    avoided_plus10_vs_65_annual =
      total_delta_y_annual_delta_x_plus10 - total_delta_y_annual_delta_x_65
  )

# Format Table
avoided_cases_long <- avoided_cases %>%
  dplyr::select(
    health_outcome,
    avoided_65_vs_70_study,
    avoided_plus10_vs_70_study,
    avoided_plus10_vs_65_study,
    avoided_65_vs_70_annual,
    avoided_plus10_vs_70_annual,
    avoided_plus10_vs_65_annual
  ) %>%
  pivot_longer(
    cols = -health_outcome,
    names_to = "comparison",
    values_to = "avoided_cases"
  )


################################################################################
# Results by exceedance days
################################################################################

# ------------------------------------------------------------------------------
# Total burden by scenario
# ------------------------------------------------------------------------------
total_burden2 <- hia_results2 %>%
  group_by(scenario, health_outcome) %>%
  summarise(
    total_delta_y = sum(mean_delta_y_study, na.rm = TRUE),
    lower_total = sum(lower_95, na.rm = TRUE),
    upper_total = sum(upper_95, na.rm = TRUE),
    total_delta_y_annual = sum(mean_delta_y_annual, na.rm = TRUE),
    lower_annual = sum(annual_lower_95),
    upper_annual = sum(annual_upper_95),
    .groups = "drop"
  )

write_csv(x = total_burden2,
          file = "output/total_burden2.csv")


total_burden_age2 <- hia_results2 %>%
  group_by(scenario, health_outcome, age) %>%
  summarise(
    total_delta_y = sum(mean_delta_y_study, na.rm = TRUE),
    lower_total = sum(lower_95, na.rm = TRUE),
    upper_total = sum(upper_95, na.rm = TRUE),
    total_delta_y_annual = sum(mean_delta_y_annual, na.rm = TRUE),
    lower_annual = sum(annual_lower_95),
    upper_annual = sum(annual_upper_95),
    .groups = "drop"
  )

write_csv(x = total_burden_age2,
          file = "output/total_burden_age2.csv")



# ------------------------------------------------------------------------------
# Avoided Cases
# Note: delta_x_65 - delta_x_70 gives the additional burden 
# captured/removed under the 65 ppb scenario compared with 70 ppb.
# ------------------------------------------------------------------------------

avoided_cases2 <- total_burden2 %>%
  dplyr::select(
    health_outcome,
    scenario,
    total_delta_y,
    lower_total,
    upper_total,
    total_delta_y_annual,
    lower_annual,
    upper_annual
  ) %>%
  pivot_wider(
    names_from = scenario,
    values_from = c(
      total_delta_y,
      lower_total,
      upper_total,
      total_delta_y_annual,
      lower_annual,
      upper_annual
    )
  ) %>%
  mutate(
    # Study period avoided/additional cases
    additional_benefit_65_vs_70_study =
      total_delta_y_delta_x_65 - total_delta_y_delta_x_70,
    
    avoided_plus10_vs_70_study =
      total_delta_y_delta_x_plus10 - total_delta_y_delta_x_70,
    
    avoided_plus10_vs_65_study =
      total_delta_y_delta_x_plus10 - total_delta_y_delta_x_65,
    
    # Annual avoided/additional cases
    avoided_65_vs_70_annual =
      total_delta_y_annual_delta_x_65 - total_delta_y_annual_delta_x_70,
    
    avoided_plus10_vs_70_annual =
      total_delta_y_annual_delta_x_plus10 - total_delta_y_annual_delta_x_70,
    
    avoided_plus10_vs_65_annual =
      total_delta_y_annual_delta_x_plus10 - total_delta_y_annual_delta_x_65
  )


write_csv(file = "output/avoided_cases2.csv",
          x = avoided_cases2)



# Format Table
avoided_cases_long2 <- avoided_cases2 %>%
  dplyr::select(
    health_outcome,
    additional_benefit_65_vs_70_study,
    avoided_plus10_vs_70_study,
    avoided_plus10_vs_65_study,
    avoided_65_vs_70_annual,
    avoided_plus10_vs_70_annual,
    avoided_plus10_vs_65_annual
  ) %>%
  pivot_longer(
    cols = -health_outcome,
    names_to = "comparison",
    values_to = "avoided_cases"
  )
























