################################################################################
# title: Obtain ACS Data 
# author: Beth Lunsford 
# date: 2026-06-08 
#
#
# Last Run: 06/08/2026 and code was in working order using R 4.5.3
#
#
# This code is to obtain population data from ACS and ALA.
# ACS 5 year estimates for 2020-2024 for total population and under 18.
# 
# https://www2.census.gov/programs-surveys/acs/tech_docs/statistical_testing/2023_Instructions_for_Stat_Testing_ACS.pdf
#   
# SE for each population age group is derived from MOE with Z = 1.645.
# Standard Error = Margin of Error / 1.645 
# MOE(A+B) = SQRT(MOE_A + MOE_B)
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
# Population at risk data from ALA 2026 Report card for Colorado. Downloaded
# as a csv file.
# ------------------------------------------------------------------------------

pop_at_risk <- readr::read_csv(file = "data/ALA_pop_at_risk.csv",
                               col_names = TRUE)

county_of_interest <- c("Adams",
                        "Arapahoe",
                        "Boulder",
                        "Broomfield",
                        "Denver",
                        "Douglas",
                        "Jefferson",
                        "Weld",
                        "Larimer")

DMNFR_pop_at_risk <- pop_at_risk %>%
  filter(County %in% county_of_interest)

#-------------------------------------------------------------------------------
# Obtain ACS population data.
#-------------------------------------------------------------------------------

# Extract ACS 
key = "47f2731165560967710bbed4138568c692bcb565"
census_api_key(key = key, install = TRUE, overwrite = TRUE)
# Load 2025 ACS5 variable list
acs5_var_list <- load_variables(2024, "acs5", cache = TRUE)

# Create new list with obtained codes
acs5_vars <- c(total_pop = "B01001_001",
               total_male = "B01001_002",
               male_under5 = "B01001_003",
               male_5_9 = "B01001_004",
               male_10_14 = "B01001_005",
               male_15_17 = "B01001_006",
               male_18_19 = "B01001_007",
               male_20 = "B01001_008",
               male_21 = "B01001_009",
               male_22_24 = "B01001_010",
               male_25_29 = "B01001_011",
               male_30_34 = "B01001_012",
               male_35_39 = "B01001_013",
               male_40_44 = "B01001_014",
               male_45_49 = "B01001_015",
               male_50_54 = "B01001_016",
               male_55_59 = "B01001_017",
               male_60_61 = "B01001_018",
               male_62_64 = "B01001_019",
               male_65_66 = "B01001_020",
               male_67_69 = "B01001_021",
               male_70_74 = "B01001_022",
               male_75_79 = "B01001_023",
               male_80_84 = "B01001_024",
               male_85 = "B01001_025",
               total_female = "B01001_026",
               female_under5 = "B01001_027",
               female_5_9 = "B01001_028",
               female_10_14 = "B01001_029",
               female_15_17 = "B01001_030",
               female_18_19 = "B01001_031",
               female_20 = "B01001_032",
               female_21 = "B01001_033",
               female_22_24 = "B01001_034",
               female_25_29 = "B01001_035",
               female_30_34 = "B01001_036",
               female_35_39 = "B01001_037",
               female_40_44 = "B01001_038",
               female_45_49 = "B01001_039",
               female_50_54 = "B01001_040",
               female_55_59 = "B01001_041",
               female_60_61 = "B01001_042",
               female_62_64 = "B01001_043",
               female_65_66 = "B01001_044",
               female_67_69 = "B01001_045",
               female_70_74 = "B01001_046",
               female_75_79 = "B01001_047",
               female_80_84 = "B01001_048",
               female_85 = "B01001_049",
               nonhisp_black_pop = "B03002_004",
               nonhisp_white_pop = "B03002_003",
               hisp_pop = "B03002_012")

# Get Denver ACS at the census 
acs5_2024 <- get_acs(geography = "county", #county
                         variables = acs5_vars, # variable list
                         year = 2024, #2020-2024 5 yr acs
                         #output = "tidy",
                         state = "Colorado",
                         geometry = TRUE,
                         key = key,
                         survey = "acs5") #default of get acs is 5

class(acs5_2024)
save(file = "data/acs5_2024.RData", acs5_2024)

rm(acs5_2024)

load(file = "data/acs5_2024.RData")

##___ save as CSV
#readr::write_csv(file = "data/acs5_2024.csv", x = acs5_2024)
# acs5_2024 <- readr::read_csv(file = "data/acs5_2024.csv")

#-------------------------------------------------------------------------------
# Match coordinate reference systems between ACS and air pollutant estimates
## EPA air data CRS is WGS84, therefore the exposure estimates are in WGS84.
#-------------------------------------------------------------------------------

acs5_2024_wgs <- acs5_2024 %>%
  # Use st_transform to change ACS from NAD83 to WGS84
    st_transform(crs = 4326) %>%
  
  # Pivot wider for age summaries.
  pivot_wider(id_cols = c(GEOID, NAME, geometry),
              names_from = variable,
              values_from = c(estimate, moe))

################################################################################
# Calculate total count and SE (male and female) by age.
################################################################################

Z = 1.645

acs_pop <- acs5_2024_wgs %>%
  mutate(
    # Total population
    total_pop = estimate_total_male + estimate_total_female,
    total_pop_se = sqrt(moe_total_male + moe_total_female) / Z,
   
    # Under 5 population
    under5 = estimate_male_under5 + estimate_female_under5,
    under5_se = sqrt(moe_male_under5 + moe_female_under5) / Z,
    
    # ages 5 to 9 population
    age_5_9 = estimate_male_5_9 + estimate_female_5_9,
    age_5_9_se = sqrt(moe_male_5_9 + moe_female_5_9) / Z,
   
    # ages 10 to 14 population
    age_10_14 = estimate_male_10_14 + estimate_female_10_14,
    age_10_14_se = sqrt(moe_male_10_14 + moe_female_10_14) / Z,
  
    # ages 15 to 17 population
    age_15_17 = estimate_male_15_17 + estimate_female_15_17,
    age_15_17_se = sqrt(moe_male_15_17 + moe_female_15_17) / Z,
  
    # age 18 and 19 population
    age_18_19 = estimate_male_18_19 + estimate_female_18_19,
    age_18_19_se = sqrt(moe_male_18_19 + moe_female_18_19) / Z,
  
    # age 20 population
    age_20 = estimate_male_20 + estimate_female_20,
    age_20_se = sqrt(moe_male_20 + moe_female_20) / Z,
   
    # age 21
    age_21 = estimate_male_21 + estimate_female_21,
    age_21_se = sqrt(moe_male_21 + moe_female_21) / Z,
    
    # age 22 to 24
    age_22_24 = estimate_male_22_24 + estimate_female_22_24,
    age_22_24_se = sqrt(moe_male_22_24 + moe_female_22_24) / Z,
    
    # age 25 to 29
    age_25_29 = estimate_male_25_29 + estimate_female_25_29,
    age_25_29_se = sqrt(moe_male_25_29 + moe_female_25_29) / Z,
    
    # age 30 to 34
    age_30_34 = estimate_male_30_34 + estimate_female_30_34,
    age_30_34_se = sqrt(moe_male_30_34 + moe_female_30_34) / Z,
    
    # age 35 to 39
    age_35_39 = estimate_male_35_39 + estimate_female_35_39,
    age_35_39_se = sqrt(moe_male_35_39 + moe_female_35_39) / Z,
    
    # age 40 to 44
    age_40_44 = estimate_male_40_44 + estimate_female_40_44,
    age_40_44_se = sqrt(moe_male_40_44 + moe_female_40_44) / Z,
    
    # age 45 to 49
    age_45_49 = estimate_male_45_49 + estimate_female_45_49,
    age_45_49_se = sqrt(moe_male_45_49 + moe_female_45_49) / Z,
    
    # age 50 to 54
    age_50_54 = estimate_male_50_54 + estimate_female_50_54,
    age_50_54_se = sqrt(moe_male_50_54 + moe_female_50_54) / Z,
    
    # age 55 to 59
    age_55_59 = estimate_male_55_59 + estimate_female_55_59,
    age_55_59_se = sqrt(moe_male_55_59 + moe_female_55_59) / Z,
    
    # age 60 to 61
    age_60_61 = estimate_male_60_61 + estimate_female_60_61,
    age_60_61_se = sqrt(moe_male_60_61 + moe_female_60_61) / Z,
    
    # age 62 to 64
    age_62_64 = estimate_male_62_64 + estimate_female_62_64,
    age_62_64_se = sqrt(moe_male_62_64 + moe_female_62_64) / Z,
   
    # age 65 to 66
    age_65_66 = estimate_male_65_66 + estimate_female_65_66,
    age_65_66_se = sqrt(moe_male_65_66 + moe_female_65_66) / Z,
 
    # age 67 to 69
    age_67_69 = estimate_male_67_69 + estimate_female_67_69,
    age_67_69_se = sqrt(moe_male_67_69 + moe_female_67_69) / Z,
   
    # age 70 to 74
    age_70_74 = estimate_male_70_74 + estimate_female_70_74,
    age_70_74_se = sqrt(moe_male_70_74 + moe_female_70_74) / Z,
  
    # age 75 to 79
    age_75_79 = estimate_male_75_79 + estimate_female_75_79,
    age_75_79_se = sqrt(moe_male_75_79 + moe_female_75_79) / Z,
  
    # age 80 to 84
    age_80_84 = estimate_male_80_84 + estimate_female_80_84,
    age_80_84_se = sqrt(moe_male_80_84 + moe_female_80_84) / Z,
  
    # age 85 plus
    age_85 = estimate_male_85 + estimate_female_85,
    age_85_se = sqrt(moe_male_85 + moe_female_85) / Z
    
  )

################################################################################
# Assuming equal distribution across all ages
################################################################################

# calculate age range stats to match BenMap or CDPHE

#-------------------------------------------------------------------------------
# CDPHE Asthma ED
#-------------------------------------------------------------------------------
asthma_ed_pop_cdphe <- acs_pop %>%
   mutate(
    
    # 0-4
    pop_0_4 = under5, 
    pop_0_4_se = under5_se,
    
    # Age 5-14
    pop_5_14 = age_5_9 + age_10_14,
    pop_5_14_se = sqrt(age_5_9_se^2 + age_10_14_se^2),
    
    # Age 15-34
    pop_15_34 = age_15_17 + age_18_19 + age_20 + 
      age_21 + age_22_24 + age_25_29 + age_30_34,
    pop_15_34_se = sqrt(age_15_17_se^2 + age_18_19_se^2 + age_20_se^2 +
                          age_21_se^2 + age_22_24_se^2 +
                          age_25_29_se^2 + age_30_34_se^2),
    # 35-64
    pop_35_64 = age_35_39 + age_40_44 + 
      age_45_49 + age_50_54 +
      age_55_59 + age_60_61 + age_62_64,
    pop_35_64_se = sqrt(age_35_39_se^2 + age_40_44_se^2 +
                          age_45_49_se^2 + age_50_54_se^2 +
                          age_55_59_se^2 + age_60_61_se^2 + age_62_64_se^2),
    
    
    # 65+
    pop_65_99 = age_65_66 + age_67_69 + age_70_74 + 
      age_75_79 + age_80_84 +
      age_85,
    pop_65_99_se = sqrt(age_65_66_se^2 + age_67_69_se^2 + age_70_74_se^2 +
                          age_75_79_se^2 + age_80_84_se^2 +
                          age_85_se^2)
  )

save(file = "output/asthma_ed_pop_cdphe.RData", asthma_ed_pop_cdphe)

#-------------------------------------------------------------------------------
# Hospital - respiratory; benmap ages
#-------------------------------------------------------------------------------

HOSP_resp_pop_benmap <- acs_pop %>%
  mutate(
    
    pop_0_1 = under5 / 5,
    pop_0_1_se = under5_se / 5,
    
    # Age 2-17
    pop_2_17 = (0.75*under5) + age_5_9 + age_10_14 + age_15_17,
    pop_2_17_se = sqrt(((0.75*under5_se)^2) + age_5_9_se^2 + age_10_14_se^2 +
                         age_15_17_se^2),
    
    # 18-24
    pop_18_24 = age_18_19 + age_20 + age_21 + age_22_24,
    pop_18_24_se = sqrt(age_18_19_se^2 + age_20_se^2 +
                          age_21_se^2 + age_22_24_se^2),
    # 25-34
    pop_25_34 = age_25_29 + age_30_34,
    pop_25_34_se = sqrt(age_25_29_se^2 + age_30_34_se^2),
    
    # 35-44
    pop_35_44 = age_35_39 + age_40_44,
    pop_35_44_se = sqrt(age_35_39_se^2 + age_40_44_se^2),
    
    # 45-54
    pop_45_54 = age_45_49 + age_50_54,
    pop_45_54_se = sqrt(age_45_49_se^2 + age_50_54_se^2),
    
    # 55-64
    pop_55_64 = age_55_59 + age_60_61 + age_62_64,
    pop_55_64_se = sqrt(age_55_59_se^2 + age_60_61_se^2 + age_62_64_se^2),
    
    # 65-74
    pop_65_74 = age_65_66 + age_67_69 + age_70_74,
    pop_65_74_se = sqrt(age_65_66_se^2 + age_67_69_se^2 + age_70_74_se^2),
    
    # 75-84
    pop_75_84 = age_75_79 + age_80_84,
    pop_75_84_se = sqrt(age_75_79_se^2 + age_80_84_se^2),
    
    # 85+
    pop_85_99 = age_85,
    pop_85_99_se = age_85_se
  )

save(file = "output/HOSP_resp_pop_benmap.RData", x = HOSP_resp_pop_benmap)

#-------------------------------------------------------------------------------
# Asthma symptom days
#-------------------------------------------------------------------------------

asthma_symp_pop_benmap <- acs_pop %>%
  mutate(
    
    # Age 5-12
    pop_5_12 = age_5_9 + (age_10_14*0.6),
    pop_5_12_se = sqrt(age_5_9_se^2 + ((age_10_14_se*0.6)^2)),
    
    # 6-17
    pop_6_17 = (age_5_9*0.8) + age_10_14 + age_15_17,
    pop_6_17_se = sqrt(((age_5_9_se*0.8)^2) + age_10_14_se^2 +
                         age_15_17_se^2),
    
    # 9-11
    pop_9_11 = (age_5_9*0.2) + (age_10_14*0.4),
    pop_9_11_se = sqrt(((age_5_9_se*0.2)^2) + ((age_10_14_se*0.4)^2))
  )

save(file = "output/asthma_symp_pop_benmap.RData", asthma_symp_pop_benmap)

#-------------------------------------------------------------------------------
# Minor restricted activity days
#-------------------------------------------------------------------------------
MRAD_pop_benmap <- acs_pop %>%
  mutate(
    
    # Age 18-64
    pop_18_64 = age_18_19 + age_20 + age_21 + 
      age_22_24 + age_25_29 + age_30_34 +
      age_35_39 + age_40_44 + age_45_49 + 
      age_50_54 + age_55_59 + age_60_61 + 
      age_62_64,
    pop_18_64_se = sqrt(age_18_19_se^2 + age_20_se^2 +
                          age_21_se^2 + age_22_24_se^2 +
                          age_25_29_se^2 + age_30_34_se^2 +
                          age_35_39_se^2 + age_40_44_se^2 +
                          age_45_49_se^2 +
                          age_50_54_se^2 + age_55_59_se^2 + 
                          age_60_61_se^2 + age_62_64_se^2)
    
  )


save(file = "output/MRAD_pop_benmap.RData", MRAD_pop_benmap)

#-------------------------------------------------------------------------------
# School Loss Days - BenMap
## Assume school-age of 6-18
#-------------------------------------------------------------------------------
SLD_pop_benmap <- acs_pop %>%
  mutate(
    
    # Age 6-18
    pop_5_18 = (age_5_9) + age_10_14 + age_15_17 + (0.5*age_18_19),
    pop_5_18_se = sqrt(age_5_9_se^2 + age_10_14_se^2 +
                         age_15_17_se^2 +
                         ((0.5 * age_18_19_se)^2))
  )

save(file = "output/SLD_pop_benmap.RData", SLD_pop_benmap)

#-------------------------------------------------------------------------------
# Non-accidental mortality
#-------------------------------------------------------------------------------
nonacc_mortality_pop_benmap <- acs_pop %>%
  mutate(
    
    pop_0_1 = under5 / 5,
    pop_0_1_se = under5_se / 5,
    
    # Age 1-17
    pop_1_17 = (0.8*under5) + age_5_9 + age_10_14 + age_15_17,
    pop_1_17_se = sqrt(((0.8*under5_se)^2) + age_5_9_se^2 + age_10_14_se^2 +
                         age_15_17_se^2),
    
    # 18-24
    pop_18_24 = age_18_19 + age_20 + age_21 + age_22_24,
    pop_18_24_se = sqrt(age_18_19_se^2 + age_20_se^2 +
                          age_21_se^2 + age_22_24_se^2),
    # 25-34
    pop_25_34 = age_25_29 + age_30_34,
    pop_25_34_se = sqrt(age_25_29_se^2 + age_30_34_se^2),
    
    # 35-44
    pop_35_44 = age_35_39 + age_40_44,
    pop_35_44_se = sqrt(age_35_39_se^2 + age_40_44_se^2),
    
    # 45-54
    pop_45_54 = age_45_49 + age_50_54,
    pop_45_54_se = sqrt(age_45_49_se^2 + age_50_54_se^2),
    
    # 55-64
    pop_55_64 = age_55_59 + age_60_61 + age_62_64,
    pop_55_64_se = sqrt(age_55_59_se^2 + age_60_61_se^2 + age_62_64_se^2),
    
    # 65-74
    pop_65_74 = age_65_66 + age_67_69 + age_70_74,
    pop_65_74_se = sqrt(age_65_66_se^2 + age_67_69_se^2 + age_70_74_se^2),
    
    # 75-84
    pop_75_84 = age_75_79 + age_80_84,
    pop_75_84_se = sqrt(age_75_79_se^2 + age_80_84_se^2),
    
    # 85+
    pop_85_99 = age_85,
    pop_85_99_se = age_85_se
    
  )

save(file = "output/nonacc_mortality_pop_benmap.RData", nonacc_mortality_pop_benmap)



################################################################################
# Format for future merging
################################################################################

# regional age table ready for rate application
library(dplyr)
library(tidyr)
library(sf)
library(stringr)

#-------------------------------------------------------------------------------
# Asthma ED
#-------------------------------------------------------------------------------
load(file = "output/asthma_ed_pop_cdphe.RData")

head(asthma_ed_pop_cdphe)

# 1. regional population + SE (correct variance rule)
pop_region <- asthma_ed_pop_cdphe %>%
  dplyr::filter(str_detect(NAME, paste(county_of_interest, collapse = "|"))) %>%
  sf::st_drop_geometry() %>%
  dplyr::summarise(
    across(
      starts_with("pop_") & !ends_with("_se"),
      ~sum(.x, na.rm = TRUE)
    ),
    across(
      ends_with("_se"),
      ~sqrt(sum(.x^2, na.rm = TRUE))
    )
  )

# 2. pivot POP only
pop_long <- pop_region %>%
  dplyr::select(starts_with("pop_"), -ends_with("_se")) %>%
  pivot_longer(
    cols = everything(),
    names_to = "age",
    values_to = "pop"
  ) %>%
  mutate(age = str_remove(age, "pop_"))

# 3. pivot SE only
se_long <- pop_region %>%
  dplyr::select(ends_with("_se")) %>%
  pivot_longer(
    cols = everything(),
    names_to = "age",
    values_to = "se"
  ) %>%
  mutate(age = str_remove(age, "pop_") %>% str_remove("_se"))

# 4. combine
dmnfr_age_pop_asthma_ed <- left_join(pop_long, se_long, by = "age")

save(file = "output/dmnfr_age_pop_asthma_ed.RData", dmnfr_age_pop_asthma_ed)


#-------------------------------------------------------------------------------
# Hospital - respiratory; benmap ages
#-------------------------------------------------------------------------------

# 1. regional population + SE (correct variance rule)
pop_region <- HOSP_resp_pop_benmap %>%
  dplyr::filter(str_detect(NAME, paste(county_of_interest, collapse = "|"))) %>%
  sf::st_drop_geometry() %>%
  dplyr::summarise(
    across(
      starts_with("pop_") & !ends_with("_se"),
      ~sum(.x, na.rm = TRUE)
    ),
    across(
      ends_with("_se"),
      ~sqrt(sum(.x^2, na.rm = TRUE))
    )
  )

# 2. pivot POP only
pop_long <- pop_region %>%
  dplyr::select(starts_with("pop_"), -ends_with("_se")) %>%
  pivot_longer(
    cols = everything(),
    names_to = "age",
    values_to = "pop"
  ) %>%
  mutate(age = str_remove(age, "pop_"))

# 3. pivot SE only
se_long <- pop_region %>%
  dplyr::select(ends_with("_se")) %>%
  pivot_longer(
    cols = everything(),
    names_to = "age",
    values_to = "se"
  ) %>%
  mutate(age = str_remove(age, "pop_") %>% str_remove("_se"))

# 4. combine
dmnfr_age_pop_hosp_resp <- left_join(pop_long, se_long, by = "age")

save(file = "output/dmnfr_age_pop_hosp_resp.RData", dmnfr_age_pop_hosp_resp)

#-------------------------------------------------------------------------------
# Asthma symptom days
#-------------------------------------------------------------------------------

# 1. regional population + SE (correct variance rule)
pop_region <- asthma_symp_pop_benmap %>%
  dplyr::filter(str_detect(NAME, paste(county_of_interest, collapse = "|"))) %>%
  sf::st_drop_geometry() %>%
  dplyr::summarise(
    across(
      starts_with("pop_") & !ends_with("_se"),
      ~sum(.x, na.rm = TRUE)
    ),
    across(
      ends_with("_se"),
      ~sqrt(sum(.x^2, na.rm = TRUE))
    )
  )

# 2. pivot POP only
pop_long <- pop_region %>%
  dplyr::select(starts_with("pop_"), -ends_with("_se")) %>%
  pivot_longer(
    cols = everything(),
    names_to = "age",
    values_to = "pop"
  ) %>%
  mutate(age = str_remove(age, "pop_"))

# 3. pivot SE only
se_long <- pop_region %>%
  dplyr::select(ends_with("_se")) %>%
  pivot_longer(
    cols = everything(),
    names_to = "age",
    values_to = "se"
  ) %>%
  mutate(age = str_remove(age, "pop_") %>% str_remove("_se"))

# 4. combine
dmnfr_age_pop_asthma_symp <- left_join(pop_long, se_long, by = "age")

save(file = "output/dmnfr_age_pop_asthma_symp.RData", dmnfr_age_pop_asthma_symp)



#-------------------------------------------------------------------------------
# Minor restricted activity days
#-------------------------------------------------------------------------------

# 1. regional population + SE (correct variance rule)
pop_region <- MRAD_pop_benmap %>%
  dplyr::filter(str_detect(NAME, paste(county_of_interest, collapse = "|"))) %>%
  sf::st_drop_geometry() %>%
  dplyr::summarise(
    across(
      starts_with("pop_") & !ends_with("_se"),
      ~sum(.x, na.rm = TRUE)
    ),
    across(
      ends_with("_se"),
      ~sqrt(sum(.x^2, na.rm = TRUE))
    )
  )

# 2. pivot POP only
pop_long <- pop_region %>%
  dplyr::select(starts_with("pop_"), -ends_with("_se")) %>%
  pivot_longer(
    cols = everything(),
    names_to = "age",
    values_to = "pop"
  ) %>%
  mutate(age = str_remove(age, "pop_"))

# 3. pivot SE only
se_long <- pop_region %>%
  dplyr::select(ends_with("_se")) %>%
  pivot_longer(
    cols = everything(),
    names_to = "age",
    values_to = "se"
  ) %>%
  mutate(age = str_remove(age, "pop_") %>% str_remove("_se"))

# 4. combine
dmnfr_age_pop_mrad <- left_join(pop_long, se_long, by = "age")

save(file = "output/dmnfr_age_pop_mrad.RData", dmnfr_age_pop_mrad)


#-------------------------------------------------------------------------------
# School Loss Days - BenMap
## Assume school-age of 6-18
#-------------------------------------------------------------------------------

# 1. regional population + SE (correct variance rule)
pop_region <- SLD_pop_benmap %>%
  dplyr::filter(str_detect(NAME, paste(county_of_interest, collapse = "|"))) %>%
  sf::st_drop_geometry() %>%
  dplyr::summarise(
    across(
      starts_with("pop_") & !ends_with("_se"),
      ~sum(.x, na.rm = TRUE)
    ),
    across(
      ends_with("_se"),
      ~sqrt(sum(.x^2, na.rm = TRUE))
    )
  )

# 2. pivot POP only
pop_long <- pop_region %>%
  dplyr::select(starts_with("pop_"), -ends_with("_se")) %>%
  pivot_longer(
    cols = everything(),
    names_to = "age",
    values_to = "pop"
  ) %>%
  mutate(age = str_remove(age, "pop_"))

# 3. pivot SE only
se_long <- pop_region %>%
  dplyr::select(ends_with("_se")) %>%
  pivot_longer(
    cols = everything(),
    names_to = "age",
    values_to = "se"
  ) %>%
  mutate(age = str_remove(age, "pop_") %>% str_remove("_se"))

# 4. combine
dmnfr_age_pop_sld <- left_join(pop_long, se_long, by = "age")

save(file = "output/dmnfr_age_pop_sld.RData", dmnfr_age_pop_sld)

#-------------------------------------------------------------------------------
# Non-accidental mortality
#-------------------------------------------------------------------------------

# 1. regional population + SE (correct variance rule)
pop_region <- nonacc_mortality_pop_benmap %>%
  dplyr::filter(str_detect(NAME, paste(county_of_interest, collapse = "|"))) %>%
  sf::st_drop_geometry() %>%
  dplyr::summarise(
    across(
      starts_with("pop_") & !ends_with("_se"),
      ~sum(.x, na.rm = TRUE)
    ),
    across(
      ends_with("_se"),
      ~sqrt(sum(.x^2, na.rm = TRUE))
    )
  )

# 2. pivot POP only
pop_long <- pop_region %>%
  dplyr::select(starts_with("pop_"), -ends_with("_se")) %>%
  pivot_longer(
    cols = everything(),
    names_to = "age",
    values_to = "pop"
  ) %>%
  mutate(age = str_remove(age, "pop_"))

# 3. pivot SE only
se_long <- pop_region %>%
  dplyr::select(ends_with("_se")) %>%
  pivot_longer(
    cols = everything(),
    names_to = "age",
    values_to = "se"
  ) %>%
  mutate(age = str_remove(age, "pop_") %>% str_remove("_se"))

# 4. combine
dmnfr_age_pop_namort <- left_join(pop_long, se_long, by = "age")

save(file = "output/dmnfr_age_pop_namort.RData", dmnfr_age_pop_namort)



#-------------------------------------------------------------------------------
# Combine them all back as one set of dmnfr population at risk
#-------------------------------------------------------------------------------
dmnfr_age_pop_asthma_ed <- dmnfr_age_pop_asthma_ed %>%
  mutate(health_outcome = "Asthma_ED")

dmnfr_age_pop_hosp_resp <- dmnfr_age_pop_hosp_resp %>%
  mutate(health_outcome = "HOSP_resp")

dmnfr_age_pop_asthma_symp <- dmnfr_age_pop_asthma_symp %>%
  mutate(health_outcome = "Asthma_symptom_day")

dmnfr_age_pop_mrad <- dmnfr_age_pop_mrad %>%
  mutate(health_outcome = "MRAD")

dmnfr_age_pop_sld <- dmnfr_age_pop_sld %>%
  mutate(health_outcome = "School_loss_day")

dmnfr_age_pop_namort <- dmnfr_age_pop_namort %>%
  mutate(health_outcome = "Non_accidental_mortality")

dmnfr_age_pop_at_risk <- rbind(dmnfr_age_pop_asthma_ed,
                               dmnfr_age_pop_hosp_resp,
                               dmnfr_age_pop_asthma_symp,
                               dmnfr_age_pop_mrad,
                               dmnfr_age_pop_sld,
                               dmnfr_age_pop_namort
                               )
dmnfr_age_pop_at_risk <- dmnfr_age_pop_at_risk %>%
  dplyr::rename(pop_se = se) %>%
  mutate(pop_l95 = pop - 1.96 * pop_se,
         pop_u95 = pop + 1.96 * pop_se)

save(file = "output/dmnfr_age_pop_at_risk.RData", dmnfr_age_pop_at_risk)



######## Summary Statistics
library(dplyr)
library(knitr)

population_summary <- dmnfr_age_pop_at_risk %>%
  group_by(health_outcome) %>%
  summarize(
    n_age_groups = n(),
    total_pop = sum(pop, na.rm = TRUE),
    total_pop_se = sqrt(sum(pop_se^2, na.rm = TRUE)),
    total_pop_l95 = total_pop - 1.96 * total_pop_se,
    total_pop_u95 = total_pop + 1.96 * total_pop_se,
    min_age_pop = min(pop, na.rm = TRUE),
    max_age_pop = max(pop, na.rm = TRUE),
    .groups = "drop"
  )

knitr::kable(
  population_summary,
  digits = 0,
  caption = "Population at risk by health outcome"
)

write_csv(file = "output/population_summary.csv", x = population_summary)

population_age_summary <- dmnfr_age_pop_at_risk %>%
  arrange(health_outcome, age) %>%
  dplyr::select(
    health_outcome,
    age,
    pop,
    pop_se,
    pop_l95,
    pop_u95
  )

knitr::kable(
  population_age_summary,
  digits = 0,
  caption = "Age-specific population at risk used in the HIA"
)

write_csv(x = population_age_summary,
          file = "output/population_age_summary.csv")

