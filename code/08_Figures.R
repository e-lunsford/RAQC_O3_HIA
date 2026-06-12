#################################################################################
# title: Figures
# author: Beth Lunsford 
# date: 2026-06-08 
#
# This code is to set the area of interest for further spatial coding and
# map making.

# Last Run: 06/09/2026 and code was in working order using R 4.5.3
#################################################################################


# ------------------------------------------------------------------------------
# Set working directory to external hard drive CEBA folder.
# ------------------------------------------------------------------------------

setwd("M:/ALA_O3_HIA/RAQC_ALA_O3_HIA/code")

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
# Set map theme.
# ------------------------------------------------------------------------------
map_theme <- theme(
  #aspect.ratio = 1,
  text  = element_text(size = 12, color = 'black'),
  #panel.spacing.y = unit(0,"cm"),
  #panel.spacing.x = unit(0.25, "lines"),
  panel.grid.minor = element_line(color = "transparent"),
  panel.grid.major = element_line(color = "transparent"),
  panel.border = element_blank(),
  panel.background=element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  # legend.position = c(0.1,0.1),
  # plot.margin = grid::unit(c(0,0,0,0), "mm"),
  legend.key = element_blank(),
  #legend.background = element_rect(fill='transparent'),
  plot.margin = unit(c(0, 0, 0, 0), "cm"),
  panel.spacing = unit(c(0, 0, 0, 0), "cm"),
  legend.position = "inside",
  legend.position.inside = c(0.84, 0.4),
  legend.title = element_text(size = 12),
  legend.background = element_rect(fill = "white"),
  legend.spacing.y = unit(0.2, 'cm'),
  legend.margin = margin(0.1,0,0,0, unit="cm")
)

den_col <- "darkgrey"
oth_col <- "#4daf4a"
hwy_col <- "darkblue"

blue <- "#377eb8"
red <- "#e41a1c"
green <- "#4daf4a"
purple <- "#984ea3"
orange <- "#ff7f00"  
yellow <- "#ffff33"

# ------------------------------------------------------------------------------
## Figure: Monitor locations.
# ------------------------------------------------------------------------------

base_map1 <- ggplot() +
  geom_sf(data = dmnfr_tract_wgs84,
          inherit.aes = F,
          fill = "lightblue",
          colour = "darkblue",
          linewidth = .5) +
  geom_sf(data = ozone_aqs2_annual,
          inherit.aes = F,
          color = "black",
          size = 2.5,
          alpha = 0.8) +
  geom_sf(data = colorado_sf_wgs84,
          inherit.aes = F,
          fill = NA,
          color = "black") +
  geom_sf(data = colorado_counties_wgs84,
          inherit.aes = F,
          fill = NA,
          color = "black") +
  labs(caption = "Source: EPA AQS") +
  annotation_north_arrow(aes(location = "tl"),
                       height = unit(1, "cm"), 
                       width = unit(1, "cm"),
                       pad_x = unit(1, "cm"),
                       pad_y = unit(1, "cm")) +
  annotation_scale(aes(location = "br"),
                   pad_x = unit(1, "cm"),
                   pad_y = unit(1, "cm")) +
  map_theme +  
  labs(title = "",
       x = "",
       y = "",
       caption = "Source: US Census Bureau") + 
  theme(plot.caption = element_text(hjust = 0, vjust = 8),
        legend.position.inside = c(0.9, 0.4),
        plot.margin = margin(t = 10,  # Top margin
                             r = 20,  # Right margin
                             b = 10,  # Bottom margin
                             l = 10)) # Left margin

base_map1

ggplot2::ggsave(base_map1,
                filename = "figures/base_map1.jpeg",
                device = "jpeg",
                units = "in",
                height = 6,
                width = 8)









base_map2 <- ggplot() +
  geom_sf(data = dmnfr_buffer_wgs,
          fill = NA,
          color = "green",
          linewidth = 1) +
  geom_sf(data = ozone_dmnfr_annual,
          color = "red",
          size = 2) +
  geom_sf(data = DMNFR_counties_wgs84,
          fill = NA,
          color = "black") + 
  map_theme +
  labs(caption = "Source: EPA AQS") +
  annotation_north_arrow(aes(location = "tl"),
                         height = unit(1, "cm"), 
                         width = unit(1, "cm"),
                         pad_x = unit(1, "cm"),
                         pad_y = unit(1, "cm")) +
  annotation_scale(aes(location = "br"),
                   pad_x = unit(.5, "cm"),
                   pad_y = unit(1, "cm")) +
  map_theme +  
  labs(title = "",
       x = "",
       y = "",
       caption = "Source: US Census Bureau") + 
  theme(plot.caption = element_text(hjust = 0, vjust = 8),
        legend.position.inside = c(0.9, 0.4),
        plot.margin = margin(t = 10,  # Top margin
                             r = 20,  # Right margin
                             b = 10,  # Bottom margin
                             l = 10)) # Left margin

base_map2

ggplot2::ggsave(base_map2,
                filename = "figures/base_map2.jpeg",
                device = "jpeg",
                units = "in",
                height = 6,
                width = 8)




################################################################################
# Small multiples
multiples_plot <- ggplot(
  total_burden,
  aes(
    x = scenario,
    y = total_delta_y,
    shape = scenario
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = lower_total,
      ymax = upper_total
    ),
    width = 0.2
  ) +
  facet_wrap(
    ~ health_outcome,
    scales = "free_y"
  ) +
  scale_shape_manual(
    name = "Scenario",
    values = c(
      "delta_x_65" = 16,      # circle
      "delta_x_70" = 17,      # triangle
      "delta_x_plus10" = 15   # square
    ),
    labels = c(
      "delta_x_65" = "Exceedance > 65 ppb",
      "delta_x_70" = "Exceedance > 70 ppb",
      "delta_x_plus10" = "Ozone +10 ppb"
    )
  ) +
  labs(
    x = NULL,
    y = expression(Delta*"Y")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  )

multiples_plot3 <- ggplot(
  total_burden3,
  aes(
    x = scenario,
    y = total_delta_y,
    shape = scenario
  )
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = lower_total,
      ymax = upper_total
    ),
    width = 0.2
  ) +
  facet_wrap(
    ~ health_outcome,
    scales = "free_y",
    labeller = as_labeller(c(
      "Asthma_ED" = "Asthma ED",
      "HOSP_resp" = "Respiratory Hospitalizations",
      "MRAD" = "Minor Restricted Activity Days",
      "School_loss_day" = "School-Loss Days",
      "Non_accidental_mortality" = "Non-Accidental Mortality"
    ))
  ) +
  scale_shape_manual(
    name = "Scenario",
    values = c(
      "delta_x_65" = 16,      # circle
      "delta_x_70" = 17,      # triangle
      "delta_x_plus10" = 15   # square
    ),
    labels = c(
      "delta_x_65" = "Exceedance > 65 ppb",
      "delta_x_70" = "Exceedance > 70 ppb",
      "delta_x_plus10" = "Ozone +10 ppb"
    )
  ) +
  labs(
    x = NULL,
    y = "Attributable health burden (5-year ΔY)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  )

ggplot2::ggsave(multiples_plot3,
                filename = "figures/multiples_plot3.jpeg",
                device = "jpeg",
                units = "in",
                height = 6,
                width = 8)


# Heat map
# Which scenario generates the largest burden?
ggplot(total_burden2,
       aes(
         x = scenario,
         y = health_outcome,
         fill = total_delta_y
       )) +
  geom_tile()


# Forest plot


ggplot(total_burden,
       aes(
         y = health_outcome,
         x = total_delta_y,
         xmin = lower_total,
         xmax = upper_total,
         color = scenario
       )) +
  geom_pointrange(
    position = position_dodge(width = 0.5)
  )

















################################################################################
# Asthma age-specific ED
################################################################################

# ------------------------------------------------------------------------------
# Figure of uncertainty
# ------------------------------------------------------------------------------

asthma_ed_hia_plot <- ggplot(
  data = asthma_ed_hia_results, 
  aes(x = age, 
      y = mean_delta_y, 
      color = scenario)) +
  
  geom_point(position = position_dodge(width = 0.4), size = 2) +
  
  geom_errorbar(
    aes(ymin = lower_95, ymax = upper_95),
    position = position_dodge(width = 0.4),
    width = 0.2
  ) +
  
  labs(
    title = "Asthma ED visits attributable to ozone exposure",
    x = "Age Group",
    y = "Attributable Cases (ΔY)",
    color = "Scenario"
  ) +
  
  theme_minimal()

ggplot2::ggsave(asthma_ed_hia_plot,
                filename = "figures/asthma_ed_hia_plot.jpeg",
                device = "jpeg",
                units = "in",
                height = 6,
                width = 8)


################################################################################
# Hospital - respiratory
################################################################################



# ------------------------------------------------------------------------------
# Figure of uncertainty
# ------------------------------------------------------------------------------

hosp_resp_hia_plot <- 
  ggplot(hosp_resp_hia_results, aes(x = age, y = mean_delta_y, color = scenario)) +
  
  geom_point(position = position_dodge(width = 0.4), size = 2) +
  
  geom_errorbar(
    aes(ymin = lower_95, ymax = upper_95),
    position = position_dodge(width = 0.4),
    width = 0.2
  ) +
  
  labs(
    title = "Respiratory hospitalizations attributable to ozone exposure",
    x = "Age Group",
    y = "Attributable Cases (ΔY)",
    color = "Scenario"
  ) +
  
  theme_minimal()

ggplot2::ggsave(hosp_resp_hia_plot,
                filename = "figures/hosp_resp_hia_plot.jpeg",
                device = "jpeg",
                units = "in",
                height = 6,
                width = 8)

################################################################################
# non-accidental mortality




# ------------------------------------------------------------------------------
# Figure of uncertainty
# ------------------------------------------------------------------------------

non_accidental_mortality_hia_plot <- 
  ggplot(Non_accidental_mortality_hia_results, 
         aes(x = age, y = mean_delta_y, color = scenario)) +
  
  geom_point(position = position_dodge(width = 0.4), size = 2) +
  
  geom_errorbar(
    aes(ymin = lower_95, ymax = upper_95),
    position = position_dodge(width = 0.4),
    width = 0.2
  ) +
  
  labs(
    title = "Non-accidental mortality attributable to ozone exposure",
    x = "Age Group",
    y = "Attributable Cases (ΔY)",
    color = "Scenario"
  ) +
  
  theme_minimal()

ggplot2::ggsave(non_accidental_mortality_hia_plot,
                filename = "figures/non_accidental_mortality_hia_plot.jpeg",
                device = "jpeg",
                units = "in",
                height = 6,
                width = 8)


################################################################################
# mrad


mrad_hia_results <- read_csv(file = "output/mrad_hia_results.csv")





# ------------------------------------------------------------------------------
# Figure of uncertainty
# ------------------------------------------------------------------------------

mrad_hia_results_plot_data <- mrad_hia_results %>%
  mutate(
    age = recode(
      age,
      "18_64" = "18-64"
    ),
    scenario = recode(
      scenario,
      "delta_x_65" = "Exceedance above 0.65 ppb",
      "delta_x_70" = "Exceedance above 0.70 ppb",
      "delta_x_plus10" = "Ozone increase by 0.10 ppb"
    )
  )


mrad_hia_results_plot <- 
  ggplot(mrad_hia_results_plot_data, aes(x = age, y = mean_delta_y, color = scenario)) +
  
  geom_point(position = position_dodge(width = 0.4), size = 2) +
  
  geom_errorbar(
    aes(ymin = lower_95, ymax = upper_95),
    position = position_dodge(width = 0.4),
    width = 0.2
  ) +
  
  labs(
    title = "Minor restricted activity days attributable to ozone exposure",
    x = "Age Group",
    y = "Attributable Cases (ΔY)",
    color = "Scenario"
  ) +
  
  theme_minimal()


ggplot2::ggsave(mrad_hia_results_plot,
                filename = "figures/mrad_hia_results_plot.jpeg",
                device = "jpeg",
                units = "in",
                height = 6,
                width = 8)

################################################################################
# sld
################################################################################



# ------------------------------------------------------------------------------
# Figure of uncertainty
# ------------------------------------------------------------------------------

ggplot(sld_hia_results, aes(x = age, y = mean_delta_y, color = scenario)) +
  
  geom_point(position = position_dodge(width = 0.4), size = 2) +
  
  geom_errorbar(
    aes(ymin = lower_95, ymax = upper_95),
    position = position_dodge(width = 0.4),
    width = 0.2
  ) +
  
  labs(
    title = "sld Attributable to Ozone Exposure",
    x = "Age Group",
    y = "Attributable Cases (ΔY)",
    color = "Scenario"
  ) +
  
  theme_minimal()


























