library(tidyverse)
library(SPEI)

srer_clim <- read_csv("data/climate_NA/srer/srer_units_latlon_1949-2024MP(in).csv")
# how this corresponds to cover data:
# id1 = PASTTRAN
# id2 = PASTURE
# Tmax01 - Tmax12 : monthly maximum temps for months 1-12 (jan-dec)
# same with Tmin01 - Tmin12; Tave01-Tave12; PPT01-PPT12; rsds01-rsds12
# last 2 characters of those column names = month

# first, pivot long?

srer_long <- srer_clim %>% 
  pivot_longer(cols = Tmax01:rsds12, names_to = "var", values_to = "value") %>% 
  mutate(month = as.numeric(str_sub(var, -2, -1)), variable = str_sub(var, 1, 3)) %>% 
  dplyr::select(-var) %>% 
  mutate(var = case_when(variable == "Tma" ~ "tmax", variable == "Tmi" ~ "tmin",
                         variable == "Tav" ~ "tave", variable == "PPT" ~ "ppt",
                         variable == "rsd" ~ "rad")) %>% 
  dplyr::select(-variable)

srer_wide <- srer_long %>% 
  pivot_wider(names_from = var, values_from = value) %>% 
  mutate(date = make_date(year = Year, month = month, day = 1),
         unit = id1, elevation = Elevation, lati = Latitude) %>% 
  dplyr::select(date, unit, elevation, tmax, tmin, tave, ppt, rad, lati)

# srer_lat <- srer

ggplot(srer_wide, aes(x = date, y = tave))+
  geom_line(aes(group = unit))+
  geom_smooth(method = "lm")


srer_pet <- srer_wide %>% 
  group_by(unit) %>% 
  mutate(pet_1 = hargreaves(tmin, tmax, lat = unique(lati)), # calculate with latitude-based radiation; one value per site
         pet_2 = hargreaves(tmin, tmax, Ra = rad)) # calculate with climateNA radiation output (need clarity from ceci on what this is from?)

ggplot(srer_pet, aes(x = date, y = pet_2))+
  geom_line(aes(group = unit, color = unit), alpha = 0.4)

srer_spei_raw <- srer_pet %>% 
  mutate(water_balance = ppt - pet_2)

srer_spei <- srer_spei_raw %>% 
  group_by(unit) %>% 
  arrange(date) %>% 
  mutate(spei12 = spei(water_balance, 12)$fitted)


ggplot(filter(srer_spei, year(date)>=2020), aes(x = date, y = spei12))+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_line(aes(group = unit, color = unit), alpha = 0.4)+
  theme_minimal(base_size = 16)+
  labs(x = "Date", y = "SPEI-12")


