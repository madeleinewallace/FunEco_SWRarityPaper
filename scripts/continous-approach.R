# ==============================================================================
# --- info ---------------------------------------------------------------------
# okay readjusting after failed climwin analysis
# new plan, use SPEI3 lol

# do rare perennial species gain cover relative to common species
# over drier survey intervals? (observational version of Watkins et al. 2026,
# Ecology 107:e70478) at Jornada (NM), Flagstaff (AZ) and Santa Rita (AZ)

# october 28 2026
# maddie wallace
# ==============================================================================



# ==============================================================================
# --- packages -----------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(readxl)
# ==============================================================================



# ==============================================================================
# --- stuff i need -------------------------------------------------------------
#data_dir    <- "~/Desktop/01-research/FunEco_SWRarityPaper/data"
#climate_dir <- file.path(data_dir, "climate")
#out_dir     <- file.path(data_dir, "rarity_drought_output")   # new folder
#dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sites <- c("jornada", "flagstaff", "santarita")

# site-year counts only if this many quadrats/transects were read
min_units <- 10

# jornada provides survey dates, it's mostly between oct-dec ishhh
# SRER source studies say late Aug to early Sept
# flag data source says Aug-Sept
fixed_survey <- list(flagstaff = c(month = 9, day = 1),
                     santarita = c(month = 9, day = 1))

# SPEI to calculate
spei_scales  <- c(1, 3, 6, 12)
primary_spei <- "spei03" # can be changed idk what to do after climwin lol

site_points <- tribble(
  ~site,       ~point,     ~lat,    ~lon,
  "jornada",   "centroid", 32.62, -106.74,
  "flagstaff", "centroid", 35.27, -111.74,
  "santarita", "centroid", 31.83, -110.85)

# drop survey intervals longer than this (to remove jornada's weird 1980-1995 gap in a min)
max_gap_years <- 5
# ==============================================================================



# ==============================================================================
# --- clean cover --------------------------------------------------------------
site_data <- function(df, site_name) {
  df |> 
    filter(!is.na(total_cover)) |> 
    group_by(year) |> 
    summarize(n_units = n_distinct(unit),
              cover = mean(total_cover),
              cover_sd = sd(total_cover),
              survey_date = median(survey_date),
              .groups = "drop") |> 
    mutate(site = site_name, .before = 1) #forces to first column in new df <3
}

# JORNADA
jorn_cover <- read_csv("data/jornada/Jornada_permanent_quads_cover.csv",
                       col_types = cols(.default = col_character()))

jrn_quad <- jorn_cover |> 
  filter(chart_status == "processed") |> 
  mutate(year = as.integer(year),
         month = as.integer(month),
         day = as.integer(day),
         area = as.numeric(area)) |> 
  filter(!is.na(month)) |> 
  mutate(day = replace_na(day, 15)) |> 
  group_by(quadrat, year, month, day) |> 
  summarize(total_cover = 100 * sum(area[cover_status == "PG/Shr present"], na.rm = T),
            n_polygons = sum(cover_status == "PG/Shr present"), .groups = "drop") |> 
  mutate(survey_date = make_date(year, month, day))

jrn_quad <- jrn_quad |> 
  group_by(quadrat) |> 
  mutate(month_middle = as.integer(names(which.max(table(month)))),
         month_dist = pmin(abs(month - month_middle), 12 - abs(month - month_middle))) |> 
  group_by(quadrat, year) |> 
  arrange(month_dist, survey_date, .by_group = TRUE) |> 
  mutate(n_charts_in_year = n()) |> 
  slice(1) |> 
  ungroup() |> 
  transmute(unit = quadrat, year, survey_date, total_cover, n_charts_in_year)

jrn_site <- site_data (jrn_quad, "jornada")

jrn_spp <- jorn_cover |> 
  filter(chart_status == "processed", cover_status == "PG/Shr present",
         scientific_name != ".") |> 
  mutate(year = as.integer(year),
         month = as.integer(month),
         day = replace_na(as.integer(day), 15),
         area = as.numeric(area),
         survey_date = make_date(year, month, day)) |> 
  semi_join(jrn_quad, by = c("quadrat" = "unit", "year", "survey_date")) |> 
  group_by(year, species = scientific_name) |> 
  summarize(total_cover = 100 * sum(area, na.rm = TRUE), .groups = "drop") |> 
  mutate(site = "jornada")

# FLAGSTAFF
flag_cover <- read_csv("data/flag/Basal_Cover_per_Species_by_Quadrat_and_Year.csv")

flag_quad <- flag_cover |> 
  mutate(area_cm2 = if_else(Species == "No Cover Species Observed", 0, Area_Sq_Cm)) |> 
  group_by(Site, Quadrat, Year) |> 
  summarize(total_cover = sum(area_cm2, na.rm = TRUE) / 100, .groups = "drop") |> 
  transmute(unit = paste(Site, Quadrat, sep = "_"), subsite = Site, year = Year,
            survey_date = make_date(year, fixed_survey$flagstaff["month"], fixed_survey$flagstaff["day"]),
            total_cover)

flag_site <- site_data(flag_quad, "flagstaff")

flag_spp <- flag_cover |> 
  filter(!Species %in% c("No Cover Species Observed", "Unknown graminoid")) |> 
  group_by(year = Year, species = Species) |> 
  summarize(total_cover = sum(Area_Sq_Cm, na.rm = TRUE) / 100, .groups = "drop") |> 
  mutate(site = "flagstaff")

# SRER
# -9999 = not measured
# 0 = not observed
srer_cover <- read_excel("data/srer/ltcover_2024Dec20.xlsx") |> 
  mutate(across(starts_with("YR"), as.numeric)) |> 
  pivot_longer(starts_with("YR"), names_to = "year", names_prefix = "YR",
               names_transform = as.integer, values_to = "cover") |> 
  mutate(cover = na_if(cover, -9999))

srer_tran <- srer_cover |>
  group_by(PASTURE, TRANSECT, year, CATEGORY) |>
  summarise(group_cover = if (all(is.na(cover))) NA_real_ else sum(cover, na.rm = TRUE),
            .groups = "drop") |>
  group_by(PASTURE, TRANSECT, year) |>
  summarise(total_cover = if (any(is.na(group_cover))) NA_real_ else sum(group_cover),
            .groups = "drop") |>
  transmute(unit = TRANSECT, year,
            survey_date = make_date(year, fixed_survey$santarita["month"],
                                    fixed_survey$santarita["day"]),
            total_cover)

srer_site <- site_data(srer_tran, "santarita")

srer_spp <- srer_cover |>
  filter(!SPECIES %in% c("OTHERPG", "OTHSHRUB")) |>
  summarize(cover = mean(cover, na.rm = TRUE), .groups = "drop") |>
  filter(!is.nan(cover)) |>
  mutate(site = "santarita")
# ==============================================================================



# ==============================================================================
# --- species x site x survey cover --------------------------------------------

