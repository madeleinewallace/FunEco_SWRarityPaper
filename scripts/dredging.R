# maddie wallace
# oct 29 2026

# okay try with the dredge function, try with MAP, MAT, PDSI, VPD
# which climate variable best predicts total plant cover at each site?
# dredge version: SPEI (all scales), MAT, MAP, VPD, PDSI

#install.packages("MuMIn") 
#remotes::install_github("Sibada/scPDSI")
library(tidyverse)
library(lubridate)
library(prism)
library(terra)
library(MuMIn)
library(scPDSI)
library(nlme)

window_months <- 12   # climate window
awc_by_site <- c(jornada = 110, santarita = 114, flagstaff = 118)  # i can gert this from the NRCS web soil survey
r_max <- 0.7  
max_terms <- 2

# get extra PRISM data (VPD)
prism_set_dl_dir("data/climate")
#get_prism_monthlys(type = "vpdmin", years = 1915:2024, mon = 1:12, keepZip = FALSE)
#get_prism_monthlys(type = "vpdmax", years = 1915:2024, mon = 1:12, keepZip = FALSE)

site_points <- tribble(
  ~site,       ~point,       ~lat,    ~lon,
  "jornada",   "centroid",   32.62, -106.74,
  "flagstaff", "centroid",   35.27, -111.74,
  "santarita", "centroid",   31.83, -110.85)

# function to extract monthly PRISM data for points
# i am proud of this i am pracitcing functions lol
extract_prism <- function(type, points) {
  # 1. list of PRISM archive dataset folder names
  pd <- prism_archive_subset(type, "monthly", resolution = "4km")
  
  # 2. get dates from folder names
  dates <- as.Date(pd_get_date(pd))
  
  # 3. get file paths ordered chronologically
  files <- pd_to_file(pd[order(dates)])
  
  # 4. load rasters and extract values at points above
  r <-rast(files)
  pts <- vect(points, geom = c("lon", "lat"), crs = "EPSG:4326")
  
  vals <- terra::extract(r, pts, ID = FALSE)
  colnames(vals) <- as.character(sort(dates))
  vals$site <- points$site
  vals$point_id <- points$point_id
  
  # 5. wide columns into long format
  pivot_longer(vals, cols = -c(site, point_id), names_to = "date", values_to = type) |> 
    mutate(date = as.Date(date))
}

site_points <- site_points %>% mutate(point_id = paste0(site, "_", point))

tmean_df <- extract_prism("tmean", site_points)
vpdmin_df <- extract_prism("vpdmin", site_points)
vpdmax_df <- extract_prism("vpdmax", site_points)
ppt_df  <- extract_prism("ppt", site_points)
tmin_df <- extract_prism("tmin", site_points)
tmax_df <- extract_prism("tmax", site_points)

spei_site <- read_csv("data/spei_monthly_by_site.csv")

climate <- spei_site |> 
  left_join(select(tmean_df,  site, date, tmean),  by = c("site", "date")) |>
  left_join(select(vpdmin_df, site, date, vpdmin), by = c("site", "date")) |>
  left_join(select(vpdmax_df, site, date, vpdmax), by = c("site", "date")) |>
  mutate(vpd = (vpdmin + vpdmax) / 2) |>
  arrange(site, date)

# PDSI per site with ppt and Hargreaves PET
climate <- climate |>
  group_by(site) |>
  mutate(pdsi = as.numeric(pdsi(P = ppt, PE = pet, AWC = awc_mm, sc = TRUE)$X)) |>
  ungroup()

climate <- climate |>
  group_by(site) |>
  mutate(pdsi = as.numeric(pdsi(P = ppt, PE = pet,
                                        AWC = awc_by_site[site[1]], sc = TRUE)$X)) |>
  ungroup()

# summarize climate in the window before each survey
cover_site <- read_csv("data/cover_site_year.csv") |> 
  filter(keep)

cover_climate <- cover_site |>
  mutate(survey_month = floor_date(survey_date, "month"),
         win_start = survey_month %m-% months(window_months),
         win_end   = survey_month %m-% months(1)) |>
  left_join(climate, by = "site", relationship = "many-to-many") |>
  filter(date >= win_start, date <= win_end) |>
  arrange(site, year, date) |>
  group_by(site, year) |>
  summarize(cover = first(cover),
            MAP = sum(ppt), # total precip in window
            MAT = mean((tmin + tmax) / 2), # mean temp in window, i think the prism download was incomplete lol
            VPD = mean(vpd), # mean VPD in window
            PDSI = last(pdsi), # PDSI month before survey
            across(starts_with("spei"), last), # each SPEI scale, month before survey
            n_months = n(),
            .groups = "drop") |>
  filter(n_months == window_months) |>
  drop_na()

head(cover_climate)

write_csv(cover_climate, "data/cover_climate.csv")

# okay now dredge
options(na.action = "na.fail")  # dredge refuses to run without this ??????
spei_var <- grep("^spei", names(cover_climate), value = TRUE) # iver always wanted to grep
clim_var <- c("MAT", "MAP", "VPD", "PDSI", spei_var)

# JORNADA
?dredge
d_jornada <- cover_climate |> 
  filter(site == "jornada") |> 
  arrange(year) |> 
  mutate(across(all_of(clim_var), ~ as.numeric(scale(.x))))

d_jornada <- d_jornada |> select(cover, all_of(clim_var))

model_jornada <- lm(cover ~ ., data = d_jornada)

dd_single_jornada <- dredge(model_jornada, m.lim = c(0, 1), rank = "AICc")
dd_pairs_jornada  <- dredge(model_jornada, m.lim = c(0, max_terms), rank = "AICc")

jorn1.best <- get.models(dd_single_jornada, delta==0)[[1]]
summary(jorn1.best)

jorn2.best <- get.models(dd_pairs_jornada, delta==0)[[1]]
summary(jorn2.best)

# FLAG
d_flag <- cover_climate |> 
  filter(site == "flagstaff") |> 
  arrange(year) |> 
  mutate(across(all_of(clim_var), ~ as.numeric(scale(.x))))

d_flag <- d_flag |> select(cover, all_of(clim_var))

model_flag <- lm(cover ~ ., data = d_flag)

dd_single_flag <- dredge(model_flag, m.lim = c(0, 1), rank = "AICc")
dd_pairs_flag  <- dredge(model_flag, m.lim = c(0, max_terms), rank = "AICc")

flag1.best <- get.models(dd_single_flag, delta==0)[[1]]
summary(flag1.best)

flag2.best <- get.models(dd_pairs_flag, delta==0)[[1]]
summary(flag2.best)

# SRER
d_srer <- cover_climate |> 
  filter(site == "santarita") |> 
  arrange(year) |> 
  mutate(across(all_of(clim_var), ~ as.numeric(scale(.x))))

d_srer <- d_srer |> select(cover, all_of(clim_var))

model_srer <- lm(cover ~ ., data = d_srer)

dd_single_srer <- dredge(model_srer, m.lim = c(0, 1), rank = "AICc")
dd_pairs_srer  <- dredge(model_srer, m.lim = c(0, max_terms), rank = "AICc")

srer1.best <- get.models(dd_single_srer, delta==0)[[1]]
summary(srer1.best)

srer2.best <- get.models(dd_pairs_srer, delta==0)[[1]]
summary(srer2.best)

dd_single_srer
dd_pairs_srer

options(na.action = "na.omit")


# HMMMMMMMM