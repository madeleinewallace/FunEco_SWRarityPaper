# trying to figure out climate window stuff across the three sites
# which SPEI timescale best predicts total plant cover?
# follows Rubio-Cuadrado et al. (2022, Dendrochronologia) methods sorta
# oct 18 2026
# maddie wallace

# packages
# ==============================================================================
#install.packages("climwin")
#install.packages("SPEI")
library(climwin)
library(SPEI)
library(readxl)
library(dplyr)
library(lme4)
library(prism)
library(terra)
library(tidyverse)
library(lubridate)
library(lmerTest)

# get climate data
# ==============================================================================
# i'm getting tmin and tmax too bcus i think hargreaves method is better for arid sites?
# need to dig up that reference
# sry it's so much data
prism_set_dl_dir("data/climate")
#get_prism_monthlys(type = "ppt", years = 1915:2024, mon = 1:12, keepZip = FALSE)
#get_prism_monthlys(type = "tmean", years = 1915:2024, mon = 1:12, keepZip = FALSE)
#get_prism_monthlys(type = "tmin", years = 1915:2024, mon = 1:12, keepZip = FALSE)
#get_prism_monthlys(type = "tmax", years = 1915:2024, mon = 1:12, keepZip = FALSE)

ppt_files <- prism_archive_subset("ppt", "monthly", resolution = "4km")
tmean_files <- prism_archive_subset("tmean", "monthly", resolution = "4km")
tmin_files <- prism_archive_subset("tmin", "monthly", resolution = "4km")
tmax_files <- prism_archive_subset("tmax", "monthly", resolution = "4km")

ppt_rast <- rast(pd_to_file(ppt_files))
tmean_rast <- rast(pd_to_file(tmean_files))
tmin_rast <- rast(pd_to_file(tmin_files))
tmax_rast <- rast(pd_to_file(tmax_files))

# stuff i need i think, come back to this section to change stuff
# ==============================================================================
#october 1st date for sites that lack daily survey dates in raw data
#jornada provides survey dates, it's mostly between oct-dec ishhh
#SRER source studies say late Aug to early Sept
#flag data source says Aug-Sept
fixed_survey <- list(flagstaff = c(month = 9, day = 1),
                     santarita = c(month = 9, day = 1))

# get cover data and organize it
# ==============================================================================
flag_cover <- read_csv("data/flag/Basal_Cover_per_Species_by_Quadrat_and_Year.csv")
jorn_cover <- read_csv("data/jornada/Jornada_permanent_quads_cover.csv")
srer_cover <- read_excel("data/srer/ltcover_2024Dec20.xlsx") |> 
  mutate(across(starts_with("YR"), as.numeric)) |> 
  pivot_longer(starts_with("YR"), names_to = "year", names_prefix = "YR",
               names_transform = as.integer, values_to = "cover") |> 
  mutate(cover = na_if(cover, -9999))

#made a function to take quad or tran data and calculate summary data 
# cus it was a lot of code
# unqiue count of units (quad or tran), mean cover, SD of cover, median survey date,
# and add site as first column
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


# jornada
# only keeps proccessed data , converst dates and area to numeric
# sums polygon and multiply by 100 to get total % cov for 1m2 quads
# adds date
# lots of missing days in this dataset, so force it to the 15th
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

# there are some years where quads were surveyed twice
# i want the ones closest to the fall aka the typical survey time
# this took forever i used the evil ocomputer machine to unstick me
# why survey twice a year??
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

?slice() #nice

# okay now use function to get cover data
jrn_site <- site_data (jrn_quad, "jornada")

#flagstaff
# sets spp cover area to 0 for empty quadrats
# sums plant area per quad per yr
# converts cm2 to percent cover
# assigns survey date as mentioned above
flag_quad <- flag_cover |> 
  mutate(area_cm2 = if_else(Species == "No Cover Species Observed", 0, Area_Sq_Cm)) |> 
  group_by(Site, Quadrat, Year) |> 
  summarize(total_cover = sum(area_cm2, na.rm = TRUE) / 100, .groups = "drop") |> 
  transmute(unit = as.character(Quadrat), subsite = Site, year = Year,
            survey_date = make_date(year, fixed_survey$flagstaff["month"], fixed_survey$flagstaff["day"]),
            total_cover)

# put in function
flag_site <- site_data(flag_quad, "flagstaff")

# SRER
# i hate our formatting at SRER :)
srer_group <- srer_cover |>
  group_by(PASTURE, TRANSECT, year, CATEGORY) |>
  summarise(n_species = n(),
            n_missing = sum(is.na(cover)),
            group_cover = if (all(is.na(cover))) NA_real_ else sum(cover, na.rm = TRUE),
            .groups = "drop")

# sums total plant category cover to get total transect cover per year
# importantly not in percent cover, but i think it doesnt matter because i dont (yet) compare cover btwn sites
srer_transect <- srer_group |>
  group_by(PASTURE, TRANSECT, year) |>
  summarise(n_groups_unmeasured = sum(is.na(group_cover)),
            n_species_missing   = sum(n_missing),
            total_cover = if (any(is.na(group_cover))) NA_real_ else sum(group_cover),
            .groups = "drop")

# assigns survey date
srer_transect <- srer_transect |>
  transmute(unit = TRANSECT, pasture = PASTURE, year,
            survey_date = make_date(year, fixed_survey$santarita["month"], 
                                    fixed_survey$santarita["day"]),
            total_cover, n_species_missing)

# put into fucntion
srer_site <- site_data(srer_transect, "santarita")

#combining all this data just in case i need it later
#importantly, also tosses out a year if there are less than 10 units per year
# cough cough jornada
# necessary for climwin!
cover_site <- bind_rows(jrn_site, flag_site, srer_site) |> 
  mutate(keep = n_units >= 10)

# puts data into a csv and rds to pass it onto a new script if needed
write_csv(cover_site, file.path("data/cover_site_year.csv"))
write_rds(list(jornada = jrn_quad, flagstaff = flag_quad, santarita = srer_transect,
               santarita_groups = srer_group),
          file.path("data/cover_unit_level.rds"))


# PRISM and SPEI stuff
# oh lord
# much more wrangling than i expected
# ==============================================================================
# prepare point IDs
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

# prepare point IDs
site_points <- site_points %>% mutate(point_id = paste0(site, "_", point))

# get climate variables
ppt_df  <- extract_prism("ppt", site_points)
tmin_df <- extract_prism("tmin", site_points)
tmax_df <- extract_prism("tmax", site_points)

# combine and calculate water balance (precip - PET using hargreves)
climate_df <- ppt_df |> 
  # combine variables from above 
  left_join(tmin_df, by = c("site", "point_id", "date")) |> 
  left_join(tmax_df, by = c("site", "point_id", "date")) |> 
  left_join(site_points |>  group_by(site) |>  summarise(lat = mean(lat)), by = "site") |> 
  group_by(site, date) |> 
  summarise(
    ppt  = mean(ppt, na.rm = TRUE),
    tmin = mean(tmin, na.rm = TRUE),
    tmax = mean(tmax, na.rm = TRUE),
    lat  = mean(lat, na.rm = TRUE),
    .groups = "drop") |> 
  arrange(site, date) |> 
  # group by site to compute Hargreaves PET per site (package needs just one lat value i guess not a vector)
  group_by(site) |> 
  mutate(pet = as.numeric(hargreaves(Tmin = tmin, Tmax = tmax, lat = lat[1]))) |> 
  mutate(wb  = ppt - pet) |> 
  ungroup()

?hargreaves

# now, calc SPEI across these scales
spei_scales <- c(1, 3, 6, 9, 12, 24, 36, 48)

?spei

# another function :D it finally worked
calc_site_spei <- function(df, scales) {
  start_ym <- c(year(min(df$date)), month(min(df$date))) 
  wb_ts <- ts(df$wb, start = start_ym, frequency = 12) #time series required to go into spei()
  
  for (s in scales) {
    col_name <- sprintf("spei%02d", s) #returns a chatacter containing this format ex. spei01
    res <- spei(wb_ts, scale = s, verbose = FALSE)$fitted # spei function and extracts numbers
    res_vec  <- as.numeric(res) # converts ts output from spei() into numeric vector
    res_vec[!is.finite(res_vec)] <- NA # wasnt working without this idk
    df[[col_name]] <- res_vec # puts spei values back into dataframe
  }
  return(df)
}

# splits dataset by site
# runs calc_site_spei() on each site chunk
# recombines results and ungroups
spei_site <- climate_df |> 
  group_by(site) |> 
  group_modify(~ calc_site_spei(.x, spei_scales)) |> 
  ungroup()
#view(spei_site)

write_csv(spei_site, "data/spei_monthly_by_site.csv")


# find SPEI timescale per site
# ==============================================================================
window_range <- c(48, 1)  # windows can open up to 48 months before the survey ?
n_rand <- 1000 # randomizations (1000 as in Rubio-Cuadrado et al)

# JORNADA
#combine climate and cover
clim_jornada <- spei_site |> 
  filter(site == "jornada") |> 
  drop_na(starts_with("spei"))

# climwin needs daily inputs??? trying to expand monthly data in my df
# idk or maybe we could download daily data but that seems like an F ton of data from 1915
clim_daily <- data.frame(date = seq(min(clim_jornada$date), max(clim_jornada$date), by = "day")) |> 
  mutate(month = floor_date(date, "month"),
         cdate = format(date, "%d/%m/%Y")) |> 
  left_join(clim_jornada, by = c("month" = "date"))

# filter plant cover data from big df
min_clim_date <- min(clim_jornada$date) # this is bcus there is missing SPEI at the beg of the record

cover_jornada <- cover_site |> 
  filter(site == "jornada", keep == TRUE,
         survey_date >= (min_clim_date %m+% months(window_range[1])), 
         survey_date <= max(clim_daily$date)) |> 
  mutate(bdate = format(survey_date, "%d/%m/%Y"))

# ok finally use climwin, at line 300 lolz
spei_list <- as.list(select(clim_daily, starts_with("spei")))

?slidingwin # so many things you can chnage or pick
sw_jornada <- slidingwin(xvar = spei_list, # climate variables we are interested in
                         cdate = clim_daily$cdate, #climate date
                         bdate = cover_jornada$bdate, #biological date
                         baseline = lm(cover ~ 1, data = cover_jornada), # can do anything here whcih is kinda cool but i think we just keep it simple
                         cinterval = "month", # no clue why i needed to put it into daily format if i can pick this but it wouldnt work otherwise
                         range = window_range, # set above..... idk how to make this decision
                         type = "relative",
                         stat = "mean",
                         func = "lin",
                         cmissing = FALSE) 

# find best SPEI scale and window from function above
# row index of model above with lowest AICc, whichmin is from climwin??
best_row <- which.min(sw_jornada$combos$DeltaAICc)
best_scale <- as.character(sw_jornada$combos$climate[best_row])
best_scale
# YIPPPPEEEEEEE in agreement with a paper i read from jornada
best_model <- sw_jornada[[best_row]]

# run climwin's randomiztion test with best scale from above
?randwin
rand_jornada <- randwin(repeats = n_rand,
                        window = "sliding",
                        xvar = list(spei = clim_daily[[best_scale]]),
                        cdate = clim_daily$cdate,
                        bdate = cover_jornada$bdate,
                        baseline = lm(cover ~ 1, data = cover_jornada),
                        cinterval = "month",
                        range = window_range,
                        type = "relative",
                        stat = "mean",
                        func = "lin",
                        cmissing = FALSE)

# calc p value and plot diagnostics from climwin package
# computes p-value for AICc relative to randomized runs
p_aicc <- as.numeric(pvalue(best_model$Dataset, rand_jornada[[1]], metric = "AIC"))
# protects against missing numeric conversions by defaulting to 0 idk wouldn't work without it
if (is.na(p_aicc)) p_aicc <- 0
# calculates randomization probability metric Pc
p_c <- pvalue(best_model$Dataset, rand_jornada[[1]], metric = "C", sample.size = nrow(cover_jornada))

# these are the output diagnostics from the climwin package and from the paper
pdf("data/climwin_output/climwin_plots_jornada.pdf", width = 12, height = 8)
plotall(dataset= best_model$Dataset,
        datasetrand = rand_jornada[[1]],
        bestmodel = best_model$BestModel,
        bestmodeldata = best_model$BestModelData)
dev.off()

# final summary for jornada climwin for lab meeting
results_jornada <- data.frame(site = "jornada",
                              spei_scale = best_scale,
                              window_open = best_model$Dataset$WindowOpen[1],
                              window_close = best_model$Dataset$WindowClose[1],
                              delta_aicc = best_model$Dataset$deltaAICc[1],
                              slope = best_model$Dataset$ModelBeta[1],
                              adj_r2 = summary(best_model$BestModel)$adj.r.squared,
                              p_delta_aicc = p_aicc,
                              p_c = p_c,
                              significant = p_aicc < 0.05,
                              n_years = nrow(cover_jornada),
                              first_year = min(cover_jornada$year),
                              last_year = max(cover_jornada$year))
results_jornada

# okay now i copy and pste cus im sick of functions

# FLAGSTAFF
#combine climate and cover
clim_flag <- spei_site |> 
  filter(site == "flagstaff") |> 
  drop_na(starts_with("spei"))

# climwin needs daily inputs? trying to expand monthly data in my df
clim_daily_flag <- data.frame(date = seq(min(clim_flag$date), max(clim_flag$date), by = "day")) |> 
  mutate(month = floor_date(date, "month"),
         cdate = format(date, "%d/%m/%Y")) |> 
  left_join(clim_flag, by = c("month" = "date"))

# filter plant cover data from big df
min_clim_date_flag <- min(clim_flag$date) # this is bcus there is missing SPEI at the beg of the record

cover_flag <- cover_site |> 
  filter(site == "flagstaff", keep == TRUE,
         survey_date >= (min_clim_date_flag %m+% months(window_range[1])), 
         survey_date <= max(clim_daily_flag$date)) |> 
  mutate(bdate = format(survey_date, "%d/%m/%Y"))

# ok finally use climwin, at line 300 lolz
spei_list_flag <- as.list(select(clim_daily_flag, starts_with("spei")))

?slidingwin 
sw_flag <- slidingwin(xvar = spei_list_flag, # climate variables we are interested in
                         cdate = clim_daily_flag$cdate, #climate date
                         bdate = cover_flag$bdate, #biological date
                         baseline = lm(cover ~ 1, data = cover_flag), # can do anything here whcih is kinda cool but i think we just keep it simple
                         cinterval = "month", # no clue why i needed to put it into daily format if i can pick this but it wouldnt work otherwise
                         range = window_range, # set above..... idk how to make this decision
                         type = "relative",
                         stat = "mean",
                         func = "lin",
                         cmissing = FALSE) 

# find best SPEI scale and window from function above
# row index of model above with lowest AICc, whichmin is from climwin
best_row_flag <- which.min(sw_flag$combos$DeltaAICc)
best_scale_flag <- as.character(sw_flag$combos$climate[best_row_flag])
best_model_flag <- sw_flag[[best_row_flag]]
best_scale_flag

# run climwin's randomiztion test with best scale from above
?randwin
rand_flag <- randwin(repeats = n_rand,
                        window = "sliding",
                        xvar = list(spei = clim_daily_flag[[best_scale_flag]]),
                        cdate = clim_daily_flag$cdate,
                        bdate = cover_flag$bdate,
                        baseline = lm(cover ~ 1, data = cover_flag),
                        cinterval = "month",
                        range = window_range,
                        type = "relative",
                        stat = "mean",
                        func = "lin",
                        cmissing = FALSE)

# calc p value and plot diagnostics from climwin package
# computes p-value for AICc relative to randomized runs
p_aicc_flag <- as.numeric(pvalue(best_model_flag$Dataset, rand_flag[[1]], metric = "AIC"))
# safeguards against missing numeric conversions by defaulting to 0
# idk from the vignette
if (is.na(p_aicc_flag)) p_aicc_flag <- 0
# calculates randomization probability metric Pc, also from vignette
p_c_flag <- pvalue(best_model_flag$Dataset, rand_flag[[1]], metric = "C", sample.size = nrow(cover_flag))

# these are the output diagnostics from the climwin package and from the paper
pdf("data/climwin_output/climwin_plots_flagstaff.pdf", width = 12, height = 8)
plotall(dataset= best_model_flag$Dataset,
        datasetrand = rand_flag[[1]],
        bestmodel = best_model_flag$BestModel,
        bestmodeldata = best_model_flag$BestModelData)
dev.off()

# final summary for jornada climwin for lab meeting
results_flagstaff <- data.frame(site = "flagstaff",
                              spei_scale = best_scale_flag,
                              window_open = best_model_flag$Dataset$WindowOpen[1],
                              window_close = best_model_flag$Dataset$WindowClose[1],
                              delta_aicc = best_model_flag$Dataset$deltaAICc[1],
                              slope = best_model_flag$Dataset$ModelBeta[1],
                              adj_r2 = summary(best_model_flag$BestModel)$adj.r.squared,
                              p_delta_aicc = p_aicc_flag,
                              p_c = p_c_flag,
                              significant = p_aicc_flag < 0.05,
                              n_years = nrow(cover_flag),
                              first_year = min(cover_flag$year),
                              last_year = max(cover_flag$year))
results_flagstaff


# SRER
#combine climate and cover
clim_srer <- spei_site |> 
  filter(site == "santarita") |> 
  drop_na(starts_with("spei"))

# climwin needs daily inputs? trying to expand monthly data in my df
clim_daily_srer <- data.frame(date = seq(min(clim_srer$date), max(clim_srer$date), by = "day")) |> 
  mutate(month = floor_date(date, "month"),
         cdate = format(date, "%d/%m/%Y")) |> 
  left_join(clim_srer, by = c("month" = "date"))

# filter plant cover data from big df
min_clim_date_srer <- min(clim_srer$date) # this is bcus there is missing SPEI at the beg of the record

cover_srer <- cover_site |> 
  filter(site == "santarita", keep == TRUE,
         survey_date >= (min_clim_date_srer %m+% months(window_range[1])), 
         survey_date <= max(clim_daily_srer$date)) |> 
  mutate(bdate = format(survey_date, "%d/%m/%Y"))

# ok finally use climwin, at line 300 lolz
spei_list_srer <- as.list(select(clim_daily_srer, starts_with("spei")))

?slidingwin 
sw_srer <- slidingwin(xvar = spei_list_srer, # climate variables we are interested in
                      cdate = clim_daily_srer$cdate, #climate date
                      bdate = cover_srer$bdate, #biological date
                      baseline = lm(cover ~ 1, data = cover_srer), # can do anything here whcih is kinda cool but i think we just keep it simple
                      cinterval = "month", # no clue why i needed to put it into daily format if i can pick this but it wouldnt work otherwise
                      range = window_range, # set above..... idk how to make this decision
                      type = "relative",
                      stat = "mean",
                      func = "lin",
                      cmissing = FALSE) 

# find best SPEI scale and window from function above
# row index of model above with lowest AICc, whichmin is from climwin
best_row_srer <- which.min(sw_srer$combos$DeltaAICc)
best_scale_srer <- as.character(sw_srer$combos$climate[best_row_srer])
best_model_srer <- sw_srer[[best_row_srer]]
best_scale_srer

# run climwin's randomiztion test with best scale from above
?randwin
rand_srer <- randwin(repeats = n_rand,
                     window = "sliding",
                     xvar = list(spei = clim_daily_srer[[best_scale_srer]]),
                     cdate = clim_daily_srer$cdate,
                     bdate = cover_srer$bdate,
                     baseline = lm(cover ~ 1, data = cover_srer),
                     cinterval = "month",
                     range = window_range,
                     type = "relative",
                     stat = "mean",
                     func = "lin",
                     cmissing = FALSE)

# calc p value and plot diagnostics from climwin package
# computes p-value for AICc relative to randomized runs
p_aicc_srer <- as.numeric(pvalue(best_model_srer$Dataset, rand_srer[[1]], metric = "AIC"))
# safeguards against missing numeric conversions by defaulting to 0
# idk from the vignette
if (is.na(p_aicc_srer)) p_aicc_srer <- 0
# calculates randomization probability metric Pc, also from vignette
p_c_srer <- pvalue(best_model_srer$Dataset, rand_srer[[1]], metric = "C", sample.size = nrow(cover_srer))

# these are the output diagnostics from the climwin package and from the paper
pdf("data/climwin_output/climwin_plots_santarita.pdf", width = 12, height = 8)
plotall(dataset= best_model_srer$Dataset,
        datasetrand = rand_srer[[1]],
        bestmodel = best_model_srer$BestModel,
        bestmodeldata = best_model_srer$BestModelData)
dev.off()

# final summary for jornada climwin for lab meeting
results_santarita <- data.frame(site = "santarita",
                                spei_scale = best_scale_srer,
                                window_open = best_model_srer$Dataset$WindowOpen[1],
                                window_close = best_model_srer$Dataset$WindowClose[1],
                                delta_aicc = best_model_srer$Dataset$deltaAICc[1],
                                slope = best_model_srer$Dataset$ModelBeta[1],
                                adj_r2 = summary(best_model_srer$BestModel)$adj.r.squared,
                                p_delta_aicc = p_aicc_srer,
                                p_c = p_c_srer,
                                significant = p_aicc_srer < 0.05,
                                n_years = nrow(cover_srer),
                                first_year = min(cover_srer$year),
                                last_year = max(cover_srer$year))
results_santarita

# jornada = spei36 over 11 month window best predicts cover
# no significant windows found for flag or srer :/
# was this useless? PERHAPS.



# okay try with the dredge function, try with MAP, MAT, PDSI, VPD
# which climate variable best predicts total plant cover at each site?
# dredge version: SPEI (all scales), MAT, MAP, VPD, PDSI
#install.packages(c("MuMIn", "scPDSI"))
library(tidyverse)
library(lubridate)
library(prism)
library(terra)
library(MuMIn)
library(scPDSI)
library(nlme)



