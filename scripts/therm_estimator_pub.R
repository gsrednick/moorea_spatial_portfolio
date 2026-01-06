# Temperature curation ####

library(tidyverse)
library(vegan)
library(ggforce)
library(stats)
library(ggnewscale)
library(lubridate)
library(parallel)
library(DHARMa)
library(ggpmisc)
library(data.table)
library(stringr)

se_fn <- function(x, na.rm = FALSE) {
  if (na.rm) x <- na.omit(x)
  sd(x) / sqrt(length(x))
}

# Derived temperature predictor variables for analyses
# DTR - diurnal/daily temperature range
# DHD - degree heating days

# Data attributes
# Fringe data are 8-minute intervals 2005-present
# Backreef + Forereef are 20-minute interpolations until 2021, then 2-minute "raw" data


# Bring in data
temp_path <- "./data/environmental/temperature/MAIN" # download all csv files (n = 12 different csvs) from knb-lter-mcr.1035 to this folder, except for LTER00. See below.
temp_files <- list.files(path = temp_path, pattern = "\\.csv$", full.names = TRUE)
temp_files

combined_temp <- data.frame()

for (file in temp_files) {
  data <- read.csv(file, header = TRUE)  # Adjust options as per your CSV file

  combined_temp <- rbind(combined_temp, data)
}


# Upload LTER00 temperature data here
LTER00_data<-read.csv('./data/environmental/temperature/LTER00/MCR_LTER00_BTM_Backreef_Forereef_20250522.csv')

unique(combined_temp$reef_type_code)

LTER00_data_forjoin<-LTER00_data %>%
  mutate(site = str_replace(site,"_","0")) %>%
  mutate(site = recode(site, "LTER00" = "LTER01"))


combined_temp_siteupd<-combined_temp %>%
  mutate(site = str_replace(site,"_","0"))# %>%

combined_temp_upd<-rbind(combined_temp_siteupd,LTER00_data_forjoin)

# Time formatting
combined_temp_upd <- combined_temp_upd %>%
  mutate(
    time_local_clean = ifelse(
      str_detect(time_local, "^\\d{4}-\\d{2}-\\d{2}$"),
      paste0(time_local, " 00:00:00"),
      time_local
    ),
    time_local_clean = ymd_hms(time_local_clean)  # lubridate parsing
  )


combined_temp_upd<-combined_temp_upd %>%
  mutate(time_local = time_local_clean) %>%
  dplyr::select(-time_local_clean)



# Average 20min resolution to 1 hour resolution
dt <- as.data.table(combined_temp_upd)
remove(combined_temp_upd) # remove to clean the space

dt[, `:=`(
  hour = format(time_local, "%H"),
  min = format(time_local, "%M"),
  day = format(time_local, "%d"),
  month = format(time_local, "%m"),
  year = format(time_local, "%Y")
)]

# Group and summarize quickly
combined_temp_corrected_hour <- dt[
  , .(temperature_c = mean(temperature_c, na.rm = TRUE)),
  by = .(site, reef_type_code, sensor_depth_m, sensor_type, year, month, day, hour)
]


# Changes to sampling interval over time series
# Starts at 20 minute

combined_temp_corrected<-combined_temp_corrected_hour %>%
  mutate(time_local = as.POSIXct(paste(year, month, day, hour, "00:00"), format = "%Y %m %d %H %M"))


# Rename levels for proper merge later
dt <- as.data.table(combined_temp_corrected)
unique(combined_temp_corrected$reef_type_code)


# Rename column
setnames(dt, "site", "Site")

# Recode `reef_type_code` to `Habitat`
dt[, Habitat := fcase(
  reef_type_code == "Backreef", "Backreef",
  reef_type_code == "Forereef", "Outer",
  reef_type_code == "Fringing", "Fringe",
  default = NA_character_
)]

# Add combined columns
dt[, Habitat_depth := paste0(Habitat, "_", sensor_depth_m)]
dt[, site_hab_depth := paste0(Site, "_", Habitat_depth)]

# Store result
combined_temp_upd_mhw <- dt

# All loggers fringe
combined_temp_upd_mhw %>%
  filter(Habitat == "Fringe") %>%
  ggplot(aes(x = time_local, y= temperature_c)) +
  geom_point(size = 0.2) +
  facet_grid(Site~interaction(Habitat,sensor_depth_m)) +
  theme_bw()

# All loggers backreef
combined_temp_upd_mhw %>%
  filter(Habitat == "Backreef") %>%
  ggplot(aes(x = time_local, y= temperature_c)) +
  geom_point(size = 0.2) +
  facet_grid(Site~interaction(Habitat,sensor_depth_m)) +
  theme_bw()

# All loggers outer
combined_temp_upd_mhw %>%
  filter(Habitat == "Outer") %>%
  ggplot(aes(x = time_local, y= temperature_c)) +
  geom_point(size = 0.2) +
  facet_grid(Site~interaction(Habitat,sensor_depth_m)) +
  theme_bw()

combined_temp_upd_mhw %>%
  group_by(Site,Habitat) %>%
  dplyr::summarise(n = n())

# Select depths of interest
combined_temp_upd_mhw_ready<-combined_temp_upd_mhw %>%
  filter(sensor_depth_m < 30,
         !(Habitat == "Fringe" & sensor_depth_m >= 3)) # this is the fringe thermistor of interest



# Steps - Fixing missing data ####
# Will predict missing data based on adjacent sites
# (1) Find specific gaps
# (2) Generate regressions for adjacent sites; verify relationships
# (3) Fill values with predicted
# (4) Fill remaining will spatially aggregated mean


## Filling date sequence ####
# Starting and ending date-times
combined_temp_upd_mhw_ready <- combined_temp_upd_mhw_ready %>%
  mutate(time_local = force_tz(time_local, tzone = "Pacific/Tahiti"))

temp_start <- force_tz(as.POSIXct("2005-12-10 00:00:00"), tzone = "Pacific/Tahiti")
temp_end <- force_tz(as.POSIXct(max(combined_temp_upd_mhw_ready$time_local, na.rm = TRUE)), tzone = "Pacific/Tahiti")
full_date_seq <- seq.POSIXt(temp_start, temp_end, by = "hour", tz = "Pacific/Tahiti")

length(full_date_seq)
length(unique(full_date_seq))

group_vars <- unique(combined_temp_upd_mhw_ready$site_hab_depth)
length(group_vars)
length(unique(group_vars))

# empty list to store each chunk
list_chunks <- list()

# Loop through each combination and generate the full sequence separately -- necessary because the dataset is so big
for(i in 1:length(group_vars)) {
  # Subset the data for this site_hab_depth combination
  current_site_hab_depth <- group_vars[i]

  # Create a subset of the data for this site_hab_depth
  temp_subset <- expand.grid(
    site_hab_depth = current_site_hab_depth,
    time_local = full_date_seq
  )

  # Add the combination data to the list
  list_chunks[[i]] <- temp_subset
}

# Combine all site/habitats together
full_seq_df_pre <- do.call(rbind, list_chunks)
dim(full_seq_df_pre)
any(duplicated(full_seq_df_pre))
full_seq_df<-full_seq_df_pre
dim(full_seq_df)


full_seq_lab_df <- full_seq_df %>%
  mutate(
  ) %>%
  separate(site_hab_depth, into = c("Site", "Habitat", "Depth"), sep = "_", remove = FALSE) %>%
  mutate(
    Habitat = case_when(
      str_detect(site_hab_depth, "10") ~ "10m",
      str_detect(site_hab_depth, "20") ~ "17m",
      TRUE ~ Habitat
    )
  ) %>%
  select(-Depth)

combined_temp_upd_mhw_ready_actual <- combined_temp_upd_mhw_ready %>%
  mutate(
  ) %>%
  separate(site_hab_depth, into = c("Site", "Habitat", "Depth"), sep = "_", remove = FALSE) %>%
  mutate(
    Habitat = case_when(
      str_detect(site_hab_depth, "10") ~ "10m",
      str_detect(site_hab_depth, "20") ~ "17m",
      TRUE ~ Habitat
    )
  ) %>%
  select(-c(Depth, Habitat_depth, reef_type_code, sensor_depth_m, sensor_type))

# Join and complete
combined_temp_upd_mhw_ready_seq <- full_seq_lab_df %>%
  left_join(combined_temp_upd_mhw_ready_actual,
            by = c("Site", "Habitat", "site_hab_depth", "time_local"#, "year", "month", "day", "hour"
                   ))
View(combined_temp_upd_mhw_ready_seq)


dim(full_seq_lab_df) # [1] 3913704      13
dim(combined_temp_upd_mhw_ready_actual) # [1] 3398168      13
dim(combined_temp_upd_mhw_ready_seq) # [1] 3913704      13

# Extract elements from datetime
combined_temp_upd_mhw_ready_seq<-combined_temp_upd_mhw_ready_seq %>%
  mutate(
    day = sprintf("%02d", day(time_local)),
    month = sprintf("%02d", month(time_local)),
    year = as.character(year(time_local)),
    hour = sprintf("%02d", hour(time_local)))


# look for reps
temp_sampling<-combined_temp_upd_mhw_ready_seq %>%
  group_by(Site,Habitat,year) %>%
  summarize(count_samp_pres = sum(!is.na(temperature_c)),
            count_samp_NA = sum(is.na(temperature_c)))

# This figure shows how much data are missing per site/habitat
temp_sampling %>%
  filter(year > 2008) %>%
  #na.omit() %>%
  ggplot(aes(x = year, y = count_samp_pres)) +
  geom_hline(yintercept = 8760/2,linetype = "dotted") +
  geom_point(size = 3, color = "blue") +
  geom_point(aes(y = count_samp_NA), size = 2, color = "red", alpha = 0.8) +
  facet_grid(Site~Habitat) +
  theme_bw()

max(temp_sampling$count_samp_pres) # there cant be more than 8760 hours per year; but there are more hours because our sampling isnt exactly constrained for each year
min(temp_sampling$count_samp_pres)
min(temp_sampling$count_samp_NA)



no_missing<-combined_temp_upd_mhw_ready_seq %>%
  group_by(Site,Habitat,year) %>%
  dplyr::summarize(count_NAs = sum(is.na(temperature_c)),
            count_obs = sum(!is.na(temperature_c)),
            total_obs = count_NAs + count_obs)

no_missing_totals <- no_missing %>%
  group_by(Site,Habitat) %>%
  dplyr::summarize(count_NAs_tot = sum(count_NAs),
            count_obs_tot = sum(count_obs),
            total_obs_tot = sum(total_obs),
            frac_missing = count_NAs_tot/total_obs_tot,
            percent_missing = frac_missing * 100)


no_missing_days<-no_missing_totals %>%
  mutate(missing_hours = count_NAs_tot/24,
         missing_days = count_NAs_tot/8760)



## Regressions for predicting missing values
gap_summary <- combined_temp_upd_mhw_ready_seq %>%
  na.omit() %>%
  arrange(Site, Habitat, time_local) %>%  # safer to arrange all grouping vars first
  group_by(Site, Habitat) %>%
  mutate(
    time_diff = as.numeric(difftime(time_local, lag(time_local), units = "hours")),
    time_diff_days = time_diff / 24
  ) %>%
  #filter(!is.na(time_diff) & time_diff > gap_threshold) %>%
  dplyr::summarise(
    mean_gap_days = mean(time_diff_days, na.rm = TRUE),
    se_gap_days = sd(time_diff_days, na.rm = TRUE) / sqrt(n()),
    max_gap_days = max(time_diff_days, na.rm = TRUE),
    n_gaps = n(),
    .groups = "drop"
  )


gap_summary_overall <- combined_temp_upd_mhw_ready_seq %>%
  na.omit() %>%
  arrange(Site, Habitat, time_local) %>%  # safer to arrange all grouping vars first
  mutate(
    time_diff = as.numeric(difftime(time_local, lag(time_local), units = "hours")),
    time_diff_days = time_diff / 24
  ) %>%
  dplyr::summarise(
    mean_gap_days = mean(time_diff_days, na.rm = TRUE),
    se_gap_days = sd(time_diff_days, na.rm = TRUE) / sqrt(n()),
    max_gap_days = max(time_diff_days, na.rm = TRUE),
    n_gaps = n(),
    .groups = "drop"
  )


gap_threshold <- (3600 * 24)*14  # threshold for gaps in seconds (e.g., 1 hour (3600) for hourly data) -- this is currently 2 weeks of missing values

temporal_gaps <- combined_temp_upd_mhw_ready_seq %>%
  na.omit() %>%
  group_by(Site,Habitat) %>%
  arrange(time_local) %>%
  mutate(
    time_diff = as.numeric(difftime(time_local, lag(time_local), units = "secs")),
    start_gap = lag(time_local),
    end_gap = time_local
  ) %>%
  filter(!is.na(start_gap) & time_diff > gap_threshold) %>%
  dplyr::select(Site, Habitat, start_gap, end_gap) %>%
  tidyr::pivot_longer(cols = c(start_gap, end_gap), names_to = "gap_type", values_to = "gap_time")

temporal_gaps <- temporal_gaps %>%
  mutate(line_type = ifelse(gap_type == "start_gap", "start", "end"))

temp_start_times<-combined_temp_upd_mhw_ready %>%
    mutate(time_local = as.POSIXct(time_local,format = "%Y-%m-%d %H:%M:%S")) %>%
    group_by(Site,reef_type_code) %>%
    dplyr::summarize(first_samp= min(time_local, na.rm = T))



# Define shore groups
shore_groups <- list(
  North = c("LTER01", "LTER02"),
  East = c("LTER03", "LTER04"),
  West = c("LTER05", "LTER06")
)

## Linear model for fit between backreef and fringe ####
nonfilled_temp_df<-combined_temp_upd_mhw_ready_seq %>%
  mutate(Site_hab = paste0(Site,"_",Habitat))  %>%
  select(Site_hab,time_local,temperature_c) %>%
  pivot_wider(id_cols = time_local, names_from = "Site_hab", values_from = "temperature_c")


## Looped for every site
sites = unique(combined_temp_upd_mhw_ready_seq$Site)
fringe_predicted_list <- list()

for (site in sites) {
  backreef_col <- paste0(site, "_Backreef")
  fringe_col <- paste0(site, "_Fringe")

  # Filter rows where both are present
  temp_df <- nonfilled_temp_df %>%
    select(time_local, backreef_col, fringe_col) %>%
    rename(Backreef = backreef_col,
           Fringe = fringe_col)

  # Fit model where both are present
  fit <- lm(Fringe ~ Backreef, data = temp_df %>% filter(!is.na(Fringe) & !is.na(Backreef)))

  # Predict full vector & get residuals
  pred_df <- temp_df %>%
    mutate(
      pred_fringe = predict(fit, newdata = temp_df),
      resid_fringe = Fringe - pred_fringe
    )

  resid_sd <- sd(pred_df$resid_fringe, na.rm = TRUE)

  # Final prediction with noise for missing only
  pred_df <- pred_df %>%
    mutate(
      predicted_fringe_missing = if_else(
        is.na(Fringe),
        predict(fit, newdata = pred_df) + rnorm(n(), mean = 0, sd = resid_sd),
        Fringe
      )
    ) %>%
    select(time_local, predicted_fringe_missing) %>%
    rename(!!paste0(site) := predicted_fringe_missing)

  fringe_predicted_list[[site]] <- pred_df
}

fringe_filled_df <- reduce(fringe_predicted_list, full_join, by = "time_local")

fringe_filled_df_long<-fringe_filled_df %>%
  pivot_longer(-time_local, names_to = "Site", values_to = "lm_filled_temp") %>%
  mutate(Habitat = "Fringe")


fringe_full_merge_filled_data<-merge(fringe_filled_df_long,combined_temp_upd_mhw_ready_seq, by = c("Site", "Habitat", "time_local")) #%>%
dim(combined_temp_upd_mhw_ready_seq)
dim(fringe_filled_df_long)
dim(fringe_full_merge_filled_data)


## Now do the same thing for backreef ####
backreef_predicted_list <- list()
for (site in sites) {
  backreef_col <- paste0(site, "_Backreef")
  fringe_col <- paste0(site, "_Fringe")

  # Filter rows where both are present
  temp_df <- nonfilled_temp_df %>%
    select(time_local, backreef_col, fringe_col) %>%
    rename(Backreef = backreef_col,
           Fringe = fringe_col)

  # Fit model where both are present
  fit <- lm(Backreef ~ Fringe, data = temp_df %>% filter(!is.na(Fringe) & !is.na(Backreef)))

  # Predict full vector & get residuals
  pred_df <- temp_df %>%
    mutate(
      pred_backreef = predict(fit, newdata = temp_df),
      resid_backreef = Backreef - pred_backreef
    )

  resid_sd <- sd(pred_df$resid_backreef, na.rm = TRUE)

  # Final prediction with noise for missing only
  pred_df <- pred_df %>%
    mutate(
      predicted_backreef_missing = if_else(
        is.na(Backreef),
        predict(fit, newdata = pred_df) + rnorm(n(), mean = 0, sd = resid_sd),
        Backreef
      )
    ) %>%
    select(time_local, predicted_backreef_missing) %>%
    rename(!!paste0(site) := predicted_backreef_missing)

  backreef_predicted_list[[site]] <- pred_df
}

backreef_filled_df <- reduce(backreef_predicted_list, full_join, by = "time_local")

backreef_filled_df_long<-backreef_filled_df %>%
  pivot_longer(-time_local, names_to = "Site", values_to = "lm_filled_temp") %>%
  mutate(Habitat = "Backreef")

backreef_full_merge_filled_data<-merge(backreef_filled_df_long,combined_temp_upd_mhw_ready_seq, by = c("Site", "Habitat", "time_local")) #%>%
dim(combined_temp_upd_mhw_ready_seq)
dim(backreef_filled_df_long)
dim(backreef_full_merge_filled_data)



## Now 10m sites ####
m10_predicted_list <- list()
for (site in sites) {
  m10_col <- paste0(site, "_10m")
  m17_col <- paste0(site, "_17m")

  # Filter rows where both are present
  temp_df <- nonfilled_temp_df %>%
    select(time_local, m10_col, m17_col) %>%
    rename("10m" = m10_col,
           "17m" = m17_col)

  # Fit model where both are present
  fit <- lm(`10m` ~ `17m`, data = temp_df %>% filter(!is.na(`10m`) & !is.na(`17m`)))

  # Predict full vector & get residuals
  pred_df <- temp_df %>%
    mutate(
      pred_10m = predict(fit, newdata = temp_df),
      resid_10m = `10m` - pred_10m
    )

  resid_sd <- sd(pred_df$resid_10m, na.rm = TRUE)

  # Final prediction with noise for missing only
  pred_df <- pred_df %>%
    mutate(
      predicted_10m_missing = if_else(
        is.na(`10m`),
        predict(fit, newdata = pred_df) + rnorm(n(), mean = 0, sd = resid_sd),
        `10m`
      )
    ) %>%
    select(time_local, predicted_10m_missing) %>%
    rename(!!paste0(site) := predicted_10m_missing)

  m10_predicted_list[[site]] <- pred_df
}

m10_filled_df <- reduce(m10_predicted_list, full_join, by = "time_local")

m10_filled_df_long<-m10_filled_df %>%
  pivot_longer(-time_local, names_to = "Site", values_to = "lm_filled_temp") %>%
  mutate(Habitat = "10m")

m10_full_merge_filled_data<-merge(m10_filled_df_long,combined_temp_upd_mhw_ready_seq, by = c("Site", "Habitat", "time_local")) #%>%
dim(combined_temp_upd_mhw_ready_seq)
dim(m10_filled_df_long)
dim(m10_full_merge_filled_data)


## Now 17m sites ####
m17_predicted_list <- list()
for (site in sites) {
  m10_col <- paste0(site, "_10m")
  m17_col <- paste0(site, "_17m")

  # Filter rows where both are present
  temp_df <- nonfilled_temp_df %>%
    select(time_local, m10_col, m17_col) %>%
    rename("10m" = m10_col,
           "17m" = m17_col)

  # Fit model where both are present
  fit <- lm(`17m` ~ `10m`, data = temp_df %>% filter(!is.na(`10m`) & !is.na(`17m`)))

  # Predict full vector & get residuals
  pred_df <- temp_df %>%
    mutate(
      pred_17m = predict(fit, newdata = temp_df),
      resid_17m = `17m` - pred_17m
    )

  resid_sd <- sd(pred_df$resid_17m, na.rm = TRUE)

  # Final prediction with noise for missing only
  pred_df <- pred_df %>%
    mutate(
      predicted_17m_missing = if_else(
        is.na(`17m`),
        predict(fit, newdata = pred_df) + rnorm(n(), mean = 0, sd = resid_sd),
        `17m`
      )
    ) %>%
    select(time_local, predicted_17m_missing) %>%
    rename(!!paste0(site) := predicted_17m_missing)

  m17_predicted_list[[site]] <- pred_df
}

m17_filled_df <- reduce(m17_predicted_list, full_join, by = "time_local")

m17_filled_df_long<-m17_filled_df %>%
  pivot_longer(-time_local, names_to = "Site", values_to = "lm_filled_temp") %>%
  mutate(Habitat = "17m")

m17_full_merge_filled_data<-merge(m17_filled_df_long,combined_temp_upd_mhw_ready_seq, by = c("Site", "Habitat", "time_local")) #%>%
dim(combined_temp_upd_mhw_ready_seq)
dim(m17_filled_df_long)
dim(m17_full_merge_filled_data)


# Make dfs ready for rbind
fringe_filled_for_merge<-fringe_full_merge_filled_data %>% select(time_local,Site,Habitat,temperature_c,lm_filled_temp)
backreef_filled_for_merge<-backreef_full_merge_filled_data %>% select(time_local,Site,Habitat,temperature_c,lm_filled_temp)
m10_filled_for_merge<-m10_full_merge_filled_data %>% select(time_local,Site,Habitat,temperature_c,lm_filled_temp)
m17_filled_for_merge<-m17_full_merge_filled_data %>% select(time_local,Site,Habitat,temperature_c,lm_filled_temp)


# Bind
combined_habs_sites_forfill<-rbind(
  fringe_filled_for_merge,
  backreef_filled_for_merge,
  m10_filled_for_merge,
  m17_filled_for_merge)



# Fill remaining gaps with seasonal averaging procedure ####

fill_gaps_combined_optimized_2ndfill <- function(df, full_data) {
  # Step 1: Calculate seasonal mean for the site of interest
  df <- df %>%
    mutate(doy = yday(time_local),
           hour = hour(time_local),
           year = year(time_local)) %>%
    group_by(year,doy, hour) %>%
    mutate(seasonal_mean = mean(lm_filled_temp, na.rm = TRUE)) %>%
    ungroup()

  # Step 2: Precompute values for related sites
  shore_group <- names(shore_groups)[sapply(shore_groups, function(group) unique(df$Site) %in% group)]
  related_sites <- shore_groups[[shore_group]]

  related_values <- full_data %>%
    filter(Site %in% related_sites & !is.na(lm_filled_temp)) %>%
    mutate(doy = yday(time_local),
           hour = hour(time_local),
           year = year(time_local)) %>%
    group_by(year,doy, hour) %>%
    summarise(related_real_value = mean(lm_filled_temp, na.rm = TRUE), .groups = "drop")

  # Join related values back into the original dataframe
  df <- df %>%
    left_join(related_values, by = c("year", "doy", "hour"))

  # Step 3: Combine seasonal mean and real values from related sites
  df <- df %>%
    mutate(
      temperature_c_filled = if_else(
        is.na(lm_filled_temp),
        rowMeans(cbind(seasonal_mean, related_real_value), na.rm = TRUE),
        lm_filled_temp
      )
    ) %>%
    select(#-seasonal_mean, -related_real_value,
      -doy, -hour)

  return(df)
}


prefilled_grouped_data <- combined_habs_sites_forfill %>%
  mutate(hour = hour(time_local),
         month = month(time_local),
         day = day(time_local),
         year = year(time_local)) %>%
  dplyr::select(Site,Habitat,time_local,day,month,year,hour,temperature_c,lm_filled_temp) %>%
  group_by(Site, Habitat) %>%
  group_split()

# Fill data in parallel
gc() # garbage collection to speed things up
gc() # again for extra ummmph
gc() # again for even more ummmph

second_fill_filled_data <- mclapply(
  prefilled_grouped_data,
  function(df) fill_gaps_combined_optimized_2ndfill(df, combined_habs_sites_forfill),
  mc.cores = 8) %>%
  dplyr::bind_rows()



second_fill_filled_data_adjust<-second_fill_filled_data %>%
  mutate(adjusted_mean = rowMeans(cbind(seasonal_mean, related_real_value), na.rm = TRUE)) # this gives the "filled" value for all cases where NA is present -- i.e., omits NA


filled_ts_plot_wgaps_2ndfill<-
  ggplot() +
  geom_rect(data = temporal_gaps, inherit.aes = F,
            aes(xmin = lag(gap_time), xmax = gap_time,
                ymin = -Inf, ymax = Inf),fill = "grey", alpha = 0.6) +
  geom_point(data = second_fill_filled_data_adjust %>% mutate(Source = "Mean filling"), aes(x = time_local, y = temperature_c_filled,color = Source), size = 0.2) +
  geom_point(data = second_fill_filled_data_adjust %>% mutate(Source = "LM prediction filled"), aes(x = time_local, y = lm_filled_temp,color = Source), alpha = 0.2, size = 0.2) + # this still isnt wuite right. still need to adjust prediction; the first mean might be too much of a reduction in variance; maybe adjust by variance
  geom_point(data = second_fill_filled_data_adjust %>% mutate(Source = "Original hourly"), aes(x = time_local, y = temperature_c, color = Source), size = 0.2,alpha = 0.01) +
  theme_bw() +
  scale_color_manual(name = "Data Source",
                     values = c("Mean filling" = "red",
                                "LM prediction filled" = "blue",
                                "Original hourly" = "green")) +
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +  # Rotate x-axis labels
  facet_grid(Site~Habitat) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  labs(x = "Time (hourly)", y = "Temp (hourly mean C)", title = "Temperature time series: >14d gap shown in grey")


ggsave('./data/environmental/summarized/Figure_S1_temperature_timeseries.jpeg',
       filled_ts_plot_wgaps_2ndfill,
       dpi = 300,
       height = 10,
       width = 16)

# Export filled data to csv
write.csv(second_fill_filled_data_adjust,'./data/environmental/summarized/interpolated_temperature_dat_JUN25_2ndfill_upd.csv')


# Check residuals ####

# Look at mean and sd of residuals for each fit
temp_data_split <- second_fill_filled_data_adjust %>%
  filter(!is.na(temperature_c)) %>%
  group_by(Site, Habitat) %>%
  group_split()

# Function to run diagnostics for each group
run_diagnostics <- function(df) {
  site <- unique(df$Site)
  habitat <- unique(df$Habitat)

  if (nrow(df) < 5) return(NULL)  # skip if too few rows

  model <- lm(adjusted_mean ~ temperature_c, data = df)
  sim <- simulateResiduals(model, plot = F)

  tibble(
    Site = site,
    Habitat = habitat,
    SD_resid = sd(sim$scaledResiduals),
    var_resid = var(sim$scaledResiduals),
    Mean_resid = mean(sim$scaledResiduals)
  )
}

# Run function on all group splits
temp_interp_summary_df <- map_dfr(temp_data_split, run_diagnostics)


temp_interp_summary_df_summarized<-temp_interp_summary_df %>%
  #mutate(mean_var = SD_resid^2) %>%
  dplyr::summarize(u_mean_var = mean(var_resid),
                   sd(var_resid),
                   se = se_fn(var_resid))


filled_data_fits_plot_2nd<-second_fill_filled_data_adjust %>%
  filter(!is.na(temperature_c)) %>%
  ggplot(aes(x = temperature_c,y = adjusted_mean)) +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  geom_point(size = 0.2) +
  stat_poly_line() +
  facet_grid(Site~Habitat) +
  theme_bw() +
  labs(x = "Temperature (°C)", y= "Interpolated temperature (°C)") +
  theme(axis.title.x = element_text(color = "green"),
        axis.title.y = element_text(color = "red")) +
  geom_text(data = temp_interp_summary_df,
            aes(
              x = 22, y = 31.5,
              label = paste0("Mean resid. = ", round(Mean_resid, 3),
                             "\nSD resid.= ", round(SD_resid, 3))
            ),
            hjust = 0, vjust = 1
  ) +
  coord_fixed(xlim = c(22,33), ylim = c(22,33))

ggsave('./data/environmental/summarized/Figure_S2_temperature_fits.jpeg',
       filled_data_fits_plot_2nd,
       dpi = 300,
       height = 16,
       width = 14
)



# DTR & DHD estimation ####
mmmt <- 28.8
samples_per_day <- 24 # hourly
window_size <- 12 * samples_per_day  # 12 days for hourly data

dhd_rollify <- tibbletime::rollify(sum, window = window_size)

mod_DHD <- second_fill_filled_data_adjust %>%
  mutate(
    time_local = as_datetime(time_local),
    time_hourly = floor_date(time_local, "1 hour")
  ) %>%
  group_by(Site, Habitat, time_hourly) %>%
  summarise(
    temperature_c = mean(temperature_c, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    hspt = temperature_c - mmmt,
    hsptm = ifelse(hspt >= 1, hspt, 0),
    dhd = dhd_rollify(hsptm) / samples_per_day,  # normalize
    year = year(time_hourly)
  )

mod_DHD_window <- mod_DHD

# Sum DHD per year, site, habitat
annual_DHD <- mod_DHD_window %>%
  group_by(Site, Habitat, year) %>%
  dplyr::summarise(DHD = sum(dhd, na.rm = TRUE)/24)


annual_DTR<-second_fill_filled_data_adjust %>%
  mutate(date = as.Date(time_local)) %>%
  filter(date > "2006-01-01" & date < "2024-06-30") %>%
  group_by(Site, Habitat, date) %>%
  dplyr::summarize(
    min_daily_temp = min(temperature_c_filled, na.rm = TRUE),
    max_daily_temp = max(temperature_c_filled, na.rm = TRUE),
    over29_8_flag = as.integer(any(temperature_c_filled > 29.8)),
    DH = sum(ifelse(temperature_c_filled >=29.8, temperature_c_filled - 28.8, 0))/24) %>%
  mutate(year = year(date),
         temp_range = max_daily_temp - min_daily_temp) %>% # calculate daily temperature range
  group_by(Site,Habitat,year) %>%
  dplyr::summarize(DTR = mean(temp_range, na.rm = TRUE),
                   Days_above_threshold = sum(over29_8_flag),
                   DHD_old = sum(DH, na.rm = TRUE)) %>%
  mutate(Days_above_threshold = ifelse(is.na(Days_above_threshold),0,Days_above_threshold))


# Output data ####
therm_dat_updated<-merge(annual_DTR,annual_DHD)
write.csv(therm_dat_updated, './data/environmental/summarized/thermal_predictors_updated.csv', row.names = F)


# END ####

