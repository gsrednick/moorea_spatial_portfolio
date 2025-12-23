# Wavelet clean script

# Packages
library(tidyverse)
library(stringr)
library(reshape2)
library(wsyn)
library(zoo)
library(patchwork)
library(ggExtra)
library(ggtext)


# Set seed for reproducibility
set.seed(42)

# Helper functions ####

# extract plotmag helper function ####
extract_plotmag_data <- function(wmpf_object, sigthresh = 0.95) {
  library(wsyn)

  # Extract magnitude, time, timescales, and significance
  wav <- Mod(get_values(wmpf_object))  # Magnitude
  times <- get_times(wmpf_object)  # Time (x-axis)
  timescales <- get_timescales(wmpf_object)  # Timescales (y-axis)
  signif <- get_signif(wmpf_object)  # Significance object

  # Create a dataframe for plotting
  plot_data <- expand.grid(Time = times, Timescale = timescales)
  plot_data$Magnitude <- as.vector(wav)

  # Initialize an empty dataframe for contours
  contour_df <- NULL

  # Extract significance contours if significance data exists
  if (!all(is.na(signif))) {
    if (signif[[1]] == "quick") {
      q <- stats::quantile(signif[[2]], sigthresh)
      contour_lines <- contourLines(
        x = times,
        y = log2(timescales),  # Convert timescales to log2
        z = wav,
        levels = q
      )
    } else if (signif[[1]] %in% c("fft", "aaft")) {
      contour_lines <- contourLines(
        x = times,
        y = log2(timescales),
        z = signif[[3]],
        levels = sigthresh
      )
    }

    # Convert contour lines to a dataframe
    if (!is.null(contour_lines) && length(contour_lines) > 0) {
      contour_df <- do.call(rbind, lapply(seq_along(contour_lines), function(i) {
        data.frame(
          Time = contour_lines[[i]]$x,
          Log2_Period = contour_lines[[i]]$y,
          Group = i  # Assign each contour segment a group ID
        )
      }))
    }
  }

  # Return a list with the extracted data
  return(list(
    plot_data = plot_data,
    contour_data = contour_df
  ))
}

extract_plotmag_wmf_data <- function(wmf_object) {
  library(wsyn)

  # Extract magnitude, time, timescales, and significance
  wav <- Mod(get_values(wmf_object))  # Magnitude
  times <- get_times(wmf_object)  # Time (x-axis)
  timescales <- get_timescales(wmf_object)  # Timescales (y-axis)

  # Create a dataframe for plotting
  plot_data <- expand.grid(Time = times, Timescale = timescales)
  plot_data$Magnitude <- as.vector(wav)

  # Return a list with the extracted data
  return(plot_data)
}



extract_plotmag_wmf_data <- function(wmf_object, sigthresh = 0.95) {
  library(wsyn)

  # Extract magnitude, time, timescales, and significance
  wav <- Mod(get_values(wmf_object))  # Magnitude
  real <- Re(get_values(wmf_object))  # Real part

  times <- get_times(wmf_object)  # Time (x-axis)
  timescales <- get_timescales(wmf_object)  # Timescales (y-axis)

  # Create a dataframe for plotting
  plot_data <- expand.grid(Time = times, Timescale = timescales)
  plot_data$Magnitude <- as.vector(wav)
  plot_data$Real <- as.vector(real)


  # Return a list with the extracted data
  return(plot_data)

}

# AR1 surrogate for generating WMF contours
generate_ar1_surrogate <- function(ts) {
  # Fit AR(1) model
  ar_fit <- tryCatch(ar(ts, aic = FALSE, order.max = 1), error = function(e) NULL)
  if (is.null(ar_fit)) return(rep(NA, length(ts)))
  phi <- ar_fit$ar
  arima.sim(model = list(ar = phi), n = length(ts))
}


# Data import ####

# Coral community structure data for species of interest -- from coral_data_processing.R
MCR_coral_data<-read.csv('./data/summarized/MCR_CTS_site_updated.csv')

# Summarized temperature data (DTR & DHD) -- from therm_estimator.R
thermal_conditions<-read.csv('./data/thermal_predictors_updated.csv') # this is with DHD estimated within year -- 2008-2024; Wyatt method, but not backwards in time

# Summarized algal data -- from algal_data_processing.R
biol_predictors<-read.csv('./data/summarized/MCR_alg_summarized.csv')

# Coordinates of each LTER site
mra_site_cords<-read.csv("./data/LTER_sites.csv")
mra_sites<-mra_site_cords %>% filter(Type == "survey", !Habitat == "Outer")


# Ocean distance between monitoring sites
ocean_dist<-read.csv('./data/ocean_dist_sites.csv')


# Environmental data curation ####
coral_years<-unique(MCR_coral_data$Date) # 2005-2007 temperature data are missing for the fringing reef. Using 2008 onward
length(coral_years) # 20 years shown here 2005-2024
years = seq(2006,max(coral_years)) # new, no 2005 data for backreef
length(years)

env_years<-unique(thermal_conditions$year) # 2005-2007 temperature data are missing for the fringing reef. Probably inappropriate to interpolate these here. Will have to use 2008 onward
length(env_years) # 17 years 2008-2024

env_years = 2008:2024


# Fill scale for all wavelet plots
#global_fill_scale <- scale_fill_viridis_c(limits = c(0, 1.4), option = "turbo", name = "Synchrony") # old, bad for Blue-Green Colorblind
global_fill_scale <- scale_fill_viridis_c(limits = c(0, 1.4), option = "magma", name = "Synchrony")

# 'wsyn' "cleandat" setting for coral data
bio_clev = 2


## Temperature data to matrix ####
DTR_matrix <- thermal_conditions %>%
  dplyr::filter(year %in% env_years) %>%
  dplyr::mutate(Site_hab = paste0(Site,"_",Habitat)) %>%
  dplyr::select(Site_hab, year, DTR) %>%
  pivot_wider(names_from = year, values_from = DTR) %>%
  column_to_rownames("Site_hab") %>%
  as.matrix()

DHD_matrix <- thermal_conditions %>%
  dplyr::filter(year %in% env_years) %>%
  dplyr::mutate(Site_hab = paste0(Site,"_",Habitat)) %>%
  dplyr::select(Site_hab, year, DHD) %>%
  pivot_wider(names_from = year, values_from = DHD) %>%
  column_to_rownames("Site_hab") %>%
  as.matrix()


## Algal data to matrix ####

alg_matrix <- biol_predictors %>%
  dplyr::mutate(Site = str_replace(Site,"_","0")) %>%
  dplyr::mutate(Site_hab = paste0(Site,"_",Habitat)) %>%
  dplyr::select(Site_hab, Year, mean_cover) %>%
  dplyr::group_by(Site_hab) %>%
  tidyr::complete(
    Year = env_years,  # Use explicit date sequence
    fill = list(mean_cover = NA)) %>%
  dplyr::filter(Year %in% env_years) %>%
  dplyr::mutate(mean_cover_approx = na.approx(mean_cover, rule = 2)) %>%
  dplyr::select(Site_hab, Year, mean_cover_approx) %>%
  pivot_wider(names_from = Year, values_from = mean_cover_approx) %>%
  column_to_rownames("Site_hab") %>%
  as.matrix()



## Clean predictors ####

dtr_clean<-cleandat(DTR_matrix, times = env_years, clev = 4) # removes mean, detrends, standardizes variance, and Box-Cox transforms
DHD_clean<-cleandat(DHD_matrix, times = env_years, clev = 4) # removes mean, detrends, standardizes variance, and Box-Cox transforms
alg_clean<-cleandat(alg_matrix, times = env_years, clev = 4) # removes mean, detrends, standardizes variance, and Box-Cox transforms

dtr_clean_cdat<-dtr_clean$cdat %>%
  as.data.frame() %>%
  as.matrix()

DHD_clean_cdat<-DHD_clean$cdat %>%
  as.data.frame() %>%
  as.matrix()

alg_clean_cdat<-alg_clean$cdat %>%
  as.data.frame() %>%
  as.matrix()


# Wavelet analysis | Corals ####

## Pocillopora ####
poc_df<-MCR_coral_data %>%
  dplyr::select(1:3,Pocillopora) %>%
  dplyr::mutate(Site = str_replace(Site," ","0")) %>%
  dplyr::mutate(Site_hab = paste0(Site,"_",Habitat)) %>%
  dplyr::select(-c(Site,Habitat))

startyear = min(years)
endyear = max(years)

poc_df_filled <- poc_df %>%
  filter(Date %in% years) %>%
  dplyr::group_by(Site_hab)

poc_df_wide<-poc_df_filled %>%
  pivot_wider(values_from = Pocillopora,names_from = Date) %>%
  as.data.frame()

row.names(poc_df_wide)<-poc_df_wide$Site_hab
poc_df_wide$Site_hab<-NULL
poc_mat<-as.matrix(poc_df_wide)

### Wavelet mean fields ####
# Wavelet Mean Field (WMF) for full timeseries (Figure 3)
poc_df_clean<-cleandat(poc_mat, times = years, clev = bio_clev) # removes mean, detrends, and standardizes variance
poc_cdat<-as.vector(poc_df_clean$cdat)
poc_fres_wmf<-wmf(poc_df_clean$cdat,years, f0 = 0.5) # for Figure 3 -- non-truncated time series

# Wavelet Mean Field (WMF) for analysis of environmental drivers of synchrony
poc_mat_short <- poc_mat[, colnames(poc_mat) %in% env_years]
poc_df_short_clean<-cleandat(poc_mat_short, times = env_years, clev = bio_clev) # removes mean, detrends, and standardizes variance
poc_fres_short_wmf<-wmf(poc_df_short_clean$cdat,env_years, f0 = 0.5) # for rest of analyses that use truncated time series


poc_detrended<-as.data.frame(poc_df_clean$cdat) %>%
  mutate(Site_hab = rownames(.)) %>%
  pivot_longer(-Site_hab, names_to = "Year",values_to = "cover") %>%
  separate(Site_hab, into = c("Site","Habitat"),sep = "_") %>%
  mutate(Year = as.numeric(Year))

# Plotting means and var for Figure 3 - modifed
poc_ts_plot_df<-poc_detrended %>%
  group_by(Year) %>%
  #filter(Date > 2007 & Date < 2020) %>%
  dplyr::summarise(mean = mean(cover,na.rm = T),
                   var = var(cover, na.rm = T),
                   sd = sd(cover, na.rm = T),
                   n = n(),
                   se = sd/sqrt(n))

poc_ts_plot<-poc_ts_plot_df %>%
  ggplot(aes(x = Year , y = mean)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(ymin = mean - se, ymax  = mean + se)) +
  geom_line() +
  theme_bw() +
  xlim(c(2005,2024))

# Adjacency matrix for pairwise synchrony comparisons
poc_mra_coords<-mra_sites %>%
  dplyr::mutate(Habitat = ifelse(Habitat == "BR", "Backreef", Habitat)) %>%
  dplyr::mutate(site_hab = paste0(Site,"_",Habitat)) %>%
  dplyr::filter(site_hab %in% row.names(poc_df_clean$cdat)) %>%
  dplyr::arrange(match(site_hab, row.names(poc_df_clean$cdat))) %>%
  dplyr::select(Lat,Long) %>%
  dplyr::rename("X" = "Long", "Y" = "Lat")


poc_wmf_plotdat_df<-extract_plotmag_wmf_data(poc_fres_wmf) # extract using helper function
poc_wmf_pretty_TS <- pretty(poc_wmf_plotdat_df$Timescale, n = 8) # "pretty" time axis labels


poc_whole_wmf_plot<-poc_wmf_plotdat_df %>%
  na.omit() %>%
  ggplot(aes(x = Time, y = log2(Timescale), fill = Magnitude)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Coherence") +
  scale_y_continuous(
    breaks = log2(poc_wmf_pretty_TS),
    labels = round(poc_wmf_pretty_TS, 2)) +
  scale_x_continuous(
    breaks = seq(floor(min(poc_wmf_plotdat_df$Time)),
                 ceiling(max(poc_wmf_plotdat_df$Time)),
                 by = 1),
    labels = seq(floor(min(poc_wmf_plotdat_df$Time)),
                 ceiling(max(poc_wmf_plotdat_df$Time)),
                 by = 1)) +
  theme_bw() +
  labs(x = "Time",y = "Time period (years)") +
  global_fill_scale


### Extracting contours ####
poc_obs_vals <- Mod(poc_fres_wmf$values)  # WMF magnitude
poc_mat= as.matrix(poc_df_wide)

# Number of permutations
nrand <- 1000

# Get matrix dimensions
n_sites <- nrow(poc_mat)
n_time <- ncol(poc_mat)

# For storing null WMFs
poc_null_vals <- array(NA, dim = c(nrow(poc_obs_vals), ncol(poc_obs_vals), nrand))
poc_cleaned <- cleandat(poc_mat, times = years, clev = bio_clev)
poc_cdat<-poc_cleaned$cdat


# Loop through permutations
for (i in 1:nrand) {
  # Ciruclar rotation of time series - preserve spectral structure within sites
  surrogates <- t(apply(poc_cdat, 1, generate_ar1_surrogate))

  # Re-center - site's time series to zero mean
  surrogate_mat <- sweep(surrogates, 1, rowMeans(surrogates, na.rm = TRUE), "-")
  colMeans(surrogate_mat)  # should all be ~0

  # WMF
  shuffled_wmf <- wmf(surrogate_mat, times = years, f0 = 0.5)

  # magnitude
  poc_null_vals[,,i] <- Mod(shuffled_wmf$values)
}

#  empirical p-values
poc_pvals <- matrix(NA, nrow = nrow(poc_obs_vals), ncol = ncol(poc_obs_vals))
for (r in 1:nrow(poc_obs_vals)) {
  for (c in 1:ncol(poc_obs_vals)) {
    poc_pvals[r,c] <- mean(poc_null_vals[r,c,] >= poc_obs_vals[r,c])
  }
}

#  significance mask
poc_sig_mask <- poc_pvals < 0.05

wmf_df <- expand.grid(
  Time = shuffled_wmf$times,
  Timescale = shuffled_wmf$timescales
)
Timescale_pretty <- pretty(wmf_df$Timescale, n = 8)

poc_wmf_df<-wmf_df
poc_wmf_df$Magnitude <- as.vector(poc_obs_vals)
poc_wmf_df$Significant <- as.vector(poc_sig_mask)
poc_wmf_contours_df<-poc_wmf_df


### Wavelet plotting ####

# Short time series
poc_wmf_short_plotdat_df<-extract_plotmag_wmf_data(poc_fres_short_wmf)
poc_wmf_short_pretty_TS <- pretty(poc_wmf_plotdat_df$Timescale, n = 8)

# Full time series
poc_wmf_plotdat_df<-extract_plotmag_wmf_data(poc_fres_wmf)
poc_wmf_pretty_TS <- pretty(poc_wmf_plotdat_df$Timescale, n = 8)

poc_whole_wmf_plot<-poc_wmf_plotdat_df %>%
  na.omit() %>%
  ggplot(aes(x = Time, y = log2(Timescale), fill = Magnitude)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Coherence") +
  scale_y_continuous(
    breaks = log2(poc_wmf_pretty_TS),
    labels = round(poc_wmf_pretty_TS, 2)) +
  scale_x_continuous(
    breaks = seq(floor(min(poc_wmf_plotdat_df$Time)),
                 ceiling(max(poc_wmf_plotdat_df$Time)),
                 by = 1),
    labels = seq(floor(min(poc_wmf_plotdat_df$Time)),
                 ceiling(max(poc_wmf_plotdat_df$Time)),
                 by = 1)) +
  theme_bw() +
  labs(
    x = "Time",
    y = "Time period (years)") +
  global_fill_scale


poc_wmf_short_plotdat_df<-extract_plotmag_wmf_data(poc_fres_short_wmf)
poc_wmf_short_pretty_TS <- pretty(poc_wmf_plotdat_df$Timescale, n = 8)

poc_whole_wmf_short_plot<-poc_wmf_short_plotdat_df %>%
  na.omit() %>%
  ggplot(aes(x = Time, y = log2(Timescale), fill = Magnitude)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", name = "Coherence") +
  scale_y_continuous(
    breaks = log2(poc_wmf_short_pretty_TS),
    labels = round(poc_wmf_short_pretty_TS, 2)) +
  scale_x_continuous(
    breaks = seq(floor(min(poc_wmf_short_plotdat_df$Time)),
                 ceiling(max(poc_wmf_short_plotdat_df$Time)),
                 by = 1),
    labels = seq(floor(min(poc_wmf_short_plotdat_df$Time)),
                 ceiling(max(poc_wmf_short_plotdat_df$Time)),
                 by = 1)) +
  theme_bw() +
  labs(
    x = "Time",
    y = "Time period (years)") +
  global_fill_scale


# Full observed wavelet
poc_whole_wmf_plot_ready_marked<-poc_whole_wmf_plot +
  geom_vline(xintercept = 2010.5,linetype = "dashed", color = "black") +
  geom_vline(xintercept = 2018.5,linetype = "dashed", color = "black") +
  annotate(geom = "text", label = "Cyclone", color = "black", x = 2012.5, y = Inf, hjust = 0.8, vjust = 2.5) +
  annotate(geom = "text", label = "Bleaching", color = "black", x = 2021, y = Inf, hjust = 0.8,vjust = 2.5)


# Short observed wavelet
poc_whole_wmf_plot_short_marked<-poc_whole_wmf_short_plot +
  geom_vline(xintercept = 2010.5,linetype = "dashed", color = "black") +
  geom_vline(xintercept = 2018.5,linetype = "dashed", color = "black") +
  annotate(geom = "text", label = "Cyclone", color = "black", x = 2012.5, y = Inf, hjust = 0.8, vjust = 2.5) +
  annotate(geom = "text", label = "Bleaching", color = "black", x = 2021, y = Inf, hjust = 0.8,vjust = 2.5)




## Porites ####
por_df<-MCR_coral_data %>%
  dplyr::select(1:3,Porites_new) %>%
  dplyr::mutate(Site = str_replace(Site," ","0")) %>%
  dplyr::mutate(Site_hab = paste0(Site,"_",Habitat)) %>%
  dplyr::select(-c(Site,Habitat))


por_df_filled<-por_df %>%
  filter(Date %in% years) %>%
  dplyr::group_by(Site_hab) %>%
  tidyr::complete(Date = startyear:endyear, fill = list(Pocillopora = NA)) %>%
  arrange(Site_hab, Date) %>%
  dplyr::mutate(Porites_new = case_when(
    is.na(Porites_new) & Date != 2011 ~ 0,
    TRUE ~ Porites_new
  )) %>%
  dplyr::group_by(Site_hab) %>%
  dplyr::mutate(Porites_new = if_else(
    Date == 2011 & is.na(Porites_new),
    na.approx(Porites_new, maxgap = 2, rule = 2),
    Porites_new
  ))

por_df_wide<-por_df_filled %>%
  pivot_wider(values_from = Porites_new,names_from = Date) %>%
  as.data.frame()

row.names(por_df_wide)<-por_df_wide$Site_hab
por_df_wide$Site_hab<-NULL
por_mat<-as.matrix(por_df_wide)
row.names(por_mat)<-row.names(por_df_wide)

### Wavelet mean fields ####
# Full time series
por_df_clean<-cleandat(por_mat, times = years, clev = bio_clev) # removes mean, detrends, and standardizes variance
por_cdat<-as.vector(por_df_clean$cdat)
por_fres_wmf<-wmf(por_df_clean$cdat,years, f0 = 0.5)

row.names(por_df_clean$cdat)<-row.names(poc_df_wide)

# Cleandat and WMF for data with shorter temporal duration for WLM
por_mat_short <- por_mat[, colnames(por_mat) %in% env_years]
por_df_short_clean<-cleandat(por_mat_short, times = env_years, clev = bio_clev) # removes mean, detrends, and standardizes variance
por_fres_short_wmf<-wmf(por_df_short_clean$cdat,env_years, f0 = 0.5)

# Detrended and ready for time series plotting
por_detrended<-as.data.frame(por_df_clean$cdat) %>%
  dplyr::mutate(Site_hab = rownames(.)) %>%
  pivot_longer(-Site_hab, names_to = "Year",values_to = "cover") %>%
  separate(Site_hab, into = c("Site","Habitat"),sep = "_") %>%
  dplyr::mutate(Year = as.numeric(Year))

# Plotting means and var for Figure 3 - modifed
por_ts_plot_df<-por_detrended %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(mean = mean(cover,na.rm = T),
                   var = var(cover, na.rm = T),
                   sd = sd(cover, na.rm = T),
                   n = n(),
                   se = sd/sqrt(n))

por_ts_plot<-por_ts_plot_df %>%
  ggplot(aes(x = Year, y = mean)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(ymin = mean - se, ymax  = mean + se)) +
  geom_line() +
  theme_bw() +
  xlim(c(2005,2024)) +
  ylim(c(-4,4))



# Adjacency matrix for pairwise synchrony comparisons
por_mra_coords<-mra_sites %>%
  dplyr::mutate(Habitat = ifelse(Habitat == "BR", "Backreef", Habitat)) %>%
  dplyr::mutate(site_hab = paste0(Site,"_",Habitat)) %>%
  dplyr::filter(site_hab %in% row.names(por_df_clean$cdat)) %>%
  arrange(match(site_hab, row.names(por_df_clean$cdat))) %>%
  mutate(row.names(.) == site_hab) %>%
  dplyr::select(Lat,Long) %>%
  dplyr::rename("X" = "Long", "Y" = "Lat")


### Wavelet plotting ####
por_wmf_plotdat_df <- extract_plotmag_wmf_data(por_fres_wmf)
por_wmf_pretty_TS <- pretty(por_wmf_plotdat_df$Timescale, n = 8)

# Full time series
por_whole_wmf_plot<-por_wmf_plotdat_df %>%
  na.omit() %>%
  ggplot(aes(x = Time, y = log2(Timescale), fill = Magnitude)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Coherence") +
  scale_y_continuous(
    breaks = log2(por_wmf_pretty_TS),
    labels = round(por_wmf_pretty_TS, 2),
  ) +
  scale_x_continuous(
    breaks = seq(floor(min(por_wmf_plotdat_df$Time)),
                 ceiling(max(por_wmf_plotdat_df$Time)),
                 by = 1),
    labels = seq(floor(min(por_wmf_plotdat_df$Time)),
                 ceiling(max(por_wmf_plotdat_df$Time)),
                 by = 1)
  ) +
  theme_bw() +
  labs(
    x = "Time",
    y = "Time period (years)") +
  global_fill_scale


### Extracting contours ####
por_obs_vals <- Mod(por_fres_wmf$values)  # WMF magnitude
por_mat= as.matrix(por_df_wide)

# Number of permutations
nrand <- 1000

# Get matrix dimensions
n_sites <- nrow(por_mat)
n_time <- ncol(por_mat)

# For storing null WMFs
por_null_vals <- array(NA, dim = c(nrow(por_obs_vals), ncol(por_obs_vals), nrand))
por_cleaned <- cleandat(por_mat, times = years, clev = bio_clev)
por_cdat<-por_cleaned$cdat


# Loop through permutations
for (i in 1:nrand) {
  # Ciruclar rotation of time series - preserve spectral structure within sites
  surrogates <- t(apply(por_cdat, 1, generate_ar1_surrogate))

  # Re-center -- site's time series (i.e., columns) to zero mean
  surrogate_mat <- sweep(surrogates, 1, rowMeans(surrogates, na.rm = TRUE), "-")
  colMeans(surrogate_mat)  # should all be ~0

  # WMF
  shuffled_wmf <- wmf(surrogate_mat, times = years, f0 = 0.5)

  # magnitude
  por_null_vals[,,i] <- Mod(shuffled_wmf$values)
}

# Calculate empirical p-values
por_pvals <- matrix(NA, nrow = nrow(por_obs_vals), ncol = ncol(por_obs_vals))
for (r in 1:nrow(por_obs_vals)) {
  for (c in 1:ncol(por_obs_vals)) {
    por_pvals[r,c] <- mean(por_null_vals[r,c,] >= por_obs_vals[r,c])
  }
}

# Create significance mask
por_sig_mask <- por_pvals < 0.05

wmf_df <- expand.grid(
  Time = shuffled_wmf$times,
  Timescale = shuffled_wmf$timescales
)
Timescale_pretty <- pretty(wmf_df$Timescale, n = 8)

por_wmf_df<-wmf_df
por_wmf_df$Magnitude <- as.vector(por_obs_vals)
por_wmf_df$Significant <- as.vector(por_sig_mask)
por_wmf_contours_df<-por_wmf_df


# Short time series
por_wmf_short_plotdat_df<-extract_plotmag_wmf_data(por_fres_short_wmf)
por_wmf_short_pretty_TS <- pretty(por_wmf_plotdat_df$Timescale, n = 8)

por_whole_wmf_short_plot<-por_wmf_short_plotdat_df %>%
  na.omit() %>%
  ggplot(aes(x = Time, y = log2(Timescale), fill = Magnitude)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Coherence") +
  scale_y_continuous(
    breaks = log2(por_wmf_short_pretty_TS),
    labels = round(por_wmf_short_pretty_TS, 2)) +
  scale_x_continuous(
    breaks = seq(floor(min(por_wmf_short_plotdat_df$Time)),
                 ceiling(max(por_wmf_short_plotdat_df$Time)),
                 by = 1),  # Whole number year breaks
    labels = seq(floor(min(por_wmf_short_plotdat_df$Time)),
                 ceiling(max(por_wmf_short_plotdat_df$Time)),
                 by = 1)) +
  theme_bw() +
  labs(
    x = "Time",
    y = "Time period (years)") +
  global_fill_scale

por_whole_wmf_plot_ready <- por_whole_wmf_plot + global_fill_scale +
  theme(axis.title.y = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank())

por_whole_wmf_short_plot <- por_whole_wmf_short_plot + global_fill_scale +
  theme(axis.title.y = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank())

por_whole_wmf_short_plot_marked<-por_whole_wmf_short_plot +
  geom_vline(xintercept = 2010.5,linetype = "dashed", color = "black") +
  geom_vline(xintercept = 2018.5,linetype = "dashed", color = "black") +
  annotate(geom = "text", label = "Cyclone", color = "black", x = 2012.5, y = Inf, hjust = 0.8, vjust = 2.5) +
  annotate(geom = "text", label = "Bleaching", color = "black", x = 2021, y = Inf, hjust = 0.8,vjust = 2.5)

por_whole_wmf_plot_ready_marked<-por_whole_wmf_plot_ready +
  geom_vline(xintercept = 2010.5,linetype = "dashed", color = "black") +
  geom_vline(xintercept = 2018.5,linetype = "dashed", color = "black") +
  annotate(geom = "text", label = "Cyclone", color = "black", x = 2012.5, y = Inf, hjust = 0.8, vjust = 2.5) +
  annotate(geom = "text", label = "Bleaching", color = "black", x = 2021, y = Inf, hjust = 0.8,vjust = 2.5)






## Montipora ####
Mont_df<-MCR_coral_data %>%
  dplyr::select(1:3,Montipora) %>%
  dplyr::mutate(Site = str_replace(Site," ","0")) %>%
  dplyr::mutate(Site_hab = paste0(Site,"_",Habitat)) %>%
  dplyr::select(-c(Site,Habitat))


Mont_df_filled<-Mont_df %>%
  filter(Date %in% years) %>%
  dplyr::group_by(Site_hab) %>%
  tidyr::complete(Date = startyear:endyear, fill = list(Montipora = NA)) %>%
  arrange(Site_hab, Date) %>%
  mutate(Montipora = case_when(
    is.na(Montipora) & Date != 2011 ~ 0,
    TRUE ~ Montipora
  )) %>%
  dplyr::group_by(Site_hab) %>%
  dplyr::mutate(Montipora = if_else(
    Date == 2011 & is.na(Montipora),
    na.approx(Montipora, maxgap = 2, rule = 2),
    Montipora
  ))


Mont_df_wide<-Mont_df_filled %>%
  pivot_wider(values_from = Montipora,names_from = Date) %>%
  as.data.frame()

row.names(Mont_df_wide)<-Mont_df_wide$Site_hab
Mont_df_wide$Site_hab<-NULL
mont_mat<-as.matrix(Mont_df_wide)
row.names(mont_mat)<-row.names(Mont_df_wide)

## Wavelet mean fields ####
# Full time series
mont_df_clean<-cleandat(mont_mat, times = years, clev = bio_clev) # removes mean, detrends, and standardizes variance
mont_cdat<-as.vector(mont_df_clean$cdat)
row.names(mont_df_clean$cdat)<-row.names(Mont_df_wide)
mont_fres_wmf<-wmf(mont_df_clean$cdat,years,f0 = 0.5)

# Cleandat and WMF for data with shorter temporal duration for WLM
mont_mat_short <- mont_mat[, colnames(mont_mat) %in% env_years]
mont_df_short_clean<-cleandat(mont_mat_short, times = env_years, clev = bio_clev) # removes mean, detrends, and standardizes variance
mont_fres_short_wmf<-wmf(mont_df_short_clean$cdat,env_years,f0 = 0.5)


mont_detrended<-as.data.frame(mont_df_clean$cdat) %>%
  dplyr::mutate(Site_hab = rownames(.)) %>%
  pivot_longer(-Site_hab, names_to = "Year",values_to = "cover") %>%
  separate(Site_hab, into = c("Site","Habitat"),sep = "_") %>%
  dplyr::mutate(Year = as.numeric(Year))

# Plotting means and var for Figure 3 - modifed
Mont_ts_plot_df<-mont_detrended %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(mean = mean(cover,na.rm = T),
                   var = var(cover, na.rm = T),
                   sd = sd(cover, na.rm = T),
                   n = n(),
                   se = sd/sqrt(n))

mont_ts_plot<-Mont_ts_plot_df %>%
  ggplot(aes(x = Year , y = mean)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(ymin = mean - se, ymax  = mean + se)) +
  geom_line() +
  theme_bw() +
  xlim(c(2005,2024))


mont_mra_coords<-mra_sites %>%
  dplyr::mutate(Habitat = ifelse(Habitat == "BR", "Backreef", Habitat)) %>%
  dplyr::mutate(site_hab = paste0(Site,"_",Habitat)) %>%
  dplyr::filter(site_hab %in% row.names(mont_df_clean$cdat)) %>%
  arrange(match(site_hab, row.names(mont_df_clean$cdat))) %>%
  dplyr::select(Lat,Long) %>%
  dplyr::rename("X" = "Long", "Y" = "Lat")



## Wavelet plotting ####
mont_wmf_plotdat_df <- extract_plotmag_wmf_data(mont_fres_wmf)
mont_wmf_pretty_TS <- pretty(mont_wmf_plotdat_df$Timescale, n = 8)

mont_whole_wmf_plot<-mont_wmf_plotdat_df %>%
  na.omit() %>%
  ggplot(aes(x = Time, y = log2(Timescale), fill = Magnitude)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Coherence") +
  scale_y_continuous(
    breaks = log2(mont_wmf_pretty_TS),
    labels = round(mont_wmf_pretty_TS, 2)) +
  scale_x_continuous(
    breaks = seq(floor(min(mont_wmf_plotdat_df$Time)),
                 ceiling(max(mont_wmf_plotdat_df$Time)),
                 by = 1),
    labels = seq(floor(min(mont_wmf_plotdat_df$Time)),
                 ceiling(max(mont_wmf_plotdat_df$Time)),
                 by = 1)) +
  theme_bw() +
  labs(x = "Time", y = "Time period (years)") +
  global_fill_scale

## Extracting contours ####
mont_obs_vals <- Mod(mont_fres_wmf$values)  # WMF magnitude
mont_mat = as.matrix(Mont_df_wide)

# Number of permutations
nrand <- 1000

# Get matrix dimensions
n_sites <- nrow(mont_mat)
n_time <- ncol(mont_mat)

# For storing null WMFs
mont_null_vals <- array(NA, dim = c(nrow(mont_obs_vals), ncol(mont_obs_vals), nrand))
mont_cleaned <- cleandat(mont_mat, times = years, clev = bio_clev)
mont_cdat<-mont_cleaned$cdat


# Loop through permutations
for (i in 1:nrand) {
  # Ciruclar rotation of time series - preserve spectral structure within sites
  surrogates <- t(apply(mont_cdat, 1, generate_ar1_surrogate))

  # Re-center -- site's time series (i.e., columns) to zero mean
  surrogate_mat <- sweep(surrogates, 1, rowMeans(surrogates, na.rm = TRUE), "-")
  colMeans(surrogate_mat)  # should all be ~0

  # WMF
  shuffled_wmf <- wmf(surrogate_mat, times = years, f0 = 0.5)

  # magnitude
  mont_null_vals[,,i] <- Mod(shuffled_wmf$values)
}

# Calculate empirical p-values
mont_pvals <- matrix(NA, nrow = nrow(mont_obs_vals), ncol = ncol(mont_obs_vals))
for (r in 1:nrow(mont_obs_vals)) {
  for (c in 1:ncol(mont_obs_vals)) {
    mont_pvals[r,c] <- mean(mont_null_vals[r,c,] >= mont_obs_vals[r,c])
  }
}

# Create significance mask
mont_sig_mask <- mont_pvals < 0.05

wmf_df <- expand.grid(
  Time = shuffled_wmf$times,
  Timescale = shuffled_wmf$timescales
)
Timescale_pretty <- pretty(wmf_df$Timescale, n = 8)

mont_wmf_df<-wmf_df
mont_wmf_df$Magnitude <- as.vector(mont_obs_vals)
mont_wmf_df$Significant <- as.vector(mont_sig_mask)
mont_wmf_contours_df<-mont_wmf_df


mont_wmf_short_plotdat_df<-extract_plotmag_wmf_data(mont_fres_short_wmf)
mont_wmf_short_pretty_TS <- pretty(mont_wmf_plotdat_df$Timescale, n = 8)

mont_whole_wmf_short_plot<-mont_wmf_short_plotdat_df %>%
  na.omit() %>%
  ggplot(aes(x = Time, y = log2(Timescale), fill = Magnitude)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Coherence") +
  scale_y_continuous(
    breaks = log2(mont_wmf_short_pretty_TS),
    labels = round(mont_wmf_short_pretty_TS, 2)) +
  scale_x_continuous(
    breaks = seq(floor(min(mont_wmf_short_plotdat_df$Time)),
                 ceiling(max(mont_wmf_short_plotdat_df$Time)),
                 by = 1),  # Whole number year breaks
    labels = seq(floor(min(mont_wmf_short_plotdat_df$Time)),
                 ceiling(max(mont_wmf_short_plotdat_df$Time)),
                 by = 1)) +
  theme_bw() +
  labs(
    x = "Time",
    y = "Time period (years)") +
  global_fill_scale


mont_whole_wmf_plot_ready <- mont_whole_wmf_plot + global_fill_scale +
  theme(axis.title.y = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank())

mont_whole_wmf_short_plot <- mont_whole_wmf_short_plot + global_fill_scale +
  theme(axis.title.y = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank())



mont_whole_wmf_plot_ready_marked<-mont_whole_wmf_plot_ready +
  geom_vline(xintercept = 2010.5,linetype = "dashed", color = "black") +
  geom_vline(xintercept = 2018.5,linetype = "dashed", color = "black") +
  annotate(geom = "text", label = "Cyclone", color = "black", x = 2012.5, y = Inf, hjust = 0.8, vjust = 2.5) +
  annotate(geom = "text", label = "Bleaching", color = "black", x = 2021, y = Inf, hjust = 0.8,vjust = 2.5)

mont_whole_wmf_short_plot_marked<-mont_whole_wmf_short_plot +
  geom_vline(xintercept = 2010.5,linetype = "dashed", color = "black") +
  geom_vline(xintercept = 2018.5,linetype = "dashed", color = "black") +
  annotate(geom = "text", label = "Cyclone", color = "black", x = 2012.5, y = Inf, hjust = 0.8, vjust = 2.5) +
  annotate(geom = "text", label = "Bleaching", color = "black", x = 2021, y = Inf, hjust = 0.8,vjust = 2.5)






## Acropora ####
Acro_df<-MCR_coral_data %>%
  dplyr::select(1:3,Acropora) %>%
  dplyr::mutate(Site = str_replace(Site," ","0")) %>%
  dplyr::mutate(Site_hab = paste0(Site,"_",Habitat)) %>%
  dplyr::select(-c(Site,Habitat))


Acro_df_filled<-Acro_df %>%
  filter(Date %in% years) %>%
  dplyr::group_by(Site_hab) %>%
  tidyr::complete(Date = startyear:endyear, fill = list(Acropora = NA)) %>%
  arrange(Site_hab, Date) %>%
  dplyr::mutate(Acropora = case_when(
    is.na(Acropora) & Date != 2011 ~ 0,
    TRUE ~ Acropora
  )) %>%
  group_by(Site_hab) %>%
  dplyr::mutate(Acropora = if_else(
    Date == 2011 & is.na(Acropora),
    na.approx(Acropora, maxgap = 2, rule = 2),
    Acropora
  ))


Acro_df_wide<-Acro_df_filled %>%
  pivot_wider(values_from = Acropora,names_from = Date) %>%
  as.data.frame()

row.names(Acro_df_wide)<-Acro_df_wide$Site_hab
Acro_df_wide$Site_hab<-NULL
acro_mat<-as.matrix(Acro_df_wide)
row.names(acro_mat)<-row.names(Acro_df_wide)

acro_zero_rows <- rowSums(acro_mat) == 0
acro_mat_nozero <- acro_mat[!acro_zero_rows, ]

### Wavelet mean fields ####
# Full time series
acro_df_clean<-cleandat(acro_mat_nozero, times = years, clev = bio_clev) # removes mean, detrends, and standardizes variance
acro_cdat<-as.vector(acro_df_clean$cdat)
acro_fres_wmf<-wmf(acro_df_clean$cdat,years, f0 = 0.5)

row.names(acro_df_clean$cdat)<-row.names(acro_mat_nozero)

# Cleandat and WMF for data with shorter temporal duration for WLM
acro_mat_short <- acro_mat[, colnames(acro_mat) %in% env_years]
acro_short_zero_rows <- rowSums(acro_mat_short) == 0
acro_mat_short_nozero <- acro_mat_short[!acro_short_zero_rows, ]

acro_df_short_clean<-cleandat(acro_mat_short_nozero, times = env_years, clev = bio_clev) # removes mean, detrends, and standardizes variance
acro_fres_short_wmf<-wmf(acro_df_short_clean$cdat,env_years, f0 = 0.5)
row.names(acro_df_short_clean$cdat)<-row.names(acro_mat_short_nozero)

# Detrended and ready for time series plotting
acro_detrended<-as.data.frame(acro_df_clean$cdat) %>%
  mutate(Site_hab = rownames(.)) %>%
  pivot_longer(-Site_hab, names_to = "Year",values_to = "cover") %>%
  separate(Site_hab, into = c("Site","Habitat"),sep = "_") %>%
  mutate(Year = as.numeric(Year))

# Plotting means and var for Figure 3 - modifed
acro_ts_plot_df<-acro_detrended %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(mean = mean(cover,na.rm = T),
                   var = var(cover, na.rm = T),
                   sd = sd(cover, na.rm = T),
                   n = n(),
                   se = sd/sqrt(n))

acro_ts_plot<-acro_ts_plot_df %>%
  ggplot(aes(x = Year , y = mean)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(ymin = mean - se, ymax  = mean + se)) +
  geom_line() +
  theme_bw() +
  xlim(c(2005,2024))



# Adjacency matrix for pairwise synchrony comparisons
acro_mra_coords<-mra_sites %>%
  dplyr::mutate(Habitat = ifelse(Habitat == "BR", "Backreef", Habitat)) %>%
  dplyr::mutate(site_hab = paste0(Site,"_",Habitat)) %>%
  dplyr::filter(site_hab %in% row.names(acro_df_clean$cdat)) %>%
  arrange(match(site_hab, row.names(acro_df_clean$cdat))) %>%
  dplyr::select(Lat,Long) %>%
  dplyr::rename("X" = "Long", "Y" = "Lat")


### Wavelet plotting ####
acro_wmf_plotdat_df<-extract_plotmag_wmf_data(acro_fres_wmf)
acro_wmf_pretty_TS <- pretty(acro_wmf_plotdat_df$Timescale, n = 8)

acro_whole_wmf_plot<-acro_wmf_plotdat_df %>%
  na.omit() %>%
  ggplot(aes(x = Time, y = log2(Timescale), fill = Magnitude)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Coherence") +
  scale_y_continuous(
    breaks = log2(acro_wmf_pretty_TS),
    labels = round(acro_wmf_pretty_TS, 2),
  ) +
  scale_x_continuous(
    breaks = seq(floor(min(acro_wmf_plotdat_df$Time)),
                 ceiling(max(acro_wmf_plotdat_df$Time)),
                 by = 1),
    labels = seq(floor(min(acro_wmf_plotdat_df$Time)),
                 ceiling(max(acro_wmf_plotdat_df$Time)),
                 by = 1)) +
  theme_bw() +
  labs(x = "Time",y = "Time period (years)") +
  global_fill_scale

## Extracting contours ####
acro_obs_vals <- Mod(acro_fres_wmf$values)  # WMF magnitude
acro_mat = as.matrix(Acro_df_wide)

# Number of permutations
nrand <- 1000

# Get matrix dimensions
n_sites <- nrow(acro_mat)
n_time <- ncol(acro_mat)

# For storing null WMFs
acro_null_vals <- array(NA, dim = c(nrow(acro_obs_vals), ncol(acro_obs_vals), nrand))
acro_cleaned <- cleandat(acro_mat, times = years, clev = bio_clev)
acro_cdat<-acro_cleaned$cdat


# Loop through permutations
for (i in 1:nrand) {
  # Ciruclar rotation of time series - preserve spectral structure within sites
  surrogates <- t(apply(acro_cdat, 1, generate_ar1_surrogate))

  # Re-center -- site's time series (i.e., columns) to zero mean
  surrogate_mat <- sweep(surrogates, 1, rowMeans(surrogates, na.rm = TRUE), "-")
  colMeans(surrogate_mat)  # should all be ~0

  # WMF
  shuffled_wmf <- wmf(surrogate_mat, times = years, f0 = 0.5)

  # magnitude
  acro_null_vals[,,i] <- Mod(shuffled_wmf$values)
}

# Calculate empirical p-values
acro_pvals <- matrix(NA, nrow = nrow(acro_obs_vals), ncol = ncol(acro_obs_vals))
for (r in 1:nrow(acro_obs_vals)) {
  for (c in 1:ncol(acro_obs_vals)) {
    acro_pvals[r,c] <- mean(acro_null_vals[r,c,] >= acro_obs_vals[r,c])
  }
}

# Create significance mask
acro_sig_mask <- acro_pvals < 0.05

wmf_df <- expand.grid(
  Time = shuffled_wmf$times,
  Timescale = shuffled_wmf$timescales
)
Timescale_pretty <- pretty(wmf_df$Timescale, n = 8)

acro_wmf_df<-wmf_df
acro_wmf_df$Magnitude <- as.vector(acro_obs_vals)
acro_wmf_df$Significant <- as.vector(acro_sig_mask)
acro_wmf_contours_df<-acro_wmf_df

acro_wmf_short_plotdat_df<-extract_plotmag_wmf_data(acro_fres_short_wmf)
acro_wmf_short_pretty_TS <- pretty(acro_wmf_plotdat_df$Timescale, n = 8)

acro_whole_wmf_short_plot<-acro_wmf_short_plotdat_df %>%
  na.omit() %>%
  ggplot(aes(x = Time, y = log2(Timescale), fill = Magnitude)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Coherence") +
  scale_y_continuous(
    breaks = log2(acro_wmf_short_pretty_TS),
    labels = round(acro_wmf_short_pretty_TS, 2)) +
  scale_x_continuous(
    breaks = seq(floor(min(acro_wmf_short_plotdat_df$Time)),
                 ceiling(max(acro_wmf_short_plotdat_df$Time)),
                 by = 1),  # Whole number year breaks
    labels = seq(floor(min(acro_wmf_short_plotdat_df$Time)),
                 ceiling(max(acro_wmf_short_plotdat_df$Time)),
                 by = 1)) +
  theme_bw() +
  labs(
    x = "Time",
    y = "Time period (years)") +
  global_fill_scale


acro_whole_wmf_plot_ready <- acro_whole_wmf_plot + global_fill_scale +
  theme(axis.title.y = element_blank())

acro_whole_wmf_short_plot <- acro_whole_wmf_short_plot + global_fill_scale +
  theme(axis.title.y = element_blank())


acro_whole_wmf_plot_ready_marked<-acro_whole_wmf_plot_ready +
  geom_vline(xintercept = 2010.5,linetype = "dashed", color = "black") +
  geom_vline(xintercept = 2018.5,linetype = "dashed", color = "black") +
  annotate(geom = "text", label = "Cyclone", color = "black", x = 2012.5, y = Inf, hjust = 0.8, vjust = 2.5) +
  annotate(geom = "text", label = "Bleaching", color = "black", x = 2021, y = Inf, hjust = 0.8,vjust = 2.5)

acro_whole_wmf_short_plot_marked<-acro_whole_wmf_short_plot +
  geom_vline(xintercept = 2010.5,linetype = "dashed", color = "black") +
  geom_vline(xintercept = 2018.5,linetype = "dashed", color = "black") +
  annotate(geom = "text", label = "Cyclone", color = "black", x = 2012.5, y = Inf, hjust = 0.8, vjust = 2.5) +
  annotate(geom = "text", label = "Bleaching", color = "black", x = 2021, y = Inf, hjust = 0.8,vjust = 2.5)





# Collective wavelet plotting ####
# Combine the three species
poc_whole_wmf_plot_ready_marked <- poc_whole_wmf_plot_ready_marked + theme(axis.title.x = element_blank(),axis.text.x = element_blank()) + labs(x= "Year")
poc_whole_wmf_plot_short_marked <- poc_whole_wmf_plot_short_marked + theme(axis.title.x = element_blank(),axis.text.x = element_blank()) + labs(x= "Year")



poc_whole_wmf_plot_short_marked_upd<-poc_whole_wmf_plot_short_marked  +
  geom_hline(yintercept = log2(2), color = "white", alpha = 0.9, size = 1) +
  geom_hline(yintercept = log2(5), color = "white", alpha = 0.9, size = 1) +
  geom_hline(yintercept = log2(10), color = "white",alpha = 0.9, size = 1) +
  annotate("text", y = log2(mean(c(2,5))), x = 2012, label = "2-5y",color = "white") +
  annotate("text", y = log2(mean(c(5,10))), x = 2014, label = "5-10y",color = "white") +
  annotate("segment", x = 2009.3, xend = 2009.3, y = log2(2), yend = log2(10), linetype = "solid", color = "black", size = 1) +
  annotate("text", y = log2(11), x = 2010.3, label = "2-10y",color = "black")


por_whole_wmf_short_plot_marked_upd<-por_whole_wmf_short_plot_marked  +
  geom_hline(yintercept = log2(2), color = "white", alpha = 0.9, size = 1) +
  geom_hline(yintercept = log2(5), color = "white", alpha = 0.9, size = 1) +
  geom_hline(yintercept = log2(10), color = "white",alpha = 0.9, size = 1) +
  annotate("text", y = log2(mean(c(2,5))), x = 2012, label = "2-5y",color = "white") +
  annotate("text", y = log2(mean(c(5,10))), x = 2014, label = "5-10y",color = "white") +
  annotate("segment", x = 2009.3, xend = 2009.3, y = log2(2), yend = log2(10), linetype = "solid", color = "black", size = 1)


mont_whole_wmf_short_plot_marked_upd<-mont_whole_wmf_short_plot_marked  +
  geom_hline(yintercept = log2(2), color = "white", alpha = 0.9, size = 1) +
  geom_hline(yintercept = log2(5), color = "white", alpha = 0.9, size = 1) +
  geom_hline(yintercept = log2(10), color = "white",alpha = 0.9, size = 1) +
  annotate("text", y = log2(mean(c(2,5))), x = 2012, label = "2-5y",color = "white") +
  annotate("text", y = log2(mean(c(5,10))), x = 2014, label = "5-10y",color = "white") +
  annotate("segment", x = 2009.3, xend = 2009.3, y = log2(2), yend = log2(10), linetype = "solid", color = "black", size = 1)


acro_whole_wmf_short_plot_marked_upd<-acro_whole_wmf_short_plot_marked  +
  geom_hline(yintercept = log2(2), color = "white", alpha = 0.9, size = 1) +
  geom_hline(yintercept = log2(5), color = "white", alpha = 0.9, size = 1) +
  geom_hline(yintercept = log2(10), color = "white",alpha = 0.9, size = 1) +
  annotate("text", y = log2(mean(c(2,5))), x = 2012, label = "2-5y",color = "white") +
  annotate("text", y = log2(mean(c(5,10))), x = 2014, label = "5-10y",color = "white") +
  annotate("segment", x = 2009.3, xend = 2009.3, y = log2(2), yend = log2(10), linetype = "solid", color = "black", size = 1)










## Figure 3 ####
poc_whole_wmf_plot_ready_marked_noflip<-poc_whole_wmf_plot +
  geom_contour(data = poc_wmf_contours_df, aes(z = as.numeric(Significant)), breaks = 0.5, color = "white", size = 0.3) +
  geom_vline(xintercept = 2010,linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = 2019,linetype = "dashed", color = "darkgrey") +
  annotate(geom = "text", label = "Oli", color = "black", x = 2009, y = Inf, vjust = 1.3) +
  annotate(geom = "text", label = "MHW", color = "black", x = 2021, y = Inf, vjust = 1.3) +
  scale_x_continuous(n.breaks = 16, limits = c(2005,2024)) +
  theme(axis.text.x = element_blank(),
        plot.title = element_markdown()) +
  global_fill_scale +
  labs(title = "(A) *Pocillopora* spp.") +
  geom_hline(yintercept = log2(2), color = "white", alpha = 0.9, size = 1) +
  geom_hline(yintercept = log2(5), color = "white", alpha = 0.9, size = 1) +
  geom_hline(yintercept = log2(10), color = "white",alpha = 0.9, size = 1) +
  annotate("text", y = log2(mean(c(2,5))), x = 2012, label = "2-5y",color = "white") +
  annotate("text", y = log2(mean(c(7,10))), x = 2014, label = "5-10y",color = "white") +
  annotate("segment", x = 2006, xend = 2006, y = log2(2), yend = log2(10), linetype = "solid", color = "black", size = 1) +
  annotate("text", y = log2(6), x = 2005.2, label = "2-10y",color = "black", angle = 90)


por_whole_wmf_plot_ready_marked_noflip<-por_whole_wmf_plot +
  geom_contour(data = por_wmf_contours_df, aes(z = as.numeric(Significant)), breaks = 0.5, color = "white", size = 0.3) +
  geom_vline(xintercept = 2010,linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = 2019,linetype = "dashed", color = "darkgrey") +
  scale_x_continuous(n.breaks = 16, limits = c(2005,2024)) +
  theme(axis.text.x = element_blank(),
        plot.title = element_markdown()) +
  global_fill_scale +
  labs(title = "(B) *Porites* spp.") +
  geom_hline(yintercept = log2(2), color = "white", alpha = 0.9, size = 1) +
  geom_hline(yintercept = log2(5), color = "white", alpha = 0.9, size = 1) +
  geom_hline(yintercept = log2(10), color = "white",alpha = 0.9, size = 1) +
  annotate("segment", x = 2006, xend = 2006, y = log2(2), yend = log2(10), linetype = "solid", color = "black", size = 1)

mont_whole_wmf_plot_ready_marked_noflip<-mont_whole_wmf_plot +
  geom_contour(data = mont_wmf_contours_df, aes(z = as.numeric(Significant)), breaks = 0.5, color = "white", size = 0.3) +
  geom_vline(xintercept = 2010,linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = 2019,linetype = "dashed", color = "darkgrey") +
  scale_x_continuous(n.breaks = 16, limits = c(2005,2024)) +
  theme(axis.text.x = element_blank(),
        plot.title = element_markdown()) +
  global_fill_scale +
  labs(title = "(C) *Montipora* spp.") +
  geom_hline(yintercept = log2(2), color = "white", alpha = 0.9, size = 1) +
  geom_hline(yintercept = log2(5), color = "white", alpha = 0.9, size = 1) +
  geom_hline(yintercept = log2(10), color = "white",alpha = 0.9, size = 1) +
  annotate("segment", x = 2006, xend = 2006, y = log2(2), yend = log2(10), linetype = "solid", color = "black", size = 1)



acro_whole_wmf_plot_ready_marked_noflip<-acro_whole_wmf_plot +
  geom_contour(data = acro_wmf_contours_df, aes(z = as.numeric(Significant)), breaks = 0.5, color = "white", size = 0.3) +
  geom_vline(xintercept = 2010,linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = 2019,linetype = "dashed", color = "darkgrey") +
  scale_x_continuous(n.breaks = 16, limits = c(2005,2024)) +
  theme(axis.text.x = element_blank(),
        plot.title = element_markdown()) +
  global_fill_scale +
  labs(title = "(D) *Acropora* spp.") +
  geom_hline(yintercept = log2(2), color = "white", alpha = 0.9, size = 1) +
  geom_hline(yintercept = log2(5), color = "white", alpha = 0.9, size = 1) +
  geom_hline(yintercept = log2(10), color = "white",alpha = 0.9, size = 1) +
  annotate("segment", x = 2006, xend = 2006, y = log2(2), yend = log2(10), linetype = "solid", color = "black", size = 1)



# Formatting plots together
Poc_1<-poc_whole_wmf_plot_ready_marked_noflip / poc_ts_plot + plot_layout(guides = "collect",axes = "collect",axis_titles = "collect", heights = c(4,1))
Por_2<-por_whole_wmf_plot_ready_marked_noflip / por_ts_plot + plot_layout(guides = "collect",axes = "collect",axis_titles = "collect", heights = c(4,1))
Mont_3<-mont_whole_wmf_plot_ready_marked_noflip / mont_ts_plot + plot_layout(guides = "collect",axes = "collect",axis_titles = "collect", heights = c(4,1))
Acro_4<-acro_whole_wmf_plot_ready_marked_noflip / acro_ts_plot + plot_layout(guides = "collect",axes = "collect",axis_titles = "collect", heights = c(4,1))


# Top panel with WMFs
fig_1_top<-(poc_whole_wmf_plot_ready_marked_noflip | por_whole_wmf_plot_ready_marked_noflip | mont_whole_wmf_plot_ready_marked_noflip | acro_whole_wmf_plot_ready_marked_noflip) + plot_layout(guides = "collect",axes = "collect", axis_titles = "collect") & theme(axis.title.x = element_blank()) & labs(y = "Timescale (years)")

poc_ts_plot_ready<-poc_ts_plot + theme(plot.margin = unit(c(0,0,0,0),"cm")) + labs(x = "Year") + geom_vline(xintercept = 2010, linetype = "dashed") + geom_vline(xintercept = 2019, linetype = "dashed") + annotate("text", x = -Inf, y = Inf, , hjust = -0.5, vjust = 1.5 , label = "(E)")
por_ts_plot_ready<-por_ts_plot + theme(plot.margin = unit(c(0,0,0,0),"cm")) + labs(x = "Year") + geom_vline(xintercept = 2010, linetype = "dashed") + geom_vline(xintercept = 2019, linetype = "dashed") + annotate("text", x = -Inf, y = Inf, , hjust = -0.5, vjust = 1.5 , label = "(F)")
mont_ts_plot_ready<-mont_ts_plot + theme(plot.margin = unit(c(0,0,0,0),"cm")) + labs(x = "Year") + geom_vline(xintercept = 2010, linetype = "dashed") + geom_vline(xintercept = 2019, linetype = "dashed") + annotate("text", x = -Inf, y = Inf, , hjust = -0.5, vjust = 1.5 , label = "(G)")
acro_ts_plot_ready<-acro_ts_plot + theme(plot.margin = unit(c(0,0,0,0),"cm")) + labs(x = "Year") + geom_vline(xintercept = 2010, linetype = "dashed") + geom_vline(xintercept = 2019, linetype = "dashed") + annotate("text", x = -Inf, y = Inf, , hjust = -0.5, vjust = 1.5 , label = "(H)")

# Bottom panel with cleaned time-series'
fig_1_bot<-(poc_ts_plot_ready | por_ts_plot_ready | mont_ts_plot_ready | acro_ts_plot_ready) + plot_layout(guides = "collect",axis_titles = "collect") & labs(y = "Detrended mean\n% cover ± se")


combined_fig3<-fig_1_top/fig_1_bot + plot_layout(guides = "collect",axes = "collect",axis_titles = "collect", heights = c(5,2))

ggsave('./figures/Figure_3_ADA.jpeg',
       combined_fig3,
       dpi = 300,
       width = 14,
       height = 4)


# Wavelet summary stats | Corals ####
poc_wmf_short_plotdat_df %>% dplyr::summarise(mean(Magnitude,na.rm = T))
por_wmf_short_plotdat_df %>% dplyr::summarise(mean(Magnitude,na.rm = T))
mont_wmf_short_plotdat_df %>% dplyr::summarise(mean(Magnitude,na.rm = T))
acro_wmf_short_plotdat_df %>% dplyr::summarise(mean(Magnitude,na.rm = T))




# Environmental data wavelets ####
#global_env_fill_scale <- scale_fill_viridis_c(limits = c(0, 1.75), option = "turbo", name = "Synchrony")  # old, bad for BG colorblind
global_env_fill_scale <- scale_fill_viridis_c(limits = c(0, 1.75), option = "magma", name = "Synchrony")

## Algae ####
alg_clean_cdat
alg_fres_wmf<-wmf(alg_clean_cdat,env_years, f0 =0.5)
plotmag(alg_fres_wmf)

# Extract
alg_wmf_plotdat_df<-extract_plotmag_wmf_data(alg_fres_wmf)
alg_wmf_pretty_TS <- pretty(alg_wmf_plotdat_df$Timescale, n = 8)

alg_whole_wmf_plot<-alg_wmf_plotdat_df %>%
  na.omit() %>%
  ggplot(aes(x = Time, y = log2(Timescale), fill = Magnitude)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Coherence") +

  scale_y_continuous(
    breaks = log2(alg_wmf_pretty_TS),
    labels = round(alg_wmf_pretty_TS, 2)) +
  scale_x_continuous(
    breaks = seq(floor(min(alg_wmf_plotdat_df$Time)),
                 ceiling(max(alg_wmf_plotdat_df$Time)),
                 by = 1),
    labels = seq(floor(min(alg_wmf_plotdat_df$Time)),
                 ceiling(max(alg_wmf_plotdat_df$Time)),
                 by = 1)) +
  theme_bw() +
  labs(x = "Time", y = "Time period (years)") +
  global_env_fill_scale +
  geom_vline(xintercept = 2010,linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = 2019,linetype = "dashed", color = "darkgrey") +
  annotate(geom = "text", label = "Oli", color = "black", x = 2009, y = Inf, vjust = 1.3) +
  annotate(geom = "text", label = "MHW", color = "black", x = 2021, y = Inf, vjust = 1.3) +
  scale_x_continuous(n.breaks = length(env_years), limits = c(min(env_years),max(env_years))) +
  theme(axis.text.x = element_blank()) +
  labs(title = "(A) Algal cover") +
  geom_hline(yintercept = log2(2), color = "white", alpha = 0.9, size = 1) +
  geom_hline(yintercept = log2(5), color = "white", alpha = 0.9, size = 1) +
  geom_hline(yintercept = log2(10), color = "white",alpha = 0.9, size = 1) +
  annotate("text", y = log2(mean(c(2,5))), x = 2012, label = "2-5y",color = "white") +
  annotate("text", y = log2(mean(c(5,10))), x = 2014, label = "5-10y",color = "white") +
  annotate("segment", x = 2009, xend = 2009, y = log2(2), yend = log2(10), linetype = "solid", color = "black", size = 1) +
  annotate("text", y = log2(6), x = 2008.2, label = "2-10y",color = "black", angle = 90)



# Underpanel
alg_detrended<-as.data.frame(alg_clean_cdat) %>%
  mutate(Site_hab = rownames(.)) %>%
  pivot_longer(-Site_hab, names_to = "Year",values_to = "cover") %>%
  separate(Site_hab, into = c("Site","Habitat"),sep = "_") %>%
  mutate(Year = as.numeric(Year))

alg_ts_plot_df<-alg_detrended %>%
  group_by(Year) %>%
  dplyr::summarise(mean = mean(cover,na.rm = T),
                   var = var(cover, na.rm = T),
                   sd = sd(cover, na.rm = T),
                   n = n(),
                   se = sd/sqrt(n))

alg_ts_plot<-alg_ts_plot_df %>%
  ggplot(aes(x = Year , y = mean)) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  geom_pointrange(aes(ymin = mean - se, ymax  = mean + se), color = "darkgreen") +
  geom_line(color = "darkgreen") +
  scale_x_continuous(n.breaks = length(env_years)/1.5, limits = c(min(env_years),max(env_years))) +
  scale_y_continuous(breaks = c(-1,0,1), limits = c(-1.5,1.5)) +
  theme_bw()


## DHD ####
DHD_clean_cdat
dhd_fres_wmf<-wmf(DHD_clean_cdat,env_years, f0 = 0.5)
plotmag(dhd_fres_wmf)

# Extract
dhd_wmf_plotdat_df<-extract_plotmag_wmf_data(dhd_fres_wmf)
dhd_wmf_pretty_TS <- pretty(dhd_wmf_plotdat_df$Timescale, n = 8)

dhd_whole_wmf_plot<-dhd_wmf_plotdat_df %>%
  na.omit() %>%
  ggplot(aes(x = Time, y = log2(Timescale), fill = Magnitude)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Coherence") +
  scale_y_continuous(
    breaks = log2(dhd_wmf_pretty_TS),
    labels = round(dhd_wmf_pretty_TS, 2),
  ) +
  scale_x_continuous(
    breaks = seq(floor(min(dhd_wmf_plotdat_df$Time)),
                 ceiling(max(dhd_wmf_plotdat_df$Time)),
                 by = 1),
    labels = seq(floor(min(dhd_wmf_plotdat_df$Time)),
                 ceiling(max(dhd_wmf_plotdat_df$Time)),
                 by = 1)
  ) +
  theme_bw() +
  labs(x = "Time", y = "Time period (years)") +
  global_env_fill_scale +
  geom_vline(xintercept = 2010,linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = 2019,linetype = "dashed", color = "darkgrey") +
  annotate(geom = "text", label = "Oli", color = "black", x = 2009, y = Inf, vjust = 1.3) +
  annotate(geom = "text", label = "MHW", color = "black", x = 2021, y = Inf, vjust = 1.3) +
  scale_x_continuous(n.breaks = length(env_years), limits = c(min(env_years),max(env_years))) +
  theme(axis.text.x = element_blank()) +
  labs(title = "(B) DHD") +
  geom_hline(yintercept = log2(2), color = "white", alpha = 0.9, size = 1) +
  geom_hline(yintercept = log2(5), color = "white", alpha = 0.9, size = 1) +
  geom_hline(yintercept = log2(10), color = "white",alpha = 0.9, size = 1) +
  annotate("segment", x = 2009, xend = 2009, y = log2(2), yend = log2(10), linetype = "solid", color = "black", size = 1)



# Underpanel
dhd_detrended<-as.data.frame(DHD_clean_cdat) %>%
  mutate(Site_hab = rownames(.)) %>%
  pivot_longer(-Site_hab, names_to = "Year",values_to = "cover") %>%
  separate(Site_hab, into = c("Site","Habitat"),sep = "_") %>%
  mutate(Year = as.numeric(Year))

dhd_ts_plot_df<-dhd_detrended %>%
  group_by(Year) %>%
  dplyr::summarise(mean = mean(cover,na.rm = T),
                   var = var(cover, na.rm = T),
                   sd = sd(cover, na.rm = T),
                   n = n(),
                   se = sd/sqrt(n))

dhd_ts_plot<-dhd_ts_plot_df %>%
  ggplot(aes(x = Year , y = mean)) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  geom_pointrange(aes(ymin = mean - se, ymax  = mean + se), color = "red") +
  geom_line(color = "red") +
  scale_x_continuous(n.breaks = length(env_years)/1.5, limits = c(min(env_years),max(env_years))) +
  scale_y_continuous(breaks = c(-1,0,1), limits = c(-1.7,1.5)) +
  theme_bw()


## DTR ####
dtr_clean_cdat
dtr_fres_wmf<-wmf(dtr_clean_cdat,env_years,f0 = 0.5)
plotmag(dtr_fres_wmf)

# Extract
dtr_wmf_plotdat_df<-extract_plotmag_wmf_data(dtr_fres_wmf)
dtr_wmf_pretty_TS <- pretty(dtr_wmf_plotdat_df$Timescale, n = 8)

dtr_whole_wmf_plot<-dtr_wmf_plotdat_df %>%
  na.omit() %>%
  ggplot(aes(x = Time, y = log2(Timescale), fill = Magnitude)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Coherence") +
  scale_y_continuous(
    breaks = log2(dtr_wmf_pretty_TS),
    labels = round(dtr_wmf_pretty_TS, 2)) +
  scale_x_continuous(
    breaks = seq(floor(min(dtr_wmf_plotdat_df$Time)),
                 ceiling(max(dtr_wmf_plotdat_df$Time)),
                 by = 1),
    labels = seq(floor(min(dtr_wmf_plotdat_df$Time)),
                 ceiling(max(dtr_wmf_plotdat_df$Time)),
                 by = 1)) +
  theme_bw() +
  labs(title = "DTR",x = "Time",y = "Time period (years)") +
  global_env_fill_scale +
  geom_vline(xintercept = 2010,linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = 2019,linetype = "dashed", color = "darkgrey") +
  annotate(geom = "text", label = "Oli", color = "black", x = 2009, y = Inf, vjust = 1.3) +
  annotate(geom = "text", label = "MHW", color = "black", x = 2021, y = Inf, vjust = 1.3) +
  scale_x_continuous(n.breaks = length(env_years), limits = c(min(env_years),max(env_years))) +
  theme(axis.text.x = element_blank()) +
  labs(title = "(C) DTR") +
  geom_hline(yintercept = log2(2), color = "white", alpha = 0.9, size = 1) +
  geom_hline(yintercept = log2(5), color = "white", alpha = 0.9, size = 1) +
  geom_hline(yintercept = log2(10), color = "white",alpha = 0.9, size = 1) +
  annotate("segment", x = 2009, xend = 2009, y = log2(2), yend = log2(10), linetype = "solid", color = "black", size = 1)




# Underpanel
dtr_detrended<-as.data.frame(dtr_clean_cdat) %>%
  mutate(Site_hab = rownames(.)) %>%
  pivot_longer(-Site_hab, names_to = "Year",values_to = "cover") %>%
  separate(Site_hab, into = c("Site","Habitat"),sep = "_") %>%
  mutate(Year = as.numeric(Year))

dtr_ts_plot_df<-dtr_detrended %>%
  group_by(Year) %>%
  dplyr::summarise(mean = mean(cover,na.rm = T),
                   var = var(cover, na.rm = T),
                   sd = sd(cover, na.rm = T),
                   n = n(),
                   se = sd/sqrt(n))

dtr_ts_plot<-dtr_ts_plot_df %>%
  ggplot(aes(x = Year , y = mean)) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  geom_pointrange(aes(ymin = mean - se, ymax  = mean + se), color = "blue") +
  geom_line(color = "blue") +
  scale_x_continuous(n.breaks = length(env_years)/1.5, limits = c(min(env_years),max(env_years))) +
  scale_y_continuous(breaks = c(-1,0,1), limits = c(-1.1,1.1)) +
  theme_bw()


## Plot predictor WMFs together ####
alg_whole_wmf_plot + dhd_whole_wmf_plot + dtr_whole_wmf_plot + plot_layout(guides = "collect", axes = "collect", axis_titles = "collect")

env_supp_top<-(alg_whole_wmf_plot | dhd_whole_wmf_plot | dtr_whole_wmf_plot) + plot_layout(guides = "collect",axes = "collect", axis_titles = "collect") & theme(axis.title.x = element_blank()) & labs(y = "Timescale (years)")

alg_ts_plot_ready<-alg_ts_plot +
  labs(x = "Time") + geom_vline(xintercept = 2010, linetype = "dashed") + geom_vline(xintercept = 2019, linetype = "dashed") + annotate("text", x = -Inf, y = Inf, , hjust = -0.5, vjust = 1.5 , label = "(D)")

dhd_ts_plot_ready<-dhd_ts_plot +
  labs(x = "Time") + geom_vline(xintercept = 2010, linetype = "dashed") + geom_vline(xintercept = 2019, linetype = "dashed") + annotate("text", x = -Inf, y = Inf, , hjust = -0.5, vjust = 1.5 , label = "(E)")

dtr_ts_plot_ready<-dtr_ts_plot +
  labs(x = "Time") + geom_vline(xintercept = 2010, linetype = "dashed") + geom_vline(xintercept = 2019, linetype = "dashed") + annotate("text", x = -Inf, y = Inf, , hjust = -0.5, vjust = 1.5 , label = "(F)")

env_supp_bot<-(alg_ts_plot_ready | dhd_ts_plot_ready | dtr_ts_plot_ready) + plot_layout(guides = "collect",axes = "collect",axis_titles = "collect") & labs(y = "cleaned trend\n in metric")

combined_env_supp<-env_supp_top/env_supp_bot + plot_layout(guides = "collect",axes = "collect",axis_titles = "collect", heights = c(5,2))

ggsave('./figures/Figure_S4_ADA.jpeg',
       combined_env_supp,
       dpi = 300,
       width = 15,
       height = 5)


## Environmental variable time series ####
DTR_df<-DTR_matrix %>%
  as.data.frame() %>%
  mutate(Site_hab = row.names(.)) %>%
  pivot_longer(-Site_hab, names_to = "Year",values_to = "DTR")

DHD_df<-DHD_matrix %>%
  as.data.frame() %>%
  mutate(Site_hab = row.names(.)) %>%
  pivot_longer(-Site_hab, names_to = "Year",values_to = "DHD")

Alg_df<-alg_matrix %>%
  as.data.frame() %>%
  mutate(Site_hab = row.names(.)) %>%
  pivot_longer(-Site_hab, names_to = "Year",values_to = "Algal cover")

env_pred_df<-DTR_df %>%
  merge(DHD_df) %>%
  merge(Alg_df) %>%
  separate(Site_hab, into = c("Site","Habitat"), sep = "_")


DHD_main_ts<-env_pred_df %>%
  ggplot(aes(x = as.numeric(Year), y = DHD, color = Habitat, group = Habitat)) +
  geom_rect(aes(xmin = 2006, xmax = 2010, ymin = -Inf, ymax = Inf), color = NA, fill = "grey", alpha = 0.1) +
  geom_point() +
  geom_path() +
  facet_grid( ~ Site) +
  theme_bw() +
  labs(x = "Year", y = "DHD") +
  scale_x_continuous(breaks=seq(2006, 2024, 2)) +
  theme(legend.position = "top",
        plot.margin = unit(c(0,0,0,0),"cm")) +
  scale_color_manual(values = c("#EC7D6D","#00BA38","#C67CFF","#6B9EF6")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



DTR_main_ts<-env_pred_df %>%
  ggplot(aes(x = as.numeric(Year), y = DTR, color = Habitat, group = Habitat)) +
  geom_rect(aes(xmin = 2006, xmax = 2010, ymin = -Inf, ymax = Inf), color = NA, fill = "grey", alpha = 0.1) +
  geom_point() +
  geom_path() +
  facet_grid( ~ Site) +
  theme_bw() +
  labs(x = "Year", y = "DTR") +
  scale_x_continuous(breaks=seq(2006, 2024, 2)) +
  theme(legend.position = "top",
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm")) +
  scale_color_manual(values = c("#EC7D6D","#00BA38","#C67CFF","#6B9EF6")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



Alg_main_ts<-env_pred_df %>%
  ggplot(aes(x = as.numeric(Year), y = `Algal cover`, color = Habitat, group = Habitat)) +
  geom_rect(aes(xmin = 2006, xmax = 2010, ymin = -Inf, ymax = Inf), color = NA, fill = "grey", alpha = 0.1) +
  geom_point() +
  geom_path() +
  facet_grid( ~ Site) +
  theme_bw() +
  labs(x = "Year", y = "Algal cover") +
  scale_x_continuous(breaks=seq(2006, 2024, 2)) +
  theme(legend.position = "top",
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm")) +
  scale_color_manual(values = c("#EC7D6D","#00BA38","#C67CFF","#6B9EF6")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



env_supp_ts<-(DHD_main_ts / DTR_main_ts / Alg_main_ts) + plot_layout(guides = "collect",axes = "collect",axis_titles = "collect") &
  theme(legend.position = "top") &
  geom_vline(xintercept = 2010, linetype = "dashed", alpha = 0.5) &
  geom_vline(xintercept = 2019, linetype = "dashed", alpha = 0.5)


ggsave("./figures/Figure_S3.png",
       env_supp_ts,
       dpi = 300,
       width = 11,
       height = 5)





# Timescale testing ####
short <- c(2,5)
medium <- c(5,10)
#long <-c(2,14) # OLD
long <-c(2,10)

bands <- list(
  short = c(2, 5),
  medium = c(5, 10),
  #long = c(2, 14) # old
  long = c(2,10)

)


## POC CLUSTERS ####
clustmethod = "ReXWT.sig.fft"
#clustmethod = "ReXWT.sig.aaft"
f0_set = 0.5

# short = 2-5y
poc_all_clust_short<-clust(poc_df_short_clean$cdat,env_years,poc_mra_coords,
                           method = clustmethod, nsurrogs = 1000,weighted = T,
                           tsrange = short,
                           f0 = f0_set)
# medium = 5-10y
poc_all_clust_medium<-clust(poc_df_short_clean$cdat,env_years,poc_mra_coords,
                            method = clustmethod, nsurrogs = 1000,weighted = T,
                            tsrange = medium,
                            f0 = f0_set)
# long = 2-10y --- changed 16OCT25 2-10y
poc_all_clust_long<-clust(poc_df_short_clean$cdat,env_years,poc_mra_coords,
                          method = clustmethod, nsurrogs = 1000,weighted = T,
                          tsrange = long,
                          f0 = f0_set)

## POR CLUSTERS ####
#clustmethod = "ReXWT.sig.fft"
por_all_clust_short<-clust(por_df_short_clean$cdat,env_years,por_mra_coords,
                           method = clustmethod, nsurrogs = 1000,weighted = T,
                           tsrange = short,
                           f0 = f0_set)

por_all_clust_medium<-clust(por_df_short_clean$cdat,env_years,poc_mra_coords,
                            method = clustmethod, nsurrogs = 1000,weighted = T,
                            tsrange = medium,
                            f0 = f0_set)

por_all_clust_long<-clust(por_df_short_clean$cdat,env_years,por_mra_coords,
                          method = clustmethod, nsurrogs = 1000,weighted = T,
                          tsrange = long,
                          f0 = f0_set)

## MONT CLUSTERS ####
#clustmethod = "ReXWT.sig.fft"
mont_all_clust_short<-clust(mont_df_short_clean$cdat,env_years,mont_mra_coords,
                            method = clustmethod, nsurrogs = 1000,weighted = T,
                            tsrange = short,
                            f0 = f0_set)

mont_all_clust_medium<-clust(mont_df_short_clean$cdat,env_years,mont_mra_coords,
                             method = clustmethod, nsurrogs = 1000,weighted = T,
                             tsrange = medium,
                             f0 = f0_set)

mont_all_clust_long<-clust(mont_df_short_clean$cdat,env_years,mont_mra_coords,
                           method = clustmethod, nsurrogs = 1000,weighted = T,
                           tsrange = long,
                           f0 = f0_set)


# ACRO CLUSTERS ####
#clustmethod = "ReXWT.sig.fft"
acro_mra_coords<-mra_sites %>%
  dplyr::mutate(Habitat = ifelse(Habitat == "BR", "Backreef", Habitat)) %>%
  dplyr::mutate(site_hab = paste0(Site,"_",Habitat)) %>%
  dplyr::filter(site_hab %in% row.names(acro_df_short_clean$cdat)) %>%
  arrange(match(site_hab, row.names(acro_df_short_clean$cdat))) %>%
  dplyr::select(Lat,Long) %>%
  dplyr::rename("X" = "Long", "Y" = "Lat")

acro_all_clust_short<-clust(acro_df_short_clean$cdat,env_years,acro_mra_coords,
                            method = clustmethod, nsurrogs = 1000,weighted = T,
                            tsrange = short,
                            f0 = f0_set)

acro_all_clust_medium<-clust(acro_df_short_clean$cdat,env_years,acro_mra_coords,
                             method = clustmethod, nsurrogs = 1000,weighted = T,
                             tsrange = medium,
                             f0 = f0_set)

acro_all_clust_long<-clust(acro_df_short_clean$cdat,env_years,acro_mra_coords,
                           method = clustmethod, nsurrogs = 1000,weighted = T,
                           tsrange = long,
                           f0 = f0_set)


# Environmental clusters ####
## DTR CLUSTERS ####
#clustmethod = "ReXWT.sig.fft"
dtr_mra_coords<-mra_sites %>%
  dplyr::mutate(Habitat = ifelse(Habitat == "BR", "Backreef", Habitat)) %>%
  dplyr::mutate(site_hab = paste0(Site,"_",Habitat)) %>%
  dplyr::filter(site_hab %in% row.names(dtr_clean$cdat)) %>%
  arrange(match(site_hab, row.names(dtr_clean$cdat))) %>%
  dplyr::select(Lat,Long) %>%
  dplyr::rename("X" = "Long", "Y" = "Lat")

dtr_clust_short<-clust(dtr_clean$cdat,env_years,dtr_mra_coords,
                       method = clustmethod, nsurrogs = 1000,weighted = T,
                       tsrange = short,
                       f0 = f0_set)

dtr_clust_medium<-clust(dtr_clean$cdat,env_years,dtr_mra_coords,
                        method = clustmethod, nsurrogs = 1000,weighted = T,
                        tsrange = medium,
                        f0 = f0_set)

dtr_clust_long<-clust(dtr_clean$cdat,env_years,dtr_mra_coords,
                      method = clustmethod, nsurrogs = 1000,weighted = T,
                      tsrange = long,
                      f0 = f0_set)



## DHD ####
#clustmethod = "ReXWT.sig.fft"
DHD_mra_coords<-mra_sites %>%
  dplyr::mutate(Habitat = ifelse(Habitat == "BR", "Backreef", Habitat)) %>%
  dplyr::mutate(site_hab = paste0(Site,"_",Habitat)) %>%
  dplyr::filter(site_hab %in% row.names(DHD_clean$cdat)) %>%
  arrange(match(site_hab, row.names(DHD_clean$cdat))) %>%
  dplyr::select(Lat,Long) %>%
  dplyr::rename("X" = "Long", "Y" = "Lat")

DHD_clust_short<-clust(DHD_clean$cdat,env_years,DHD_mra_coords,
                       method = clustmethod, nsurrogs = 1000,weighted = T,
                       tsrange = short,
                       f0 = f0_set)

DHD_clust_medium<-clust(DHD_clean$cdat,env_years,DHD_mra_coords,
                        method = clustmethod, nsurrogs = 1000,weighted = T,
                        tsrange = medium,
                        f0 = f0_set)

DHD_clust_long<-clust(DHD_clean$cdat,env_years,DHD_mra_coords,
                      method = clustmethod, nsurrogs = 1000,weighted = T,
                      tsrange = long,
                      f0 = f0_set)




## Algal cover ####
#clustmethod = "ReXWT.sig.fft"
alg_mra_coords<-mra_sites %>%
  dplyr::mutate(Habitat = ifelse(Habitat == "BR", "Backreef", Habitat)) %>%
  dplyr::mutate(site_hab = paste0(Site,"_",Habitat)) %>%
  dplyr::filter(site_hab %in% row.names(alg_clean$cdat)) %>%
  arrange(match(site_hab, row.names(alg_clean$cdat))) %>%
  dplyr::select(Lat,Long) %>%
  dplyr::rename("X" = "Long", "Y" = "Lat")

alg_clust_short<-clust(alg_clean$cdat,env_years,alg_mra_coords,
                       method = clustmethod, nsurrogs = 1000,weighted = T,
                       tsrange = short,
                       f0 = f0_set)

alg_clust_medium<-clust(alg_clean$cdat,env_years,alg_mra_coords,
                        method = clustmethod, nsurrogs = 1000,weighted = T,
                        tsrange = medium,
                        f0 = f0_set)

alg_clust_long<-clust(alg_clean$cdat,env_years,alg_mra_coords,
                      method = clustmethod, nsurrogs = 1000,weighted = T,
                      tsrange = long,
                      f0 = f0_set)



# END #

