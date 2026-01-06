# Script for running full analysis

# Make directories for placing relevant data ####

# Main data directory - place all coral and macroalgal data in this directory
ifelse(!dir.exists("./data"),
       dir.create("./data"), "Folder exists already")

# Directory for main temperature data -- place all temperature data (n = 12 csvs) except for LTER0 in this directory
ifelse(!dir.exists("./data/environmental/temperature/MAIN"),
       dir.create("./data/environmental/temperature/MAIN"), "Folder exists already")

# Directory for LTER0 temperature data
ifelse(!dir.exists("./data/environmental/temperature/LTER00"),
       dir.create("./data/environmental/temperature/LTER00"), "Folder exists already")

# Summarized temperature data -- output from script
ifelse(!dir.exists("./data/environmental/summarized"),
       dir.create("./data/environmental/summarized"), "Folder exists already")

# Exported results tables
ifelse(!dir.exists("./tables"),
       dir.create("./tables"), "Folder exists already")

# Exported figures
ifelse(!dir.exists("./figures"),
       dir.create("./figures"), "Folder exists already")


# Running analyses ####
# (1) Curating temperature data and estimating DHD & DTR
source("./scripts/therm_estimator_pub.R")

# (2) Curating and processing MCR coral data
source("./scripts/1_coral_data_processing.R")

# (3) Curating and processing MCR macroalgal (benthic cover) data
source("./scripts/2_algal_data_processing.R")

# (4) Running wavelet analyses with processed data
source("./scripts/3_wavelet_analyses.R")

# (5) Running analyses for predicting population synchrony
source("./scripts/4_synchrony_predictors.R")

# (6) Running analyses on portfolio effects
source("./scripts/5_portfolio_effects.R")

# (7) Running analyses on portfolio effects
source("./scripts/6_summary_stats.R")


# END ####
