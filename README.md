# Code in support of manuscript

Title: **Spatial portfolios in coral metapopulations are shaped by spatiotemporal asynchrony in environmental conditions**

Authors: **Griffin Srednick, Kristen Davis, Peter Edmunds**

In Consideration: Ecology Letters

**Folders:**

Scripts: contains scripts for:

- 'master_script.R' --> builds directories for downloaded data and runs all scripts for full analysis.

- 'therm_estimator_pub.R --> Curation of [MCR LTER Benthic Temperature data](https://doi.org/10.6073/pasta/02e0fa99c6fca29a1bdd26c46013f0f7). This interpolates missing values and provides tests for interpolation performance. Also estimates Degree Heating Days (DHD) and Diurnal Temperature Range (DTR) for use in '3_wavelet_analyses.R'.

- '1_coral_data_processing.R' --> Processing [coral cover data](https://portal.edirepository.org/nis/mapbrowse?scope=knb-lter-mcr&identifier=4) from the NSF Moorea Coral Reef Long-Term Ecological Research Program (MCR LTER)

- '2_algal_data_processing.R' --> Processing [macroalgal cover data](https://portal.edirepository.org/nis/mapbrowse?scope=knb-lter-mcr&identifier=8) from the NSF Moorea Coral Reef Long-Term Ecological Research Program (MCR LTER)

- '3_wavelet_analyses.R' --> Wavelet synchrony analyses of coral and predictor variables. Generates figures for main text and supplement.

- '4_synchrony_predictors.R' --> Assessing relationships between coral population synchrony and predictor variables. Generates Figure 4. 

- '5_portfolio_effects.R' --> Estimates portfolio effects. Assessing relationships between spatial synchrony and portfolio effects. Generates Figure 5.

- '6_summary_stats.R' --> Summary statistics from models. Generates Figure 2.



Data: Location of exported data for products of scripts above. Access raw data via the databases above.

Tables: output from modeling and summaries.

Figures: figure outputs.



**Note**: Some analyses (e.g., wavelet synchrony estimation using surrogates via wsyn::synmat()) employ surrogate generation. Although random seeds are set, minor variation in non-significant p-values may occur in analyses; all effect sizes, significance patterns, and conclusions are fully reproducible.
