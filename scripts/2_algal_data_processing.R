# Algae data processing - MCR LTER dataset: knb-lter-mcr.8 ####
# Curate and format for wavelet analyses

# Packages
library(tidyverse)

# Data load - knb-lter-mcr.8 - DOI: https://doi.org/10.6073/pasta/44d263de042710767ed770986fd12f5c
algal_data<-read.csv('./data/MCR_LTER_Annual_Survey_Benthic_Cover_20241219.csv')

# estimate algal cover
algal_data_wide <- algal_data %>%
  pivot_wider(names_from = Taxonomy_Substrate_Functional_Group, values_from = Percent_Cover, values_fill = 0) %>%
  filter(!`No data` == -1) %>%
  mutate(algal_cover = rowSums(across(!c(Year:Quadrat, Coral, Sand, Sponge, `Soft Coral`, `Coral Rubble`, `Bare Space`,`No data`), ~ ., .names = "selected_columns"), na.rm = TRUE)) %>%
  mutate(Habitat = ifelse(Habitat == "Fringing", "Fringe", Habitat))

algal_data %>% filter(Taxonomy_Substrate_Functional_Group == "No data")


alg_data_ready<-algal_data_wide %>% filter(algal_cover <= 100) # There are three rows (quadrats) where the totals are greater than 100.

# summarize to site x year
algal_data_summarized<-alg_data_ready %>%
  mutate(Habitat = ifelse(Depth == 10, '10m',
                                  ifelse(Depth == 17, '17m',
                                         Habitat))) %>%
  group_by(Year,Site,Habitat,Transect) %>%
  dplyr::summarise(mean_cover = mean(algal_cover,na.rm = T)) %>%
  group_by(Year,Site,Habitat) %>%
  dplyr::summarise(mean_cover = mean(mean_cover,na.rm = T)) %>%
  mutate(Habitat = case_when(Habitat == "Outer 10" ~ "10m",
                             Habitat == "Outer 17" ~ "17m",
                             TRUE ~ Habitat)) %>%
  ungroup()

ggplot(algal_data_summarized,aes(x = Year, y= mean_cover, color = Habitat)) +
  geom_point() +
  geom_line() +
  facet_grid(Habitat~Site) +
  theme_bw()

# Export for use in wavelet analyses
write.csv(algal_data_summarized,
          './data/summarized/MCR_alg_summarized.csv',
          row.names = F)

# END #
