# Coral data curation

# Packages
library(tidyverse)
library(reshape2)

# Forereef (10m & 17m) fringe data ####
# knb-lter-mcr.4 - Data DOI: https://doi.org/10.6073/pasta/cc668aff4dde756810631da46ff8fff8

# Data import
MCR_coral<-read.csv('./data/knb-lter-mcr.4_2_20250409.csv') # wide format data

# Data cleaning
MCR_coral$Location <- as.factor(MCR_coral$Location)
MCR_coral_num<-MCR_coral %>% mutate_if(is.character,as.numeric)


# Need to aggregate Porites due to changes in taxonomy
# this issue occurs prior to 2011
MCR_post_2011<-MCR_coral_num %>% filter(Date >= 2011)
MCR_pre_2011<-MCR_coral_num %>% filter(Date < 2011)

MCR_pre_2011_upd<-MCR_pre_2011 %>%
  mutate(Porites_new = rowSums(select(., Porites, Porites_irregularis, Porites_rus, Porites_spp_Massive), na.rm = TRUE))

MCR_pre_2011_upd_check<-MCR_pre_2011_upd %>%
  mutate(diff_por = if_else(Porites == Porites_new,"good","NO!")) %>% # all good
  filter(diff_por == "NO!")

MCR_post_2011_upd<-MCR_post_2011 %>%
  mutate(Porites_new = rowSums(select(., Porites, Porites_irregularis, Porites_rus, Porites_spp_Massive), na.rm = TRUE))

MCR_pre_2011_upd_check<-MCR_pre_2011_upd %>%
  mutate(diff_por = if_else(Porites == Porites_new,"good","NO!")) %>% # all good
  filter(diff_por == "NO!")

# Now we can drop all Porites columns except for Porites_new, then rbind back
MCR_post_2011_upd$Porites <- NULL
MCR_post_2011_upd$Porites_irregularis <- NULL
MCR_post_2011_upd$Porites_rus <- NULL
MCR_post_2011_upd$Porites_spp_Massive <- NULL

MCR_pre_2011_upd$Porites <- NULL
MCR_pre_2011_upd$Porites_irregularis <- NULL
MCR_pre_2011_upd$Porites_rus<- NULL
MCR_pre_2011_upd$Porites_spp_Massive<- NULL


MCR_almost_ready<-rbind(MCR_pre_2011_upd,MCR_post_2011_upd)

# Remove a few extra classifiers that arent relevant here
MCR_almost_ready$Non_Coralline_Crustose_Algae <- NULL
MCR_almost_ready$Sand <- NULL
MCR_almost_ready$CTB <- NULL

# Add factors for habitat, transect, and quadrat
MCR_almost_almost_ready<-MCR_almost_ready %>%
  mutate(
    Site = case_when( # habitat deliniator
      str_detect(Location, "LTER 1") ~ "LTER 1",
      str_detect(Location, "LTER 2") ~ "LTER 2",
      str_detect(Location, "LTER 3") ~ "LTER 3",
      str_detect(Location, "LTER 4") ~ "LTER 4",
      str_detect(Location, "LTER 5") ~ "LTER 5",
      str_detect(Location, "LTER 6") ~ "LTER 6",
      TRUE ~ NA_character_),
    Habitat = case_when( # habitat deliniator
      str_detect(Location, "Fringing") ~ "Fringe",
      str_detect(Location, "10 m") ~ "10m",
      str_detect(Location, "17 m") ~ "17m",
      TRUE ~ NA_character_),
    Transect = case_when( # habitat deliniator
      str_detect(Location, "1-2") ~ "1",
      str_detect(Location, "2-3") ~ "2",
      str_detect(Location, "3-4") ~ "3",
      str_detect(Location, "4-5") ~ "4",
      str_detect(Location, "5-6") ~ "5",
      TRUE ~ NA_character_),
    Quadrat = case_when( # habitat deliniator
      str_detect(Location, "Quadrat 1") ~ "1",
      str_detect(Location, "Quadrat 2") ~ "2",
      str_detect(Location, "Quadrat 3") ~ "3",
      str_detect(Location, "Quadrat 4") ~ "4",
      str_detect(Location, "Quadrat 5") ~ "5",
      str_detect(Location, "Quadrat 6") ~ "6",
      str_detect(Location, "Quadrat 7") ~ "7",
      str_detect(Location, "Quadrat 8") ~ "8",
      TRUE ~ NA_character_)
  )


MCR_ready<-MCR_almost_almost_ready %>%
  select(-Location) %>%
  mutate_all(~ifelse(is.na(.), 0, .))


MCR_ready_site_level <- MCR_ready %>%
  group_by(Site,Date,Transect, Habitat) %>%
  summarize_all(mean,na.rm =T) %>%
  group_by(Site,Date,Habitat) %>%
  summarize_all(mean,na.rm =T) %>%
  select(-c(Transect,Quadrat))

MCR_ready_tran_level <- MCR_ready %>%
  group_by(Site,Date,Transect, Habitat) %>%
  summarize_all(mean,na.rm =T) %>%
  select(-Quadrat)

MCR_ready_rep_level <- MCR_ready %>%
  group_by(Site,Date,Transect, Quadrat, Habitat) %>%
  summarize_all(mean,na.rm =T)



# Backreef data ####
knb_1038_lagoon_srednick<-read.csv('./data/knb-lter-mcr.1038_srednick_manual_fullres.csv')


# Joining and summarizing ####
# Join knb-lter-mcr.4 and knb-lter-mcr.1038 raw photos annotated by G. Srednick
# Filtered for genera of interest: Pocillopora, Porites, Acropora, Montipora

lagoon_forjoin<-knb_1038_lagoon_srednick %>%
  dplyr::mutate(Habitat = "Backreef") %>%
  dplyr::select(Site,Habitat,Date,Poles,Transect,Pocillopora,Porites_new,Acropora,Montipora) %>%
  dplyr::group_by(Site,Habitat,Date,Poles,Transect) %>%
  dplyr::summarize_if(is.numeric,mean,na.rm=T) %>%
  dplyr::group_by(Site,Habitat,Date,Poles) %>%
  dplyr::summarize_if(is.numeric,mean,na.rm=T) %>%
  dplyr::group_by(Site,Habitat,Date) %>%
  dplyr::summarize_if(is.numeric,mean,na.rm=T)

other_sites_forjoin<-MCR_ready %>%
  dplyr::select(Site,Habitat,Date,Transect,Quadrat,Pocillopora,Porites_new,Acropora,Montipora) %>%
  dplyr::group_by(Site,Date,Transect, Habitat) %>%
  dplyr::summarize_all(mean,na.rm =T) %>%
  dplyr::group_by(Site,Date,Habitat) %>%
  dplyr::summarize_all(mean,na.rm =T) %>%
  dplyr::select(-c(Transect,Quadrat)) %>%
  dplyr::mutate(Site = str_replace(Site," ","0"))

complete_MCR_coral_data_site<-rbind(other_sites_forjoin,lagoon_forjoin)

# Export -- for wavelet analyses
write.csv(complete_MCR_coral_data_site,'./data/MCR_CTS_site_updated.csv',row.names = F)


# Plotting time series ####
lagoon_forjoin_rep<-knb_1038_lagoon_srednick %>%
  dplyr::group_by(Site,Date,Poles,Transect) %>%
  dplyr::mutate(Quad_id = row_number()) %>%
  dplyr::mutate(Habitat = "Backreef",
         bommie_tran = paste0(Poles,"_",Transect)) %>%
  dplyr::select(Site,Habitat,Date,Poles,bommie_tran,Transect,Quad_id,
                Pocillopora,Porites_new,Acropora,Montipora) %>%
  dplyr::mutate(Transect = bommie_tran)

lagoon_forjoin_rep$bommie_tran<-NULL

# Time series rep check
lag_rep_check<-lagoon_forjoin_rep %>%
  group_by(Site,Date) %>%
  dplyr::summarize(n())


other_sites_forjoin_rep<-MCR_ready %>%
  dplyr::select(Site,Habitat,Date,Transect,Quadrat, Pocillopora,Porites_new,Acropora,Montipora) %>%
  mutate(Site = str_replace(Site," ","0"))

complete_MCR_coral_data_full<-plyr::rbind.fill(other_sites_forjoin_rep,lagoon_forjoin_rep)
MCR_cover_time_new<-complete_MCR_coral_data_full

MCR_cover_time_cur<-MCR_cover_time_new %>%
  dplyr::ungroup() %>%  # Remove any grouping before summarizing
  dplyr::rename("Pocillopora spp." = Pocillopora,
                "Acropora spp." = Acropora,
                "Porites spp." = Porites_new,
                "Montipora spp." = Montipora)



MCR_cover_time_cur_long<-MCR_cover_time_cur %>%
  pivot_longer(col = -c(Transect,Quadrat,Date,Habitat,Site,Poles), names_to = "species",values_to = "cover")

MCR_cover_time_cur_long$Habitat<-factor(MCR_cover_time_cur_long$Habitat, levels = c("Fringe","10m","17m","Backreef"))

se_fn <- function(x, na.rm = FALSE) {
  if (na.rm) x <- na.omit(x)
  sd(x) / sqrt(length(x))
}

MCR_cover_time_summarized<-MCR_cover_time_cur_long %>%
  group_by(Site,Habitat,Date,species) %>%
  dplyr::summarize(mean = mean(cover,na.rm = T),
                   se = se_fn(cover),
                   n = n()) %>%
  filter(Date > 2005)

MCR_cover_time_corespp<-MCR_cover_time_summarized %>% # we want full time series for plotting
  filter(species %in% c("Pocillopora spp.", "Acropora spp.", "Montipora spp.","Porites spp.","other corals"))

MCR_cover_time_corespp$species<-factor(MCR_cover_time_corespp$species, levels = c("Acropora spp.","Pocillopora spp.",  "Porites spp.","Montipora spp.","other corals"))

coral_colors <- c("Acropora spp." = "#31688e", "Pocillopora spp." = "#f8765c", "Porites spp." = "#982d80","Montipora spp." = "pink", "other corals" = "#5DC963")

# Save products to csv -- plotting in '6_summary_stats.R"
write.csv(MCR_cover_time_corespp,'./data/summarized/MCR_CTS_replicate_updated.csv',row.names = F)

min_reps<-MCR_cover_time_corespp %>%
  group_by(Site,Date,Habitat) %>%
  summarize(min_rep = min(n))


# END #
