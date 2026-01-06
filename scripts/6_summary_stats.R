# Summary stats for results

# Pocillopora ####
max_poc_sync<-max(poc_wmf_plotdat_df$Magnitude,na.rm = T)

# 2-5
poc_wmf_plotdat_df %>%
  filter(Timescale >= 2 & Timescale <= 5) %>%
  filter(Time < 2019) %>%
  dplyr::summarize(mean_mag = mean(Magnitude,na.rm = T)) %>%
  mutate(perc_sync = mean_mag/max_poc_sync*100)

poc_wmf_plotdat_df %>%
  filter(Timescale >= 2 & Timescale <= 5) %>%
  filter(Time >= 2019) %>%
  dplyr::summarize(mean_mag = mean(Magnitude,na.rm = T)) %>%
  mutate(perc_sync = mean_mag/max_poc_sync*100)

# 5-10
poc_wmf_plotdat_df %>%
  filter(Timescale >= 5 & Timescale <= 10) %>%
  filter(Time < 2019) %>%
  dplyr::summarize(mean_mag = mean(Magnitude,na.rm = T)) %>%
  mutate(perc_sync = mean_mag/max_poc_sync*100)

# 2-10
poc_wmf_short_plotdat_df %>%
  filter(Timescale >= 2 & Timescale <= 10) %>%
  dplyr::summarize(mean_mag = mean(Magnitude,na.rm = T)) %>%
  mutate(perc_sync = mean_mag/max_poc_sync*100)


# Porites ####
max_por_sync<-max(por_wmf_plotdat_df$Magnitude,na.rm = T)

# 2-5
por_wmf_plotdat_df %>%
  filter(Timescale >= 2 & Timescale <= 5) %>%
  filter(Time >= 2008 & Time <= 2010) %>%
  dplyr::summarize(mean_mag = mean(Magnitude,na.rm = T)) %>%
  mutate(perc_sync = mean_mag/max_por_sync*100)


por_wmf_plotdat_df %>%
  filter(Time > 2010) %>%
  dplyr::summarize(mean_mag = mean(Magnitude,na.rm = T)) %>%
  mutate(perc_sync = mean_mag/max_por_sync*100)

# 5-10
por_wmf_plotdat_df %>%
  filter(Timescale >= 5 & Timescale <= 10) %>%
  dplyr::summarize(mean_mag = mean(Magnitude,na.rm = T)) %>%
  mutate(perc_sync = mean_mag/max_por_sync*100)

# 2-10
por_wmf_plotdat_df %>%
  filter(Timescale >= 2 & Timescale <= 10) %>%
  dplyr::summarize(mean_mag = mean(Magnitude,na.rm = T)) %>%
  mutate(perc_sync = mean_mag/max_por_sync*100)


# Montipora ####
max_mont_sync<-max(mont_wmf_plotdat_df$Magnitude,na.rm = T)

# all timescales
mont_wmf_plotdat_df %>%
  filter(Time >= 2006 & Time <= 2018) %>%
  dplyr::summarize(mean_mag = mean(Magnitude,na.rm = T)) %>%
  mutate(perc_sync = mean_mag/max_mont_sync*100)

mont_wmf_plotdat_df %>%
  filter(Time >= 2019 & Time <= 2024) %>%
  dplyr::summarize(mean_mag = mean(Magnitude,na.rm = T)) %>%
  mutate(perc_sync = mean_mag/max_mont_sync*100)

# 2-5
mont_wmf_plotdat_df %>%
  filter(Time >= 2006 & Time <= 2018) %>%
  filter(Timescale >= 2 & Timescale <= 5) %>%
  dplyr::summarize(mean_mag = mean(Magnitude,na.rm = T)) %>%
  mutate(perc_sync = mean_mag/max_mont_sync*100)

# 5-10
mont_wmf_plotdat_df %>%
  filter(Timescale >= 5 & Timescale <= 10) %>%
  dplyr::summarize(mean_mag = mean(Magnitude,na.rm = T)) %>%
  mutate(perc_sync = mean_mag/max_mont_sync*100)

# 2-10
mont_wmf_plotdat_df %>%
  filter(Timescale >= 2 & Timescale <= 10) %>%
  dplyr::summarize(mean_mag = mean(Magnitude,na.rm = T)) %>%
  mutate(perc_sync = mean_mag/max_mont_sync*100)



# Acropora ####
max_acro_sync<-max(acro_wmf_plotdat_df$Magnitude,na.rm = T)

# 2-5
acro_wmf_plotdat_df %>%
  filter(Timescale >= 2 & Timescale <= 5) %>%
  filter(Time >= 2006 & Time <= 2010) %>%
  dplyr::summarize(mean_mag = mean(Magnitude,na.rm = T)) %>%
  mutate(perc_sync = mean_mag/max_acro_sync*100)

# 5-10
acro_wmf_plotdat_df %>%
  filter(Timescale >= 5 & Timescale <= 10) %>%
  dplyr::summarize(mean_mag = mean(Magnitude,na.rm = T)) %>%
  mutate(perc_sync = mean_mag/max_acro_sync*100)

# 2-10
acro_wmf_plotdat_df %>%
  filter(Timescale >= 2 & Timescale <= 10) %>%
  dplyr::summarize(mean_mag = mean(Magnitude,na.rm = T)) %>%
  mutate(perc_sync = mean_mag/max_acro_sync*100)






# Figure 2 Plot ####
# Coral community structure data for species of interest -- summarized for plotting time-series
MCR_cover_time_corespp<-read.csv('./data/summarized/MCR_CTS_replicate_updated.csv')

library(png)
library(grid)

cot <- readPNG("./figures/symbols/cots.png")
cot_grob <- rasterGrob(cot, interpolate = TRUE)

cyclone <- readPNG("./figures/symbols/cyclone.png")
cyclone_grob <- rasterGrob(cyclone, interpolate = TRUE)

therm <- readPNG("./figures/symbols/therm.png")
therm_grob <- rasterGrob(therm, interpolate = TRUE)


# Poc plot
poc_main_ts<-MCR_cover_time_corespp %>%
  filter(species %in% c("Pocillopora spp.")) %>%
  ggplot(aes(x = Date, y = mean, color = Habitat)) +
  geom_rect(aes(xmin = 2006, xmax = 2010, ymin = -Inf, ymax = Inf), color = NA, fill = "grey", alpha = 0.1) +
  annotation_custom(cot_grob, xmin = 2008 - 5, xmax = 2008 + 5, ymin = 50, ymax = 70) + # COTS outbreak in 2008
  annotation_custom(cyclone_grob, xmin = 2010 - 5, xmax = 2010 + 5, ymin = 50, ymax = 70) + # Cyclone Oli in 2010
  annotation_custom(therm_grob, xmin = 2019 - 5, xmax = 2019 + 5, ymin = 50, ymax = 70) + # bleaching in 2019
  geom_point() +
  geom_path() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0) +
  facet_grid(species ~ Site) +
  theme_bw() +
  labs(x = "Year", y = "Percent cover") +
  theme(legend.position = "top",
        plot.margin = unit(c(0,0,0,0),"cm")) +
  scale_color_manual(values = c("#EC7D6D","#00BA38","#C67CFF","#6B9EF6"))


por_main_ts<-MCR_cover_time_corespp %>%
  filter(species %in% c("Porites spp.")) %>%
  ggplot(aes(x = Date, y = mean, color = Habitat)) +
  geom_rect(aes(xmin = 2006, xmax = 2010, ymin = -Inf, ymax = Inf), color = NA, fill = "grey", alpha = 0.1) +
  geom_point() +
  geom_path() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0) +
  facet_grid(species ~ Site) +
  theme_bw() +
  labs(x = "Year", y = "Percent cover") +
  theme(legend.position = "top",
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm")) +
  scale_color_manual(values = c("#EC7D6D","#00BA38","#C67CFF","#6B9EF6"))



mont_main_ts<-MCR_cover_time_corespp %>%
  filter(species %in% c("Montipora spp.")) %>%
  ggplot(aes(x = Date, y = mean, color = Habitat)) +
  geom_rect(aes(xmin = 2006, xmax = 2010, ymin = -Inf, ymax = Inf), color = NA, fill = "grey", alpha = 0.1) +
  geom_point() +
  geom_path() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0) +
  facet_grid(species ~ Site) +
  theme_bw() +
  labs(x = "Year", y = "Percent cover") +
  theme(legend.position = "top",
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm")) +
  scale_color_manual(values = c("#EC7D6D","#00BA38","#C67CFF","#6B9EF6"))



acro_main_ts<-MCR_cover_time_corespp %>%
  filter(species %in% c("Acropora spp.")) %>%
  ggplot(aes(x = Date, y = mean, color = Habitat)) +
  geom_rect(aes(xmin = 2006, xmax = 2010, ymin = -Inf, ymax = Inf), color = NA, fill = "grey", alpha = 0.1) +
  geom_point() +
  geom_path() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0) +
  facet_grid(species ~ Site) +
  theme_bw() +
  labs(x = "Year", y = "Percent cover") +
  theme(legend.position = "top",
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm")) +
  scale_color_manual(values = c("#EC7D6D","#00BA38","#C67CFF","#6B9EF6"))



Figure2<-poc_main_ts / por_main_ts / mont_main_ts / acro_main_ts + plot_layout(guides = "collect",axes = "collect",axis_titles = "collect") & theme(legend.position = "top") & geom_vline(xintercept = 2010, linetype = "dashed", alpha = 0.5) & geom_vline(xintercept = 2019, linetype = "dashed", alpha = 0.5) #& geom_vline(xintercept = 2008, linetype = "dashed", alpha = 0.5)


ggsave("./figures/Figure_2.jpeg",
       Figure2,
       dpi = 300,
       width = 11,
       height = 5)

# END #
