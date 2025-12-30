# Portfolio effect analyses and plotting

# Packages
library(effectsize)
library(parameters)


# Helper functions ####
cv <- function(x){
  sd(x)/abs(mean(x))
}

CV_summary <- function(site1, site2, genus,df) {
  df1 <- df %>%
    filter(Site_hab == site1) %>%
    dplyr::select(Date, genus) %>%
    rename(Poc1 = genus)

  df2 <- df %>%
    filter(Site_hab == site2) %>%
    dplyr::select(Date, genus) %>%
    rename(Poc2 = genus)

  # Join by date
  joined <- inner_join(df1, df2, by = "Date")

  if (nrow(joined) > 1) {
    # Individual CVs
    cv1 <- sd(joined$Poc1, na.rm = TRUE) / mean(joined$Poc1, na.rm = TRUE)
    cv2 <- sd(joined$Poc2, na.rm = TRUE) / mean(joined$Poc2, na.rm = TRUE)
    mean_cv <- mean(c(cv1, cv2), na.rm = TRUE)

    # CV of the summed signal
    summed <- joined$Poc1 + joined$Poc2
    cv_sum <- sd(summed, na.rm = TRUE) / mean(summed, na.rm = TRUE)

    return(tibble(cv1 = cv1, cv2 = cv2, mean_cv = mean_cv, cv_sum = cv_sum))
  } else {
    return(tibble(cv1 = NA_real_, cv2 = NA_real_, mean_cv = NA_real_, cv_sum = NA_real_))
  }
}



# Model run for PE~synchrony
get_pairwise_comparisons <- function(model, hierarchies = c("habitat_pair","site_pair","shore_pair"), timescale_val = NULL) {

  map_dfr(hierarchies, function(h) {

    emm <- emmeans(model, as.formula(paste0("~ ", h)))
    pairs <- contrast(emm, method = "revpairwise") %>%
      summary(infer = TRUE) %>%
      data.frame()

    pairs %>%
      mutate(
        hierarchy = gsub("_pair", "", h),  # clean up name
        timescale = timescale_val
      )
  })
}


# Pocillopora ####
## Estimate pairwise CV ####
poc_cv_results <- poc_predictor_sync_mats %>%
  rowwise() %>%
  mutate(tmp = list(CV_summary(Site1, Site2, poc_df,genus = "Pocillopora"))) %>%
  unnest(cols = c(tmp)) %>%
  mutate(PE = mean_cv/cv_sum) %>%
  ungroup()



poc_cv_results_grouped <- poc_cv_results %>%
  mutate(
    hab1 = str_extract(Site1, "Fringe|Backreef|10m|17m"),
    hab2 = str_extract(Site2, "Fringe|Backreef|10m|17m"),
    site1 = str_extract(Site1, "LTER01|LTER02|LTER03|LTER04|LTER05|LTER06"),
    site2 = str_extract(Site2, "LTER01|LTER02|LTER03|LTER04|LTER05|LTER06"),
    shore1 = case_when(site1 %in% c("LTER01", "LTER02") ~ "North",
                       site1 %in% c("LTER03", "LTER04") ~ "East",
                       site1 %in% c("LTER05", "LTER06") ~ "West"),
    shore2 = case_when(site2 %in% c("LTER01", "LTER02") ~ "North",
                       site2 %in% c("LTER03", "LTER04") ~ "East",
                       site2 %in% c("LTER05", "LTER06") ~ "West"),

    habitat_pair = case_when(
      hab1 == hab2 ~ "within",
      hab1 != hab2 ~ "between",
      TRUE ~ NA_character_
    ),
    site_pair = case_when(
      site1 == site2 ~ "within",
      site1 != site2 ~ "between",
      TRUE ~ NA_character_
    ),
    shore_pair = case_when(
      shore1 == shore2 ~ "within",
      shore1 != shore2 ~ "between",
      TRUE ~ NA_character_
    )
  )
poc_sync_base_lm <- lmer(log1p(PE) ~ POC_sync + (1 | Site1) + (1 | Site2), data = poc_cv_results_grouped)
summary(poc_sync_base_lm)
poc_sync_base_lm_df<-ggpredict(poc_sync_base_lm, terms = c("POC_sync"), back_transform = F) %>% as.data.frame()

poc_cv_results_grouped_resid <- poc_cv_results_grouped %>%
  mutate(resid_sync_only = residuals(poc_sync_base_lm))

poc_sync_base_lm_r2_vals <- r.squaredGLMM(poc_sync_base_lm)

poc_sync_base_lm_r2_vals[1]
poc_sync_base_lm_r2_vals[2]

poc_annotate_label <- paste0(
  "<span style='font-size:10pt'><i>p</i> &lt; 0.001<br>",
  "<i>R</i><sup>2<sub>m</sub></sup> = 0.08<br>",
  "<i>R</i><sup>2<sub>c</sub></sup> = 0.75</span>"
)



# habitat
summary(aov(resid_sync_only ~ habitat_pair*timescale, data = poc_cv_results_grouped_resid))

# Site pair
summary(aov(resid_sync_only ~ site_pair*timescale, data = poc_cv_results_grouped_resid))

# Shore pair
summary(aov(resid_sync_only ~ shore_pair*timescale, data = poc_cv_results_grouped_resid))

eta_squared(aov(resid_sync_only ~ habitat_pair*timescale, data = poc_cv_results_grouped_resid))
eta_squared(aov(resid_sync_only ~ site_pair*timescale, data = poc_cv_results_grouped_resid))
eta_squared(aov(resid_sync_only ~ shore_pair*timescale, data = poc_cv_results_grouped_resid))

# Estimating pairwise differences
poc_hab_pairwise_lm<-lm(resid_sync_only ~ habitat_pair*timescale, data = poc_cv_results_grouped_resid)
poc_site_pairwise_lm<-lm(resid_sync_only ~ site_pair*timescale, data = poc_cv_results_grouped_resid)
poc_shore_pairwise_lm<-lm(resid_sync_only ~ shore_pair*timescale, data = poc_cv_results_grouped_resid)

emm_poc_hab_resid <- emmeans(poc_hab_pairwise_lm, pairwise ~ habitat_pair | timescale)
emm_poc_site_resid <- emmeans(poc_site_pairwise_lm, pairwise ~ site_pair | timescale)
emm_poc_shore_resid <- emmeans(poc_shore_pairwise_lm, pairwise ~ shore_pair | timescale)

poc_hab_resid_pairs <- contrast(emm_poc_hab_resid, method = "revpairwise") %>%
  summary(infer = TRUE) %>%
  data.frame() %>%
  mutate(hierarchy = "habitat")

poc_site_resid_pairs <- contrast(emm_poc_site_resid, method = "revpairwise") %>%
  summary(infer = TRUE) %>%
  data.frame() %>%
  mutate(hierarchy = "site")

poc_shore_resid_pairs <- contrast(emm_poc_shore_resid, method = "revpairwise") %>%
  summary(infer = TRUE) %>%
  data.frame() %>%
  mutate(hierarchy = "shore")


poc_sig_comparisons<-rbind(
  poc_hab_resid_pairs,
  poc_site_resid_pairs,
  poc_shore_resid_pairs) %>%
  mutate(sig  = p.value < 0.05) %>%
  mutate(timescale_written = case_when(
    timescale == "short"  ~ "2-5y",
    timescale == "medium" ~ "5-10y",
    timescale == "long"   ~ "2-10y"
  ))



poc_resid_long_hab<-poc_cv_results_grouped_resid %>%
  select(timescale, habitat_pair, PE, resid_sync_only) %>%
  rename(habitat = habitat_pair) %>%
  pivot_longer(cols = c("habitat"), names_to = "hierarchy", values_to = "level")

poc_resid_long_site<-poc_cv_results_grouped_resid %>%
  select(timescale, site_pair, PE, resid_sync_only) %>%
  rename(site = site_pair) %>%
  pivot_longer(cols = c("site"), names_to = "hierarchy", values_to = "level")

poc_resid_long_shore<-poc_cv_results_grouped_resid %>%
  select(timescale, shore_pair, PE, resid_sync_only) %>%
  rename(shore = shore_pair) %>%
  pivot_longer(cols = c("shore"), names_to = "hierarchy", values_to = "level")

poc_resid_long_complete<-rbind(poc_resid_long_hab,poc_resid_long_site,poc_resid_long_shore) %>%
  mutate(timescale_written = case_when(
    timescale == "short"  ~ "2-5y",
    timescale == "medium" ~ "5-10y",
    timescale == "long"   ~ "2-10y"
  ))

poc_resid_long_complete$timescale_written = factor(poc_resid_long_complete$timescale_written,
                                                   levels=c("2-5y","5-10y","2-10y"))


poc_asterisks_df <- poc_resid_long_complete %>%
  group_by(hierarchy, timescale_written) %>%
  dplyr::summarise(y = mean(resid_sync_only) + 0.02, .groups = "drop") %>%
  left_join(poc_sig_comparisons %>%
              filter(sig == TRUE) %>%
              select(hierarchy, timescale_written, sig),
            by = c("hierarchy", "timescale_written")) %>%
  na.omit()

poc_sync_base_lm_r2_vals
poc_annotate_label <- paste0(
  "<span style='font-size:10pt'><i>p</i> &lt; 0.001<br>",
  "<i>R</i><sup>2<sub>m</sub></sup> = 0.08<br>",
  "<i>R</i><sup>2<sub>c</sub></sup> = 0.75</span>"
)



poc_asterisks_df$timescale_written = factor(poc_asterisks_df$timescale_written,
                                             levels=c("2-5y","5-10y","2-10y"))

poc_sync_base_lm_plot<-ggplot(poc_sync_base_lm_df, aes(x = x, y = predicted)) +
  geom_point(data = poc_cv_results_grouped, aes(x = POC_sync, y = log1p(PE)), alpha = 0.1) +
  geom_ribbon(aes(ymin = conf.low,ymax = conf.high), alpha = 0.4, fill = "blue") +
  geom_line(color = "blue") +
  theme_bw() +
  facet_grid(~"PE~synchrony") +
  labs(x = "Wavelet synchrony", y = "ln(Portfolio effects + 1)") +
  ggtext::geom_richtext(
    aes(x = 1, y = 1.5, label = poc_annotate_label),
    fill = NA, label.color = NA,
    hjust = 1, size = 4, lineheight = 1.2
  )

poc_resid_plot<-ggplot(poc_resid_long_complete, aes(x = hierarchy, y = resid_sync_only, color = level)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_summary(
    fun.data = mean_se,
    geom = "pointrange",
    position = position_dodge(width = 0.6),
    size = 0.5
  ) +
  facet_grid("Pocillopora"~timescale_written) +
  theme_bw() +
  theme(strip.text.y = element_text(face = "italic")) +
  scale_y_continuous(position = "right") +
  labs(y = "PE ~ synchrony residuals", x = "Hierarchy Level") +
  scale_color_manual(values = c("within" = "black", "between" = "grey")) +
  geom_text(
    data = poc_asterisks_df,
    aes(x = hierarchy, y = y, label = "*"),
    color = "red",
    size = 6,
    inherit.aes = FALSE
  )



poc_combined_pe_lm_plot<-poc_sync_base_lm_plot + poc_resid_plot + plot_layout(widths = c(2,3))


## Porites ####
por_cv_results <- POR_predictor_sync_mats %>%
  rowwise() %>%
  mutate(tmp = list(CV_summary(Site1, Site2, por_df,genus = "Porites_new"))) %>%
  unnest(cols = c(tmp)) %>%
  mutate(PE = mean_cv/cv_sum) %>%
  ungroup()



por_cv_results_grouped <- por_cv_results %>%
  mutate(
    hab1 = str_extract(Site1, "Fringe|Backreef|10m|17m"),
    hab2 = str_extract(Site2, "Fringe|Backreef|10m|17m"),
    site1 = str_extract(Site1, "LTER01|LTER02|LTER03|LTER04|LTER05|LTER06"),
    site2 = str_extract(Site2, "LTER01|LTER02|LTER03|LTER04|LTER05|LTER06"),
    shore1 = case_when(site1 %in% c("LTER01", "LTER02") ~ "North",
                       site1 %in% c("LTER03", "LTER04") ~ "East",
                       site1 %in% c("LTER05", "LTER06") ~ "West"),
    shore2 = case_when(site2 %in% c("LTER01", "LTER02") ~ "North",
                       site2 %in% c("LTER03", "LTER04") ~ "East",
                       site2 %in% c("LTER05", "LTER06") ~ "West"),

    habitat_pair = case_when(
      hab1 == hab2 ~ "within",
      hab1 != hab2 ~ "between",
      TRUE ~ NA_character_
    ),
    site_pair = case_when(
      site1 == site2 ~ "within",
      site1 != site2 ~ "between",
      TRUE ~ NA_character_
    ),
    shore_pair = case_when(
      shore1 == shore2 ~ "within",
      shore1 != shore2 ~ "between",
      TRUE ~ NA_character_
    )
  )

por_sync_base_lm <- lmer(log1p(PE) ~ POR_sync + (1 | Site1) + (1 | Site2), data = por_cv_results_grouped)
summary(por_sync_base_lm)
por_sync_base_lm_df<-ggpredict(por_sync_base_lm, terms = c("POR_sync"), back_transform = F) %>% as.data.frame()
por_cv_results_grouped_resid <- por_cv_results_grouped %>%
  mutate(resid_sync_only = residuals(por_sync_base_lm))

por_sync_base_lm_r2_vals <- r.squaredGLMM(por_sync_base_lm)

por_sync_base_lm_r2_vals[1]
por_sync_base_lm_r2_vals[2]

por_annotate_label <- paste0(
  "<span style='font-size:10pt'><i>p</i> &lt; 0.001<br>",
  "<i>R</i><sup>2<sub>m</sub></sup> = 0.02<br>",
  "<i>R</i><sup>2<sub>c</sub></sup> = 0.64</span>"
)


# habitat
summary(aov(resid_sync_only ~ habitat_pair*timescale, data = por_cv_results_grouped_resid))

# Site pair
summary(aov(resid_sync_only ~ site_pair*timescale, data = por_cv_results_grouped_resid))

# Shore pair
summary(aov(resid_sync_only ~ shore_pair*timescale, data = por_cv_results_grouped_resid))

eta_squared(aov(resid_sync_only ~ habitat_pair*timescale, data = por_cv_results_grouped_resid))
eta_squared(aov(resid_sync_only ~ site_pair*timescale, data = por_cv_results_grouped_resid))
eta_squared(aov(resid_sync_only ~ shore_pair*timescale, data = por_cv_results_grouped_resid))

por_resid_long_hab<-por_cv_results_grouped_resid %>%
  select(timescale, habitat_pair, PE, resid_sync_only) %>%
  rename(habitat = habitat_pair) %>%
  pivot_longer(cols = c("habitat"), names_to = "hierarchy", values_to = "level")

por_resid_long_site<-por_cv_results_grouped_resid %>%
  select(timescale, site_pair, PE, resid_sync_only) %>%
  rename(site = site_pair) %>%
  pivot_longer(cols = c("site"), names_to = "hierarchy", values_to = "level")

por_resid_long_shore<-por_cv_results_grouped_resid %>%
  select(timescale, shore_pair, PE, resid_sync_only) %>%
  rename(shore = shore_pair) %>%
  pivot_longer(cols = c("shore"), names_to = "hierarchy", values_to = "level")

por_resid_long_complete<-rbind(por_resid_long_hab,por_resid_long_site,por_resid_long_shore) %>%
  mutate(timescale_written = case_when(
    timescale == "short"  ~ "2-5y",
    timescale == "medium" ~ "5-10y",
    timescale == "long"   ~ "2-10y"
  ))

por_resid_long_complete$timescale_written = factor(por_resid_long_complete$timescale_written,
                                       levels=c("2-5y","5-10y","2-10y"))

# Estimating pairwise differences
por_hab_pairwise_lm<-lm(resid_sync_only ~ habitat_pair*timescale, data = por_cv_results_grouped_resid)
por_site_pairwise_lm<-lm(resid_sync_only ~ site_pair*timescale, data = por_cv_results_grouped_resid)
por_shore_pairwise_lm<-lm(resid_sync_only ~ shore_pair*timescale, data = por_cv_results_grouped_resid)

emm_por_hab_resid <- emmeans(por_hab_pairwise_lm, pairwise ~ habitat_pair | timescale)
emm_por_site_resid <- emmeans(por_site_pairwise_lm, pairwise ~ site_pair | timescale)
emm_por_shore_resid <- emmeans(por_shore_pairwise_lm, pairwise ~ shore_pair | timescale)

por_hab_resid_pairs <- contrast(emm_por_hab_resid, method = "revpairwise") %>%
  summary(infer = TRUE) %>%
  data.frame() %>%
  mutate(hierarchy = "habitat")

por_site_resid_pairs <- contrast(emm_por_site_resid, method = "revpairwise") %>%
  summary(infer = TRUE) %>%
  data.frame() %>%
  mutate(hierarchy = "site")

por_shore_resid_pairs <- contrast(emm_por_shore_resid, method = "revpairwise") %>%
  summary(infer = TRUE) %>%
  data.frame() %>%
  mutate(hierarchy = "shore")

por_sig_comparisons<-rbind(
  por_hab_resid_pairs,
  por_site_resid_pairs,
  por_shore_resid_pairs) %>%
  mutate(sig  = p.value < 0.05) %>%
  mutate(timescale_written = case_when(
    timescale == "short"  ~ "2-5y",
    timescale == "medium" ~ "5-10y",
    timescale == "long"   ~ "2-10y"
  ))

por_asterisks_df <- por_resid_long_complete %>%
  group_by(hierarchy, timescale_written) %>%
  dplyr::summarise(y = mean(resid_sync_only) + 0.02, .groups = "drop") %>%
  left_join(por_sig_comparisons %>%
              filter(sig == TRUE) %>%
              select(hierarchy, timescale_written, sig),
            by = c("hierarchy", "timescale_written")) %>%
  na.omit()

por_asterisks_df$timescale_written = factor(por_asterisks_df$timescale_written,
                                             levels=c("2-5y","5-10y","2-10y"))

por_sync_base_lm_plot<-ggplot(por_sync_base_lm_df, aes(x = x, y = predicted)) +
  geom_point(data = por_cv_results_grouped, aes(x = POR_sync, y = log1p(PE)), alpha = 0.1) +
  geom_ribbon(aes(ymin = conf.low,ymax = conf.high), alpha = 0.4, fill = "blue") +
  geom_line(color = "blue") +
  theme_bw() +
  facet_grid(~"PE~synchrony") +
  labs(x = "Wavelet synchrony", y = "ln(Portfolio effects + 1)") +
  ggtext::geom_richtext(
    aes(x = 1, y = 1.6, label = por_annotate_label),
    fill = NA, label.color = NA,
    hjust = 1, size = 4, lineheight = 1.2
  )


por_resid_plot<-ggplot(por_resid_long_complete, aes(x = hierarchy, y = resid_sync_only, color = level)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_summary(
    fun.data = mean_se,
    geom = "pointrange",
    position = position_dodge(width = 0.6),
    size = 0.5
  ) +
  facet_grid("Porites"~timescale_written) +
  theme_bw() +
  theme(strip.text.y = element_text(face = "italic")) +
  scale_y_continuous(position = "right") +
  labs(y = "PE ~ synchrony residuals", x = "Hierarchy Level") +
  scale_color_manual(values = c("within" = "black", "between" = "grey")) +
  geom_text(
    data = por_asterisks_df,
    aes(x = hierarchy, y = y, label = "*"),
    color = "red",
    size = 6,
    inherit.aes = FALSE
  )


por_combined_pe_lm_plot<-por_sync_base_lm_plot + por_resid_plot + plot_layout(widths = c(2,3))

poc_combined_pe_lm_plot/por_combined_pe_lm_plot



## Montipora ####
mont_cv_results <- MONT_predictor_sync_mats %>%
  rowwise() %>%
  mutate(tmp = list(CV_summary(Site1, Site2, Mont_df,genus = "Montipora"))) %>%
  unnest(cols = c(tmp)) %>%
  mutate(PE = mean_cv/cv_sum) %>%
  ungroup()



mont_cv_results_grouped <- mont_cv_results %>%
  mutate(
    hab1 = str_extract(Site1, "Fringe|Backreef|10m|17m"),
    hab2 = str_extract(Site2, "Fringe|Backreef|10m|17m"),
    site1 = str_extract(Site1, "LTER01|LTER02|LTER03|LTER04|LTER05|LTER06"),
    site2 = str_extract(Site2, "LTER01|LTER02|LTER03|LTER04|LTER05|LTER06"),
    shore1 = case_when(site1 %in% c("LTER01", "LTER02") ~ "North",
                       site1 %in% c("LTER03", "LTER04") ~ "East",
                       site1 %in% c("LTER05", "LTER06") ~ "West"),
    shore2 = case_when(site2 %in% c("LTER01", "LTER02") ~ "North",
                       site2 %in% c("LTER03", "LTER04") ~ "East",
                       site2 %in% c("LTER05", "LTER06") ~ "West"),

    habitat_pair = case_when(
      hab1 == hab2 ~ "within",
      hab1 != hab2 ~ "between",
      TRUE ~ NA_character_
    ),
    site_pair = case_when(
      site1 == site2 ~ "within",
      site1 != site2 ~ "between",
      TRUE ~ NA_character_
    ),
    shore_pair = case_when(
      shore1 == shore2 ~ "within",
      shore1 != shore2 ~ "between",
      TRUE ~ NA_character_
    )
  )
mont_sync_base_lm <- lmer(log1p(PE) ~ MONT_sync + (1 | Site1) + (1 | Site2), data = mont_cv_results_grouped)
summary(mont_sync_base_lm)
mont_sync_base_lm_df<-ggpredict(mont_sync_base_lm, terms = c("MONT_sync"), back_transform = F) %>% as.data.frame()
mont_cv_results_grouped_resid <- mont_cv_results_grouped %>%
  mutate(resid_sync_only = residuals(mont_sync_base_lm))

mont_sync_base_lm_r2_vals <- r.squaredGLMM(mont_sync_base_lm)

mont_sync_base_lm_r2_vals[1]
mont_sync_base_lm_r2_vals[2]

mont_annotate_label <- paste0(
  "<span style='font-size:10pt'><i>p</i> &lt; 0.001<br>",
  "<i>R</i><sup>2<sub>m</sub></sup> = 0.04<br>",
  "<i>R</i><sup>2<sub>c</sub></sup> = 0.61</span>"
)



# habitat
summary(aov(resid_sync_only ~ habitat_pair*timescale, data = mont_cv_results_grouped_resid))

# Site pair
summary(aov(resid_sync_only ~ site_pair*timescale, data = mont_cv_results_grouped_resid))

# Shore pair
summary(aov(resid_sync_only ~ shore_pair*timescale, data = mont_cv_results_grouped_resid))

eta_squared(aov(resid_sync_only ~ habitat_pair*timescale, data = mont_cv_results_grouped_resid))
eta_squared(aov(resid_sync_only ~ site_pair*timescale, data = mont_cv_results_grouped_resid))
eta_squared(aov(resid_sync_only ~ shore_pair*timescale, data = mont_cv_results_grouped_resid))

mont_resid_long_hab<-mont_cv_results_grouped_resid %>%
  select(timescale, habitat_pair, PE, resid_sync_only) %>%
  rename(habitat = habitat_pair) %>%
  pivot_longer(cols = c("habitat"), names_to = "hierarchy", values_to = "level")

mont_resid_long_site<-mont_cv_results_grouped_resid %>%
  select(timescale, site_pair, PE, resid_sync_only) %>%
  rename(site = site_pair) %>%
  pivot_longer(cols = c("site"), names_to = "hierarchy", values_to = "level")

mont_resid_long_shore<-mont_cv_results_grouped_resid %>%
  select(timescale, shore_pair, PE, resid_sync_only) %>%
  rename(shore = shore_pair) %>%
  pivot_longer(cols = c("shore"), names_to = "hierarchy", values_to = "level")

mont_resid_long_complete<-rbind(mont_resid_long_hab,mont_resid_long_site,mont_resid_long_shore) %>%
  mutate(timescale_written = case_when(
    timescale == "short"  ~ "2-5y",
    timescale == "medium" ~ "5-10y",
    timescale == "long"   ~ "2-10y"
  ))

mont_resid_long_complete$timescale_written = factor(mont_resid_long_complete$timescale_written,
                                                   levels=c("2-5y","5-10y","2-10y"))

# Estimating pairwise differences
mont_hab_pairwise_lm<-lm(resid_sync_only ~ habitat_pair*timescale, data = mont_cv_results_grouped_resid)
mont_site_pairwise_lm<-lm(resid_sync_only ~ site_pair*timescale, data = mont_cv_results_grouped_resid)
mont_shore_pairwise_lm<-lm(resid_sync_only ~ shore_pair*timescale, data = mont_cv_results_grouped_resid)

emm_mont_hab_resid <- emmeans(mont_hab_pairwise_lm, pairwise ~ habitat_pair | timescale)
emm_mont_site_resid <- emmeans(mont_site_pairwise_lm, pairwise ~ site_pair | timescale)
emm_mont_shore_resid <- emmeans(mont_shore_pairwise_lm, pairwise ~ shore_pair | timescale)

mont_hab_resid_pairs <- contrast(emm_mont_hab_resid, method = "revpairwise") %>%
  summary(infer = TRUE) %>%
  data.frame() %>%
  mutate(hierarchy = "habitat")

mont_site_resid_pairs <- contrast(emm_mont_site_resid, method = "revpairwise") %>%
  summary(infer = TRUE) %>%
  data.frame() %>%
  mutate(hierarchy = "site")

mont_shore_resid_pairs <- contrast(emm_mont_shore_resid, method = "revpairwise") %>%
  summary(infer = TRUE) %>%
  data.frame() %>%
  mutate(hierarchy = "shore")

mont_sig_comparisons<-rbind(
  mont_hab_resid_pairs,
  mont_site_resid_pairs,
  mont_shore_resid_pairs) %>%
  mutate(sig  = p.value < 0.05) %>%
  mutate(timescale_written = case_when(
    timescale == "short"  ~ "2-5y",
    timescale == "medium" ~ "5-10y",
    timescale == "long"   ~ "2-10y"
  ))

mont_asterisks_df <- mont_resid_long_complete %>%
  group_by(hierarchy, timescale_written) %>%
  dplyr::summarise(y = mean(resid_sync_only) + 0.04, .groups = "drop") %>%
  left_join(mont_sig_comparisons %>%
              filter(sig == TRUE) %>%
              select(hierarchy, timescale_written, sig),
            by = c("hierarchy", "timescale_written")) %>%
  na.omit()

mont_asterisks_df$timescale_written = factor(mont_asterisks_df$timescale_written,
                                            levels=c("2-5y","5-10y","2-10y"))


mont_sync_base_lm_plot<-ggplot(mont_sync_base_lm_df, aes(x = x, y = predicted)) +
  geom_point(data = mont_cv_results_grouped, aes(x = MONT_sync, y = log1p(PE)), alpha = 0.1) +
  geom_ribbon(aes(ymin = conf.low,ymax = conf.high), alpha = 0.4, fill = "blue") +
  geom_line(color = "blue") +
  theme_bw() +
  facet_grid(~"PE~synchrony") +
  labs(x = "Wavelet synchrony", y = "ln(Portfolio effects + 1)") +
  ggtext::geom_richtext(
    aes(x = 1, y = 1.6, label = mont_annotate_label),
    fill = NA, label.color = NA,
    hjust = 1, size = 4, lineheight = 1.2
  )
mont_resid_plot<-ggplot(mont_resid_long_complete, aes(x = hierarchy, y = resid_sync_only, color = level)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_summary(
    fun.data = mean_se,
    geom = "pointrange",
    position = position_dodge(width = 0.6),
    size = 0.5
  ) +
  facet_grid("Montipora"~timescale_written) +
  theme_bw() +
  theme(strip.text.y = element_text(face = "italic")) +
  scale_y_continuous(position = "right") +
  labs(y = "PE ~ synchrony residuals", x = "Hierarchy Level") +
  scale_color_manual(values = c("within" = "black", "between" = "grey")) +
  geom_text(
    data = mont_asterisks_df,
    aes(x = hierarchy, y = y, label = "*"),
    color = "red",
    size = 6,
    inherit.aes = FALSE
  )



mont_combined_pe_lm_plot<-mont_sync_base_lm_plot + mont_resid_plot + plot_layout(widths = c(2,3))




## Acropora ####
acro_cv_results <- ACRO_predictor_sync_mats %>%
  rowwise() %>%
  mutate(tmp = list(CV_summary(Site1, Site2, Acro_df_filled,genus = "Acropora"))) %>%
  unnest(cols = c(tmp)) %>%
  mutate(PE = mean_cv/cv_sum) %>%
  ungroup()



acro_cv_results_grouped <- acro_cv_results %>%
  mutate(
    hab1 = str_extract(Site1, "Fringe|Backreef|10m|17m"),
    hab2 = str_extract(Site2, "Fringe|Backreef|10m|17m"),
    site1 = str_extract(Site1, "LTER01|LTER02|LTER03|LTER04|LTER05|LTER06"),
    site2 = str_extract(Site2, "LTER01|LTER02|LTER03|LTER04|LTER05|LTER06"),
    shore1 = case_when(site1 %in% c("LTER01", "LTER02") ~ "North",
                       site1 %in% c("LTER03", "LTER04") ~ "East",
                       site1 %in% c("LTER05", "LTER06") ~ "West"),
    shore2 = case_when(site2 %in% c("LTER01", "LTER02") ~ "North",
                       site2 %in% c("LTER03", "LTER04") ~ "East",
                       site2 %in% c("LTER05", "LTER06") ~ "West"),

    habitat_pair = case_when(
      hab1 == hab2 ~ "within",
      hab1 != hab2 ~ "between",
      TRUE ~ NA_character_
    ),
    site_pair = case_when(
      site1 == site2 ~ "within",
      site1 != site2 ~ "between",
      TRUE ~ NA_character_
    ),
    shore_pair = case_when(
      shore1 == shore2 ~ "within",
      shore1 != shore2 ~ "between",
      TRUE ~ NA_character_
    )
  )
acro_sync_base_lm <- lmer(log1p(PE) ~ ACRO_sync + (1 | Site1) + (1 | Site2), data = acro_cv_results_grouped)
summary(acro_sync_base_lm)
acro_sync_base_lm_df<-ggpredict(acro_sync_base_lm, terms = c("ACRO_sync"), back_transform = F) %>% as.data.frame()
acro_cv_results_grouped_resid <- acro_cv_results_grouped %>%
  mutate(resid_sync_only = residuals(acro_sync_base_lm))

acro_sync_base_lm_r2_vals <- r.squaredGLMM(acro_sync_base_lm)
acro_sync_base_lm_r2_vals[1]
acro_sync_base_lm_r2_vals[2]

acro_annotate_label <- paste0(
  "<span style='font-size:10pt'><i>p</i> &lt; 0.001<br>",
  "<i>R</i><sup>2<sub>m</sub></sup> = 0.01<br>",
  "<i>R</i><sup>2<sub>c</sub></sup> = 0.62</span>"
)


# habitat
summary(aov(resid_sync_only ~ habitat_pair*timescale, data = acro_cv_results_grouped_resid))

# Site pair
summary(aov(resid_sync_only ~ site_pair*timescale, data = acro_cv_results_grouped_resid))

# Shore pair
summary(aov(resid_sync_only ~ shore_pair*timescale, data = acro_cv_results_grouped_resid))

eta_squared(aov(resid_sync_only ~ habitat_pair*timescale, data = acro_cv_results_grouped_resid))
eta_squared(aov(resid_sync_only ~ site_pair*timescale, data = acro_cv_results_grouped_resid))
eta_squared(aov(resid_sync_only ~ shore_pair*timescale, data = acro_cv_results_grouped_resid))

acro_resid_long_hab<-acro_cv_results_grouped_resid %>%
  select(timescale, habitat_pair, PE, resid_sync_only) %>%
  rename(habitat = habitat_pair) %>%
  pivot_longer(cols = c("habitat"), names_to = "hierarchy", values_to = "level")

acro_resid_long_site<-acro_cv_results_grouped_resid %>%
  select(timescale, site_pair, PE, resid_sync_only) %>%
  rename(site = site_pair) %>%
  pivot_longer(cols = c("site"), names_to = "hierarchy", values_to = "level")

acro_resid_long_shore<-acro_cv_results_grouped_resid %>%
  select(timescale, shore_pair, PE, resid_sync_only) %>%
  rename(shore = shore_pair) %>%
  pivot_longer(cols = c("shore"), names_to = "hierarchy", values_to = "level")

acro_resid_long_complete<-rbind(acro_resid_long_hab,acro_resid_long_site,acro_resid_long_shore) %>%
  mutate(timescale_written = case_when(
    timescale == "short"  ~ "2-5y",
    timescale == "medium" ~ "5-10y",
    timescale == "long"   ~ "2-10y"
  ))

acro_resid_long_complete$timescale_written = factor(acro_resid_long_complete$timescale_written,
                                                    levels=c("2-5y","5-10y","2-10y"))

# Estimating pairwise differences
acro_hab_pairwise_lm<-lm(resid_sync_only ~ habitat_pair*timescale, data = acro_cv_results_grouped_resid)
acro_site_pairwise_lm<-lm(resid_sync_only ~ site_pair*timescale, data = acro_cv_results_grouped_resid)
acro_shore_pairwise_lm<-lm(resid_sync_only ~ shore_pair*timescale, data = acro_cv_results_grouped_resid)

emm_acro_hab_resid <- emmeans(acro_hab_pairwise_lm, pairwise ~ habitat_pair | timescale)
emm_acro_site_resid <- emmeans(acro_site_pairwise_lm, pairwise ~ site_pair | timescale)
emm_acro_shore_resid <- emmeans(acro_shore_pairwise_lm, pairwise ~ shore_pair | timescale)

acro_hab_resid_pairs <- contrast(emm_acro_hab_resid, method = "revpairwise") %>%
  summary(infer = TRUE) %>%
  data.frame() %>%
  mutate(hierarchy = "habitat")

acro_site_resid_pairs <- contrast(emm_acro_site_resid, method = "revpairwise") %>%
  summary(infer = TRUE) %>%
  data.frame() %>%
  mutate(hierarchy = "site")

acro_shore_resid_pairs <- contrast(emm_acro_shore_resid, method = "revpairwise") %>%
  summary(infer = TRUE) %>%
  data.frame() %>%
  mutate(hierarchy = "shore")

acro_sig_comparisons<-rbind(
  acro_hab_resid_pairs,
  acro_site_resid_pairs,
  acro_shore_resid_pairs) %>%
  mutate(sig  = p.value < 0.05) %>%
  mutate(timescale_written = case_when(
    timescale == "short"  ~ "2-5y",
    timescale == "medium" ~ "5-10y",
    timescale == "long"   ~ "2-10y"
  ))

acro_asterisks_df <- acro_resid_long_complete %>%
  group_by(hierarchy, timescale_written) %>%
  dplyr::summarise(y = mean(resid_sync_only) + 0.02, .groups = "drop") %>%
  left_join(acro_sig_comparisons %>%
              filter(sig == TRUE) %>%
              select(hierarchy, timescale_written, sig),
            by = c("hierarchy", "timescale_written")) %>%
  na.omit()

acro_asterisks_df$timescale_written = factor(acro_asterisks_df$timescale_written,
                                                    levels=c("2-5y","5-10y","2-10y"))


acro_sync_base_lm_plot<-ggplot(acro_sync_base_lm_df, aes(x = x, y = predicted)) +
  geom_point(data = acro_cv_results_grouped, aes(x = ACRO_sync, y = log1p(PE)), alpha = 0.1) +
  geom_ribbon(aes(ymin = conf.low,ymax = conf.high), alpha = 0.4, fill = "blue") +
  geom_line(color = "blue") +
  theme_bw() +
  facet_grid(~"PE~synchrony") +
  labs(x = "Wavelet synchrony", y = "ln(Portfolio effects + 1)") +
  ggtext::geom_richtext(
    aes(x = 1, y = 1.25, label = acro_annotate_label),
    fill = NA, label.color = NA,
    hjust = 1, size = 4, lineheight = 1.2
  )

acro_resid_plot<-ggplot(acro_resid_long_complete, aes(x = hierarchy, y = resid_sync_only, color = level)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_summary(
    fun.data = mean_se,
    geom = "pointrange",
    position = position_dodge(width = 0.6),
    size = 0.5
  ) +
  facet_grid("Acropora"~timescale_written) +
  theme_bw() +
  theme(strip.text.y = element_text(face = "italic")) +
  scale_y_continuous(position = "right") +
  labs(y = "PE ~ synchrony residuals", x = "Hierarchy Level") +
  scale_color_manual(values = c("within" = "black", "between" = "grey")) +
  geom_text(
    data = acro_asterisks_df,
    aes(x = hierarchy, y = y, label = "*"),
    color = "red",
    size = 6,
    inherit.aes = FALSE
  )


acro_combined_pe_lm_plot<-acro_sync_base_lm_plot + acro_resid_plot + plot_layout(widths = c(2,3))




# plot together
poc_sync_base_lm_plot_ready<-poc_sync_base_lm_plot + theme(axis.text.x = element_blank()) +
  ylim(c(0.5,1.8)) +
  annotate("text",x = -Inf, y = Inf,hjust = -0.05, vjust = 1.5, label = "(A)")

por_sync_base_lm_plot_ready<-por_sync_base_lm_plot + theme(axis.text.x = element_blank(), strip.background.x = element_blank(),strip.text.x = element_blank()) +
  ylim(c(0.5,2)) +
  annotate("text",x = -Inf, y = Inf,hjust = -0.05, vjust = 1.5, label = "(E)")

mont_sync_base_lm_plot_ready<-mont_sync_base_lm_plot + theme(axis.text.x = element_blank(), strip.background.x = element_blank(),strip.text.x = element_blank()) +
  ylim(c(0.5,2)) +
  annotate("text",x = -Inf, y = Inf,hjust = -0.05, vjust = 1.5, label = "(I)")

acro_sync_base_lm_plot_ready<-acro_sync_base_lm_plot + theme(strip.background.x = element_blank(),strip.text.x = element_blank()) +
  ylim(c(0.5,1.5)) +
  annotate("text",x = -Inf, y = Inf,hjust = -0.05, vjust = 1.5, label = "(M)")


# Labels for facet plots
poc_dat_text <- data.frame(
  label = c("(B)", "(C)", "(D)"),
  timescale_written   = c("2-5y", "5-10y", "2-10y"))

poc_dat_text$timescale_written = factor(poc_dat_text$timescale_written,
                                             levels=c("2-5y","5-10y","2-10y"))

por_dat_text <- data.frame(
  label = c("(F)", "(G)", "(H)"),
  timescale_written   = c("2-5y", "5-10y", "2-10y"))

por_dat_text$timescale_written = factor(por_dat_text$timescale_written,
                                        levels=c("2-5y","5-10y","2-10y"))

mont_dat_text <- data.frame(
  label = c("(J)", "(K)", "(L)"),
  timescale_written   = c("2-5y", "5-10y", "2-10y"))

mont_dat_text$timescale_written = factor(mont_dat_text$timescale_written,
                                        levels=c("2-5y","5-10y","2-10y"))

acro_dat_text <- data.frame(
  label = c("(N)", "(O)", "(P)"),
  timescale_written   = c("2-5y", "5-10y", "2-10y"))

acro_dat_text$timescale_written = factor(acro_dat_text$timescale_written,
                                         levels=c("2-5y","5-10y","2-10y"))

# Plot together
poc_resid_plot_ready<-poc_resid_plot +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(position = "right") +
  coord_cartesian(ylim = c(-0.08, 0.06)) +
  geom_text(
    data    = poc_dat_text,
    mapping = aes(x = -Inf, y = Inf, label = label),
    color = "black",
    hjust   = -0.05,
    vjust   = 1.5
  )

por_resid_plot_ready<-por_resid_plot +
  theme(axis.text.x = element_blank(), strip.background.x = element_blank(),strip.text.x = element_blank())+
  scale_y_continuous(position = "right") +
  coord_cartesian(ylim = c(-0.08, 0.06)) +
  geom_text(
    data    = por_dat_text,
    mapping = aes(x = -Inf, y = Inf, label = label),
    color = "black",
    hjust   = -0.05,
    vjust   = 1.5
  )
mont_resid_plot_ready<-mont_resid_plot +
  theme(axis.text.x = element_blank(), strip.background.x = element_blank(),strip.text.x = element_blank())+
  scale_y_continuous(position = "right") +
  coord_cartesian(ylim = c(-0.12, 0.06)) +
  geom_text(
    data    = mont_dat_text,
    mapping = aes(x = -Inf, y = Inf, label = label),
    color = "black",
    hjust   = -0.05,
    vjust   = 1.5
  )
acro_resid_plot_ready<-acro_resid_plot +
  scale_y_continuous(position = "right") +
  coord_cartesian(ylim = c(-0.04, 0.02)) +
  theme(strip.background.x = element_blank(),strip.text.x = element_blank()) +
  geom_text(
    data    = acro_dat_text,
    mapping = aes(x = -Inf, y = Inf, label = label),
    color = "black",
    hjust   = -0.05,
    vjust   = 1.5
  )

poc_resid_plot

pe_left <- poc_sync_base_lm_plot_ready /
  por_sync_base_lm_plot_ready /
  mont_sync_base_lm_plot_ready /
  acro_sync_base_lm_plot_ready +
  plot_layout(guides = "collect", axis_titles = "collect") &
  theme(legend.position = "bottom")

# Right column (residual effects)
sync_right <- poc_resid_plot_ready /
  por_resid_plot_ready /
  mont_resid_plot_ready /
  acro_resid_plot_ready +
  plot_layout(guides = "collect", axis_titles = "collect") &
  theme(legend.position = "bottom")

# Combine both columns
final_plot <- (pe_left | sync_right) +
  plot_layout(widths = c(1.5, 3), guides = "collect") &
  theme(legend.position = "bottom")

final_plot

ggsave('./figures/Figure_5_updline.jpeg',
       final_plot,
       dpi = 500,
       height = 7,
       width = 7)


# Show differences in PE and synchrony across groups
poc_cv_results_grouped_merge  <- poc_cv_results_grouped %>% dplyr::rename(sync = POC_sync) %>% mutate(genus = "Pocillopora")
por_cv_results_grouped_merge  <- por_cv_results_grouped %>% dplyr::rename(sync = POR_sync) %>% mutate(genus = "Porites")
mont_cv_results_grouped_merge <- mont_cv_results_grouped %>% dplyr::rename(sync = MONT_sync) %>% mutate(genus = "Montipora")
acro_cv_results_grouped_merge <- acro_cv_results_grouped %>% dplyr::rename(sync = ACRO_sync) %>% mutate(genus = "Acropora")

cv_results_all_merged<-rbind(poc_cv_results_grouped_merge,
                             por_cv_results_grouped_merge,
                             mont_cv_results_grouped_merge,
                             acro_cv_results_grouped_merge)
cv_results_all_merged$genus = factor(cv_results_all_merged$genus, levels=c("Pocillopora","Porites","Montipora","Acropora"))


write.csv(cv_results_all_merged,'./data/summarized/pe_output_data.csv', row.names = F)

sync_long <- cv_results_all_merged %>%
  select(timescale, habitat_pair, site_pair, shore_pair, sync, genus) %>%
  rename(habitat = habitat_pair,
         site = site_pair,
         shore = shore_pair) %>%
  pivot_longer(
    cols = c(habitat, site, shore),
    names_to = "hierarchy",
    values_to = "group"
  ) %>%
  mutate(timescale_written = case_when(
    timescale == "short"  ~ "2-5y",
    timescale == "medium" ~ "5-10y",
    timescale == "long"   ~ "2-10y"
  ))

sync_long$timescale_written = factor(sync_long$timescale_written,
                                         levels=c("2-5y","5-10y","2-10y"))

sync_diffs_plot<-ggplot(sync_long, aes(x = hierarchy, y = sync, color = group)) +
  geom_boxplot() +
  facet_grid(timescale_written~genus) +
  theme_bw() +
  scale_fill_manual(values = c("black","grey")) +
  scale_color_manual(values = c("black","grey")) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0,0,0,0),"cm")) +
  labs(y = "Wavelet synchrony", color = "level")

PE_long <- cv_results_all_merged %>%
  select(timescale, habitat_pair, site_pair, shore_pair, PE, genus) %>%
  rename(habitat = habitat_pair,
         site = site_pair,
         shore = shore_pair) %>%
  pivot_longer(
    cols = c(habitat, site, shore),
    names_to = "hierarchy",
    values_to = "group"
  )

pe_boxplot<-ggplot(PE_long, aes(x = hierarchy, y = log1p(PE), color = group)) +
  geom_boxplot() +
  facet_grid("full timeseries"~genus) +
  theme_bw() +
  scale_fill_manual(values = c("black","grey")) +
  scale_color_manual(values = c("black","grey")) +
  theme(
    strip.background.x = element_blank(),
    strip.text.x = element_blank(),
    plot.margin = unit(c(0,0,0,0),"cm")) +
  labs(y = "Portfolio effect (PE)", x = "spatial hierarchy", color = "level")

Figure_S6<-sync_diffs_plot / pe_boxplot  + plot_layout(guides = "collect",axes = "collect",heights = c(4,1.2)) &
  theme(legend.position = "bottom")

ggsave('./figures/Figure_S6_comparisons.jpeg',
       Figure_S6,
       dpi = 500,
       height = 8,
       width = 8)


# Table 1: results of Multiple linear regression

complete_MLR_slopes<-rbind(
  poc_MLR_slopes_all %>% mutate(genus = "Pocillopora"),
  POR_MLR_slopes_all %>% mutate(genus = "Porites"),
  MONT_MLR_slopes_all %>% mutate(genus = "Montipora"),
  ACRO_MLR_slopes_all %>% mutate(genus = "Acropora"))


complete_MLR_slopes<-complete_MLR_slopes %>%
  select(genus,Predictor,timescale,estimate,std.error,p.value) %>%
  arrange(genus,timescale,Predictor)

write.csv(complete_MLR_slopes,'./tables/MLR_tabled_data.csv',row.names = F)




# Table 2: results of Linear modeling and residual analyses
complete_LM_comparisons<-rbind(poc_sig_comparisons %>% mutate(genus = "Pocillopora"),
                           por_sig_comparisons %>% mutate(genus = "Porites"),
                           mont_sig_comparisons  %>% mutate(genus = "Montipora"),
                           acro_sig_comparisons  %>% mutate(genus = "Acropora"))

complete_LM_comparisons<-complete_LM_comparisons %>%
  select(genus,hierarchy,timescale,estimate,SE,df,p.value) %>%
  arrange(genus,timescale)

write.csv(complete_LM_comparisons,'./tables/complete_LM_comparisons.csv',row.names = F)



# Table 3: LM model results
poc_sync_base_lm_table<-data.frame(model_parameters(poc_sync_base_lm)) %>% mutate(genus = "Pocillopora")
por_sync_base_lm_table<-data.frame(model_parameters(por_sync_base_lm)) %>% mutate(genus = "Porites")
mont_sync_base_lm_table<-data.frame(model_parameters(mont_sync_base_lm)) %>% mutate(genus = "Montipora")
acro_sync_base_lm_table<-data.frame(model_parameters(acro_sync_base_lm)) %>% mutate(genus = "Acropora")


complete_lmer_table<-rbind(poc_sync_base_lm_table,
        por_sync_base_lm_table,
        mont_sync_base_lm_table,
      acro_sync_base_lm_table)

write.csv(complete_lmer_table,'./tables/complete_lmer_table.csv',row.names = F)






# Show spatial synchrony in environmental conditions
predictor_sync_long <- cv_results_all_merged %>%
  select(timescale, habitat_pair, site_pair, shore_pair, DTR_sync,DHD_sync,alg_sync) %>%
  rename(habitat = habitat_pair,
         site = site_pair,
         shore = shore_pair) %>%
  pivot_longer(
    cols = c(habitat, site, shore),
    names_to = "hierarchy",
    values_to = "group"
  ) %>%
  mutate(timescale_written = case_when(
    timescale == "short"  ~ "2-5y",
    timescale == "medium" ~ "5-10y",
    timescale == "long"   ~ "2-10y"
  ))

cv_results_all_merged %>%
  group_by(timescale,genus) %>%
  dplyr::summarize(count = n())

predictor_sync_long$timescale_written = factor(sync_long$timescale_written,
                                     levels=c("2-5y","5-10y","2-10y"))

DTR_sync_lm <- lm(DTR_sync~hierarchy*group*timescale_written, data = predictor_sync_long)
anova(DTR_sync_lm)

DTR_pred_data <- ggpredict(DTR_sync_lm, terms = c("hierarchy", "group", "timescale_written"))

DTR_predict_plot <- ggplot(DTR_pred_data, aes(x = x, y = predicted, color = group)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.6)) +
  facet_grid(facet~"DTR") +
  scale_color_manual(values = c("blue","lightblue")) +
  theme_bw() +
  labs(x = "Spatial hierarchy", y = "Pairwise wavelet synchrony") +
  theme(plot.title = element_blank(),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.78, 0.95)) +
  ylim(0.4, 1)


DHD_sync_lm <- lm(DHD_sync~hierarchy*group*timescale_written, data = predictor_sync_long)
anova(DHD_sync_lm)
DHD_pred_data <- ggpredict(DHD_sync_lm, terms = c("hierarchy", "group", "timescale_written"))


DHD_predict_plot <- ggplot(DHD_pred_data, aes(x = x, y = predicted, color = group)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.6)) +
  facet_grid(facet~"DHD") +
  scale_color_manual(values = c("red","pink")) +
  theme_bw() +
  labs(x = "Spatial hierarchy", y = "Pairwise wavelet synchrony") +
  theme(plot.title = element_blank(),
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.78, 0.75)) +
  ylim(0.4, 1)



alg_sync_lm <- lm(alg_sync~hierarchy*group*timescale_written, data = predictor_sync_long)
anova(alg_sync_lm)
alg_pred_data <- ggpredict(alg_sync_lm, terms = c("hierarchy", "group", "timescale_written"))

# Plot manually with custom dodge width
alg_predict_plot <- ggplot(alg_pred_data, aes(x = x, y = predicted, color = group)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.6)) +
  facet_grid(facet~"Macroalgae") +
  scale_color_manual(values = c("darkgreen", "green")) +
  theme_bw() +
  labs(x = "Spatial hierarchy", y = "Pairwise wavelet synchrony") +
  theme(plot.title = element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.78, 0.95)) +
  ylim(0.4, 1)



complete_predictor_model_plot<-DHD_predict_plot + DTR_predict_plot + alg_predict_plot + plot_layout(axes = "collect",axis_titles = "collect")

complete_predictor_model_plot

ggsave('./figures/Figure_S5_predictor_hierarchies.jpeg',
       complete_predictor_model_plot,
       dpi = 500,
       height = 7,
       width = 7)


# END #
