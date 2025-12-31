# Synchrony predictors script ####

# Packages
library(lme4)
library(lmerTest)
library(ggeffects)
library(emmeans)
library(MuMIn)

# Environmental variables #####
## DTR ####
DTR_clust_short<-synmat(dtr_clean$cdat,env_years,
                       method = clustmethod, nsurrogs = 1000,weighted = T,
                       tsrange = short,
                       f0 = f0_set)

DTR_clust_medium<-synmat(dtr_clean$cdat,env_years,
                        method = clustmethod, nsurrogs = 1000,weighted = T,
                        tsrange = medium,
                        f0 = f0_set)

DTR_clust_long<-synmat(dtr_clean$cdat,env_years,
                      method = clustmethod, nsurrogs = 1000,weighted = T,
                      tsrange = long,
                      f0 = f0_set)

## DHD ####
DHD_clust_short<-synmat(DHD_clean$cdat,env_years,
                       method = clustmethod, nsurrogs = 1000,weighted = T,
                       tsrange = short,
                       f0 = f0_set)

DHD_clust_medium<-synmat(DHD_clean$cdat,env_years,
                        method = clustmethod, nsurrogs = 1000,weighted = T,
                        tsrange = medium,
                        f0 = f0_set)

DHD_clust_long<-synmat(DHD_clean$cdat,env_years,
                      method = clustmethod, nsurrogs = 1000,weighted = T,
                      tsrange = long,
                      f0 = f0_set)

## Alg ####
alg_clust_short<-synmat(alg_clean$cdat,env_years,
                       method = clustmethod, nsurrogs = 1000,weighted = T,
                       tsrange = short,
                       f0 = f0_set)

alg_clust_medium<-synmat(alg_clean$cdat,env_years,
                        method = clustmethod, nsurrogs = 1000,weighted = T,
                        tsrange = medium,
                        f0 = f0_set)

alg_clust_long<-synmat(alg_clean$cdat,env_years,
                      method = clustmethod, nsurrogs = 1000,weighted = T,
                      tsrange = long,
                      f0 = f0_set)





## Flatten ####
### DTR ####
DTR_short_sync_mat <- as.matrix(DTR_clust_short) # synmat
DTR_site_names <- row.names(dtr_clean$cdat) #synmat

rownames(DTR_short_sync_mat) <- DTR_site_names
colnames(DTR_short_sync_mat) <- DTR_site_names

DTR_short_lower_idx <- which(lower.tri(DTR_short_sync_mat), arr.ind = TRUE)

DTR_clust_short_df <- data.frame(
  Site1 = DTR_site_names[DTR_short_lower_idx[, 1]],
  Site2 = DTR_site_names[DTR_short_lower_idx[, 2]],
  DTR_sync = DTR_short_sync_mat[lower.tri(DTR_short_sync_mat)],
  timescale = "short"
)


DTR_medium_sync_mat <- as.matrix(DTR_clust_medium) # synmat
DTR_site_names <- row.names(dtr_clean$cdat) #synmat

rownames(DTR_medium_sync_mat) <- DTR_site_names
colnames(DTR_medium_sync_mat) <- DTR_site_names

DTR_medium_lower_idx <- which(lower.tri(DTR_medium_sync_mat), arr.ind = TRUE)

DTR_clust_medium_df <- data.frame(
  Site1 = DTR_site_names[DTR_medium_lower_idx[, 1]],
  Site2 = DTR_site_names[DTR_medium_lower_idx[, 2]],
  DTR_sync = DTR_medium_sync_mat[lower.tri(DTR_medium_sync_mat)],
  timescale = "medium"
)


DTR_long_sync_mat <- as.matrix(DTR_clust_long) # synmat
DTR_site_names <- row.names(dtr_clean$cdat) #synmat

rownames(DTR_long_sync_mat) <- DTR_site_names
colnames(DTR_long_sync_mat) <- DTR_site_names

DTR_long_lower_idx <- which(lower.tri(DTR_long_sync_mat), arr.ind = TRUE)

DTR_clust_long_df <- data.frame(
  Site1 = DTR_site_names[DTR_long_lower_idx[, 1]],
  Site2 = DTR_site_names[DTR_long_lower_idx[, 2]],
  DTR_sync = DTR_long_sync_mat[lower.tri(DTR_long_sync_mat)],
  timescale = "long"
)


sort_sites <- function(df) {
  df %>%
    mutate(Site_min = pmin(Site1, Site2),
           Site_max = pmax(Site1, Site2)) %>%
    select(-Site1, -Site2) %>%
    rename(Site1 = Site_min, Site2 = Site_max)
}

DTR_clust_short_df <- sort_sites(DTR_clust_short_df)
DTR_clust_medium_df <- sort_sites(DTR_clust_medium_df)
DTR_clust_long_df <- sort_sites(DTR_clust_long_df)

DTR_mats<-rbind(DTR_clust_long_df,DTR_clust_medium_df,DTR_clust_short_df)

### DHD ####
DHD_short_sync_mat <- as.matrix(DHD_clust_short) # synmat
DHD_site_names <- row.names(DHD_clean$cdat) #synmat

rownames(DHD_short_sync_mat) <- DHD_site_names
colnames(DHD_short_sync_mat) <- DHD_site_names

DHD_short_lower_idx <- which(lower.tri(DHD_short_sync_mat), arr.ind = TRUE)

DHD_clust_short_df <- data.frame(
  Site1 = DHD_site_names[DHD_short_lower_idx[, 1]],
  Site2 = DHD_site_names[DHD_short_lower_idx[, 2]],
  DHD_sync = DHD_short_sync_mat[lower.tri(DHD_short_sync_mat)],
  timescale = "short"
)


DHD_medium_sync_mat <- as.matrix(DHD_clust_medium) # synmat
DHD_site_names <- row.names(DHD_clean$cdat) #synmat

rownames(DHD_medium_sync_mat) <- DHD_site_names
colnames(DHD_medium_sync_mat) <- DHD_site_names

DHD_medium_lower_idx <- which(lower.tri(DHD_medium_sync_mat), arr.ind = TRUE)

DHD_clust_medium_df <- data.frame(
  Site1 = DHD_site_names[DHD_medium_lower_idx[, 1]],
  Site2 = DHD_site_names[DHD_medium_lower_idx[, 2]],
  DHD_sync = DHD_medium_sync_mat[lower.tri(DHD_medium_sync_mat)],
  timescale = "medium"
)


DHD_long_sync_mat <- as.matrix(DHD_clust_long) # synmat
DHD_site_names <- row.names(DHD_clean$cdat) #synmat

rownames(DHD_long_sync_mat) <- DHD_site_names
colnames(DHD_long_sync_mat) <- DHD_site_names

DHD_long_lower_idx <- which(lower.tri(DHD_long_sync_mat), arr.ind = TRUE)

DHD_clust_long_df <- data.frame(
  Site1 = DHD_site_names[DHD_long_lower_idx[, 1]],
  Site2 = DHD_site_names[DHD_long_lower_idx[, 2]],
  DHD_sync = DHD_long_sync_mat[lower.tri(DHD_long_sync_mat)],
  timescale = "long"
)

DHD_clust_short_df <- sort_sites(DHD_clust_short_df)
DHD_clust_medium_df <- sort_sites(DHD_clust_medium_df)
DHD_clust_long_df <- sort_sites(DHD_clust_long_df)

DHD_mats<-rbind(DHD_clust_long_df,DHD_clust_medium_df,DHD_clust_short_df)

### alg ####
alg_short_sync_mat <- as.matrix(alg_clust_short) # synmat
alg_site_names <- row.names(alg_clean$cdat) #synmat

rownames(alg_short_sync_mat) <- alg_site_names
colnames(alg_short_sync_mat) <- alg_site_names

alg_short_lower_idx <- which(lower.tri(alg_short_sync_mat), arr.ind = TRUE)

alg_clust_short_df <- data.frame(
  Site1 = alg_site_names[alg_short_lower_idx[, 1]],
  Site2 = alg_site_names[alg_short_lower_idx[, 2]],
  alg_sync = alg_short_sync_mat[lower.tri(alg_short_sync_mat)],
  timescale = "short"
)


alg_medium_sync_mat <- as.matrix(alg_clust_medium) # synmat
alg_site_names <- row.names(alg_clean$cdat) #synmat

rownames(alg_medium_sync_mat) <- alg_site_names
colnames(alg_medium_sync_mat) <- alg_site_names

alg_medium_lower_idx <- which(lower.tri(alg_medium_sync_mat), arr.ind = TRUE)

alg_clust_medium_df <- data.frame(
  Site1 = alg_site_names[alg_medium_lower_idx[, 1]],
  Site2 = alg_site_names[alg_medium_lower_idx[, 2]],
  alg_sync = alg_medium_sync_mat[lower.tri(alg_medium_sync_mat)],
  timescale = "medium"
)


alg_long_sync_mat <- as.matrix(alg_clust_long) # synmat
alg_site_names <- row.names(alg_clean$cdat) #synmat

rownames(alg_long_sync_mat) <- alg_site_names
colnames(alg_long_sync_mat) <- alg_site_names

alg_long_lower_idx <- which(lower.tri(alg_long_sync_mat), arr.ind = TRUE)

alg_clust_long_df <- data.frame(
  Site1 = alg_site_names[alg_long_lower_idx[, 1]],
  Site2 = alg_site_names[alg_long_lower_idx[, 2]],
  alg_sync = alg_long_sync_mat[lower.tri(alg_long_sync_mat)],
  timescale = "long"
)

alg_clust_short_df <- sort_sites(alg_clust_short_df)
alg_clust_medium_df <- sort_sites(alg_clust_medium_df)
alg_clust_long_df <- sort_sites(alg_clust_long_df)

alg_mats<-rbind(alg_clust_long_df,alg_clust_medium_df,alg_clust_short_df)


# Merge all of the predictors
predictor_sync_mats<-merge(DTR_mats,DHD_mats) %>% merge(alg_mats)


# Coral data ####
### POC ####
poc_short_sync_mat <- as.matrix(poc_all_clust_short) # synmat
poc_site_names <- row.names(poc_cleaned$cdat) #synmat


rownames(poc_short_sync_mat) <- poc_site_names
colnames(poc_short_sync_mat) <- poc_site_names

poc_short_lower_idx <- which(lower.tri(poc_short_sync_mat), arr.ind = TRUE)

POC_clust_short_df <- data.frame(
  Site1 = poc_site_names[poc_short_lower_idx[, 1]],
  Site2 = poc_site_names[poc_short_lower_idx[, 2]],
  POC_sync = poc_short_sync_mat[lower.tri(poc_short_sync_mat)],
  timescale = "short"
)


poc_medium_sync_mat <- as.matrix(poc_all_clust_medium) # synmat
poc_site_names <- row.names(poc_cleaned$cdat) #synmat

rownames(poc_medium_sync_mat) <- poc_site_names
colnames(poc_medium_sync_mat) <- poc_site_names

poc_medium_lower_idx <- which(lower.tri(poc_medium_sync_mat), arr.ind = TRUE)

POC_clust_medium_df <- data.frame(
  Site1 = poc_site_names[poc_medium_lower_idx[, 1]],
  Site2 = poc_site_names[poc_medium_lower_idx[, 2]],
  POC_sync = poc_medium_sync_mat[lower.tri(poc_medium_sync_mat)],
  timescale = "medium"
)


poc_long_sync_mat <- as.matrix(poc_all_clust_long) # synmat
poc_site_names <- row.names(poc_cleaned$cdat) #synmat

rownames(poc_long_sync_mat) <- poc_site_names
colnames(poc_long_sync_mat) <- poc_site_names

poc_long_lower_idx <- which(lower.tri(poc_long_sync_mat), arr.ind = TRUE)

POC_clust_long_df <- data.frame(
  Site1 = poc_site_names[poc_long_lower_idx[, 1]],
  Site2 = poc_site_names[poc_long_lower_idx[, 2]],
  POC_sync = poc_long_sync_mat[lower.tri(poc_long_sync_mat)],
  timescale = "long"
)

POC_clust_short_df <- sort_sites(POC_clust_short_df)
POC_clust_medium_df <- sort_sites(POC_clust_medium_df)
POC_clust_long_df <- sort_sites(POC_clust_long_df)

POC_mats<-rbind(POC_clust_long_df,POC_clust_medium_df,POC_clust_short_df)


# Now merge POC with predictors
poc_predictor_sync_mats<-merge(POC_mats,predictor_sync_mats) %>%
  mutate(site_pair = paste(pmin(Site1, Site2), pmax(Site1, Site2), sep = " - "))

poc_predictor_sync_mats_anti<-anti_join(POC_mats,predictor_sync_mats)

# Then run multiple linear regression to see how predictors
poc_sync_MLR <- lm(POC_sync ~ (alg_sync + DTR_sync + DHD_sync) * timescale - timescale,
                   data = poc_predictor_sync_mats)

POC_sync_rand_MLR <- lmer(POC_sync ~ (alg_sync + DTR_sync + DHD_sync) * timescale - timescale + (1|site_pair),
                          data = poc_predictor_sync_mats)

summary(poc_sync_MLR)
summary(POC_sync_rand_MLR)
anova(poc_sync_MLR)
anova(POC_sync_rand_MLR)

poc_sync_tidy <- broom::tidy(poc_sync_MLR, conf.int = TRUE)
poc_sync_tidy

get_trends <- function(model, var, label) {
  em <- emtrends(model, ~ timescale, var = var)
  sm <- summary(em, infer = c(TRUE, TRUE)) %>% as_tibble()

  sm %>%
    transmute(timescale,
              estimate = !!sym(paste0(var, ".trend")),
              std.error = SE,
              conf.low = lower.CL,
              conf.high = upper.CL,
              p.value = p.value,
              sig = p.value < 0.05,
              Predictor = label)
}

# collect slopes for all three predictors
poc_MLR_slopes_all <- bind_rows(
  get_trends(POC_sync_rand_MLR, "alg_sync", "Macroalgae"),
  get_trends(POC_sync_rand_MLR, "DTR_sync", "DTR"),
  get_trends(POC_sync_rand_MLR, "DHD_sync", "DHD")
) %>%
  mutate(timescale = factor(timescale, levels = c("short","medium","long"))) %>%
  mutate(timescale_written = case_when(
    timescale == "short"  ~ "2-5y",
    timescale == "medium" ~ "5-10y",
    timescale == "long"   ~ "2-10y"
  ))

poc_MLR_slopes_all$timescale_written = factor(poc_MLR_slopes_all$timescale_written,
                                              levels=c("2-5y","5-10y","2-10y"))

poc_r2_val <- summary(poc_sync_MLR)$r.squared
poc_r2_lab <- paste0("R² = ", round(poc_r2_val, 2))

poc_lmer_r2_vals <- r.squaredGLMM(POC_sync_rand_MLR)

poc_r2_lab <- bquote(
  atop(
    R[m]^2 == .(round(poc_lmer_r2_vals[1], 2)),
    R[c]^2 == .(round(poc_lmer_r2_vals[2], 2))
  )
)

poc_r2_lab <- sprintf(
  "atop(R[m]^2==%.2f, R[c]^2==%.2f)",
  round(poc_lmer_r2_vals[1], 2),
  round(poc_lmer_r2_vals[2], 2)
)


poc_MLR_slopes_all$timescale_written = factor(poc_MLR_slopes_all$timescale_written,
                                              levels=c("2-5y","5-10y","2-10y"))

# Plot standardized coefficients
poc_MLR_plot<-ggplot(poc_MLR_slopes_all, aes(x = estimate, y = timescale_written, color = Predictor)) +
  geom_rect(aes(xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf),
            fill = "lightblue", alpha = 0.2, color = NA) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "mistyrose", alpha = 0.2, color = NA) +
  geom_point(aes(fill = ifelse(sig, Predictor, NA_character_), group = Predictor),
             size = 3, shape = 21, position = position_dodge(width = 0.5), na.rm = TRUE) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high,group = Predictor),
                 height = 0.01, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Macroalgae" = "darkgreen",
                                "DTR" = "blue",
                                "DHD" = "red")) +
  scale_fill_manual(values = c("Macroalgae" = "darkgreen",
                               "DTR" = "blue",
                               "DHD" = "red"),
                    na.value = "white") +
  labs(x = "Std. Effect (β) on pairwise pop. synchrony", y = "Timescale",
       title = "(A) *Pocillopora*",
       color = "Predictor") +
  guides(fill = "none") +
  theme_bw() +
  theme(plot.title = element_markdown()) +
  coord_cartesian(xlim = c(-0.4,0.7)) +
  annotate("text",
           x = -0.38,
           y = 3.3,
           label = poc_r2_lab,
           parse = T,
           hjust = 0,
           size = 3)


### POR ####
por_short_sync_mat <- as.matrix(por_all_clust_short) # synmat
por_site_names <- row.names(por_cleaned$cdat) #synmat

rownames(por_short_sync_mat) <- por_site_names
colnames(por_short_sync_mat) <- por_site_names

por_short_lower_idx <- which(lower.tri(por_short_sync_mat), arr.ind = TRUE)

POR_clust_short_df <- data.frame(
  Site1 = por_site_names[por_short_lower_idx[, 1]],
  Site2 = por_site_names[por_short_lower_idx[, 2]],
  POR_sync = por_short_sync_mat[lower.tri(por_short_sync_mat)],
  timescale = "short"
)


por_medium_sync_mat <- as.matrix(por_all_clust_medium) # synmat
por_site_names <- row.names(por_cleaned$cdat) #synmat

rownames(por_medium_sync_mat) <- por_site_names
colnames(por_medium_sync_mat) <- por_site_names

por_medium_lower_idx <- which(lower.tri(por_medium_sync_mat), arr.ind = TRUE)

POR_clust_medium_df <- data.frame(
  Site1 = por_site_names[por_medium_lower_idx[, 1]],
  Site2 = por_site_names[por_medium_lower_idx[, 2]],
  POR_sync = por_medium_sync_mat[lower.tri(por_medium_sync_mat)],
  timescale = "medium"
)


por_long_sync_mat <- as.matrix(por_all_clust_long) # synmat
por_site_names <- row.names(por_cleaned$cdat) #synmat

rownames(por_long_sync_mat) <- por_site_names
colnames(por_long_sync_mat) <- por_site_names

por_long_lower_idx <- which(lower.tri(por_long_sync_mat), arr.ind = TRUE)

POR_clust_long_df <- data.frame(
  Site1 = por_site_names[por_long_lower_idx[, 1]],
  Site2 = por_site_names[por_long_lower_idx[, 2]],
  POR_sync = por_long_sync_mat[lower.tri(por_long_sync_mat)],
  timescale = "long"
)


POR_clust_short_df <- sort_sites(POR_clust_short_df)
POR_clust_medium_df <- sort_sites(POR_clust_medium_df)
POR_clust_long_df <- sort_sites(POR_clust_long_df)

POR_mats<-rbind(POR_clust_long_df,POR_clust_medium_df,POR_clust_short_df)


# Now merge POR with predictors
POR_predictor_sync_mats<-merge(POR_mats,predictor_sync_mats) %>%
  mutate(site_pair = paste(pmin(Site1, Site2), pmax(Site1, Site2), sep = " - "))



# Then run multiple linear regression to see how predictors
POR_sync_MLR <- lm(POR_sync ~ (alg_sync + DTR_sync + DHD_sync) * timescale - timescale,
                   data = POR_predictor_sync_mats)

POR_sync_rand_MLR <- lmer(POR_sync ~ (alg_sync + DTR_sync + DHD_sync) * timescale - timescale + (1|site_pair),
                          data = POR_predictor_sync_mats)

summary(POR_sync_MLR)
anova(POR_sync_MLR)
POR_sync_tidy <- broom::tidy(POR_sync_MLR, conf.int = TRUE)
POR_sync_tidy
coefs <- coef(POR_sync_MLR)

# Clean term names and label interaction type
POR_sync_effects <- POR_sync_tidy %>%
  filter(term != "(Intercept)") %>%
  mutate(
    # interaction terms
    effect = if_else(str_detect(term, "timescalemedium"), "medium",
                     if_else(str_detect(term, "timescaleshort"),"short", "long")),

    # Clean term names
    term_clean = term %>%
      str_replace("timescaleshort", "") %>%
      str_replace("timescalemedium", "") %>%
      str_replace(":", "") %>%
      str_replace_all("_sync", "") %>%
      str_replace("alg", "Macroalgae") %>%
      str_replace("DTR", "DTR") %>%
      str_replace("DHD", "DHD"),

    # Order for plotting
    term_clean = fct_reorder(term_clean, estimate)
  )

get_trends <- function(model, var, label) {
  em <- emtrends(model, ~ timescale, var = var)
  sm <- summary(em, infer = c(TRUE, TRUE)) %>% as_tibble()

  sm %>%
    transmute(timescale,
              estimate = !!sym(paste0(var, ".trend")),
              std.error = SE,
              conf.low = lower.CL,
              conf.high = upper.CL,
              p.value = p.value,
              sig = p.value < 0.05,
              Predictor = label)
}

# collect slopes for all three predictors
POR_MLR_slopes_all <- bind_rows(
  get_trends(POR_sync_rand_MLR, "alg_sync", "Macroalgae"),
  get_trends(POR_sync_rand_MLR, "DTR_sync", "DTR"),
  get_trends(POR_sync_rand_MLR, "DHD_sync", "DHD")
) %>%
  mutate(timescale = factor(timescale, levels = c("short","medium","long"))) %>%
  mutate(timescale_written = case_when(
    timescale == "short"  ~ "2-5y",
    timescale == "medium" ~ "5-10y",
    timescale == "long"   ~ "2-10y"
  ))

POR_MLR_slopes_all$timescale_written = factor(POR_MLR_slopes_all$timescale_written,
                                              levels=c("2-5y","5-10y","2-10y"))



por_r2_val <- summary(POR_sync_MLR)$r.squared
por_r2_lab <- paste0("R² = ", round(por_r2_val, 2))

por_lmer_r2_vals <- r.squaredGLMM(POR_sync_rand_MLR)

por_r2_lab <- bquote(
  atop(
    R[m]^2 == .(round(por_lmer_r2_vals[1], 2)),
    R[c]^2 == .(round(por_lmer_r2_vals[2], 2))
  )
)

por_r2_lab <- sprintf(
  "atop(R[m]^2==%.2f, R[c]^2==%.2f)",
  round(por_lmer_r2_vals[1], 2),
  round(por_lmer_r2_vals[2], 2)
)

# Plot standardized coefficients
por_MLR_plot<-ggplot(POR_MLR_slopes_all, aes(x = estimate, y = timescale_written, color = Predictor)) +
  geom_rect(aes(xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf),
            fill = "lightblue", alpha = 0.2, color = NA) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "mistyrose", alpha = 0.2, color = NA) +
  geom_point(aes(fill = ifelse(sig, Predictor, NA_character_), group = Predictor),
             size = 3, shape = 21, position = position_dodge(width = 0.5), na.rm = TRUE) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high,group = Predictor),
                 height = 0.01, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Macroalgae" = "darkgreen",
                                "DTR" = "blue",
                                "DHD" = "red")) +
  scale_fill_manual(values = c("Macroalgae" = "darkgreen",
                               "DTR" = "blue",
                               "DHD" = "red"),
                    na.value = "white") +
  labs(x = "Std. Effect (β) on pairwise pop. synchrony", y = "Timescale",
       title = "(B) *Porites*",
       color = "Predictor") +
  guides(fill = "none") +
  theme_bw() +
  theme(plot.title = element_markdown()) +
  coord_cartesian(xlim = c(-0.4,0.7)) +
  annotate("text",
           x = -0.38,
           y = 3.3,
           label = por_r2_lab,
           parse = T,
           hjust = 0,
           size = 3)


### MONT ####
mont_short_sync_mat <- as.matrix(mont_all_clust_short) # synmat
mont_site_names <- row.names(mont_cleaned$cdat) #synmat

rownames(mont_short_sync_mat) <- mont_site_names
colnames(mont_short_sync_mat) <- mont_site_names

mont_short_lower_idx <- which(lower.tri(mont_short_sync_mat), arr.ind = TRUE)

MONT_clust_short_df <- data.frame(
  Site1 = mont_site_names[mont_short_lower_idx[, 1]],
  Site2 = mont_site_names[mont_short_lower_idx[, 2]],
  MONT_sync = mont_short_sync_mat[lower.tri(mont_short_sync_mat)],
  timescale = "short"
)

mont_medium_sync_mat <- as.matrix(mont_all_clust_medium) # synmat
mont_site_names <- row.names(mont_cleaned$cdat) #synmat

rownames(mont_medium_sync_mat) <- mont_site_names
colnames(mont_medium_sync_mat) <- mont_site_names

mont_medium_lower_idx <- which(lower.tri(mont_medium_sync_mat), arr.ind = TRUE)

MONT_clust_medium_df <- data.frame(
  Site1 = mont_site_names[mont_medium_lower_idx[, 1]],
  Site2 = mont_site_names[mont_medium_lower_idx[, 2]],
  MONT_sync = mont_medium_sync_mat[lower.tri(mont_medium_sync_mat)],
  timescale = "medium"
)


mont_long_sync_mat <- as.matrix(mont_all_clust_long) # synmat
mont_site_names <- row.names(mont_cleaned$cdat) #synmat

rownames(mont_long_sync_mat) <- mont_site_names
colnames(mont_long_sync_mat) <- mont_site_names

mont_long_lower_idx <- which(lower.tri(mont_long_sync_mat), arr.ind = TRUE)

MONT_clust_long_df <- data.frame(
  Site1 = mont_site_names[mont_long_lower_idx[, 1]],
  Site2 = mont_site_names[mont_long_lower_idx[, 2]],
  MONT_sync = mont_long_sync_mat[lower.tri(mont_long_sync_mat)],
  timescale = "long"
)

MONT_clust_short_df <- sort_sites(MONT_clust_short_df)
MONT_clust_medium_df <- sort_sites(MONT_clust_medium_df)
MONT_clust_long_df <- sort_sites(MONT_clust_long_df)

MONT_mats<-rbind(MONT_clust_long_df,MONT_clust_medium_df,MONT_clust_short_df)


# Now merge MONT with predictors
MONT_predictor_sync_mats<-merge(MONT_mats,predictor_sync_mats) %>%
  mutate(site_pair = paste(pmin(Site1, Site2), pmax(Site1, Site2), sep = " - "))



# Then run multiple linear regression to see how predictors
MONT_sync_MLR <- lm(MONT_sync ~ (alg_sync + DTR_sync + DHD_sync) * timescale - timescale,
                    data = MONT_predictor_sync_mats)

MONT_sync_rand_MLR <- lmer(MONT_sync ~ (alg_sync + DTR_sync + DHD_sync) * timescale - timescale + (1|site_pair),
                           data = MONT_predictor_sync_mats)

summary(MONT_sync_MLR)
anova(MONT_sync_MLR)
MONT_sync_tidy <- broom::tidy(MONT_sync_MLR, conf.int = TRUE)


get_trends <- function(model, var, label) {
  em <- emtrends(model, ~ timescale, var = var)
  sm <- summary(em, infer = c(TRUE, TRUE)) %>% as_tibble()

  sm %>%
    transmute(timescale,
              estimate = !!sym(paste0(var, ".trend")),
              std.error = SE,
              conf.low = lower.CL,
              conf.high = upper.CL,
              p.value = p.value,
              sig = p.value < 0.05,
              Predictor = label)
}

# collect slopes for all three predictors
MONT_MLR_slopes_all <- bind_rows(
  get_trends(MONT_sync_rand_MLR, "alg_sync", "Macroalgae"),
  get_trends(MONT_sync_rand_MLR, "DTR_sync", "DTR"),
  get_trends(MONT_sync_rand_MLR, "DHD_sync", "DHD")
) %>%
  mutate(timescale = factor(timescale, levels = c("short","medium","long"))) %>%
  mutate(timescale_written = case_when(
    timescale == "short"  ~ "2-5y",
    timescale == "medium" ~ "5-10y",
    timescale == "long"   ~ "2-10y"
  ))

MONT_MLR_slopes_all$timescale_written = factor(MONT_MLR_slopes_all$timescale_written,
                                               levels=c("2-5y","5-10y","2-10y"))


mont_r2_val <- summary(MONT_sync_MLR)$r.squared
mont_r2_lab <- paste0("R² = ", round(mont_r2_val, 2))

mont_lmer_r2_vals <- r.squaredGLMM(MONT_sync_rand_MLR)

mont_r2_lab <- bquote(
  atop(
    R[m]^2 == .(round(mont_lmer_r2_vals[1], 2)),
    R[c]^2 == .(round(mont_lmer_r2_vals[2], 2))
  )
)

mont_r2_lab <- sprintf(
  "atop(R[m]^2==%.2f, R[c]^2==%.2f)",
  round(mont_lmer_r2_vals[1], 2),
  round(mont_lmer_r2_vals[2], 2)
)

# Plot standardized coefficients
mont_MLR_plot<-ggplot(MONT_MLR_slopes_all, aes(x = estimate, y = timescale_written, color = Predictor)) +
  geom_rect(aes(xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf),
            fill = "lightblue", alpha = 0.2, color = NA) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "mistyrose", alpha = 0.2, color = NA) +
  annotate("text", x = -0.23, y = 0.5, label = "Desynchronized",
           hjust = 0.5, size = 3.5, fontface = "italic") +
  annotate("text", x = 0.4, y = 0.5, label = "Moran effect",
           hjust = 0.5, size = 3.5, fontface = "italic") +
  geom_point(aes(fill = ifelse(sig, Predictor, NA_character_), group = Predictor),
             size = 3, shape = 21, position = position_dodge(width = 0.5), na.rm = TRUE) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high,group = Predictor),
                 height = 0.01, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Macroalgae" = "darkgreen",
                                "DTR" = "blue",
                                "DHD" = "red")) +
  scale_fill_manual(values = c("Macroalgae" = "darkgreen",
                               "DTR" = "blue",
                               "DHD" = "red"),
                    na.value = "white") +
  labs(x = "Std. Effect (β) on pairwise pop. synchrony", y = "Timescale",
       title = "(C) *Montipora*",
       color = "Predictor") +
  guides(fill = "none") +
  theme_bw() +
  theme(plot.title = element_markdown()) +
  coord_cartesian(xlim = c(-0.4,0.7)) +
  annotate("text",
           x = -0.38,
           y = 3.3,
           label = mont_r2_lab,
           hjust = 0,
           parse = T,
           size = 3)


### ACRO ####
acro_short_sync_mat <- as.matrix(acro_all_clust_short) # synmat
acro_site_names <- row.names(acro_df_short_clean$cdat) #synmat

rownames(acro_short_sync_mat) <- acro_site_names
colnames(acro_short_sync_mat) <- acro_site_names

acro_short_lower_idx <- which(lower.tri(acro_short_sync_mat), arr.ind = TRUE)

ACRO_clust_short_df <- data.frame(
  Site1 = acro_site_names[acro_short_lower_idx[, 1]],
  Site2 = acro_site_names[acro_short_lower_idx[, 2]],
  ACRO_sync = acro_short_sync_mat[lower.tri(acro_short_sync_mat)],
  timescale = "short"
)


acro_medium_sync_mat <- as.matrix(acro_all_clust_medium) # synmat
acro_site_names <- row.names(acro_df_short_clean$cdat) #synmat

rownames(acro_medium_sync_mat) <- acro_site_names
colnames(acro_medium_sync_mat) <- acro_site_names

acro_medium_lower_idx <- which(lower.tri(acro_medium_sync_mat), arr.ind = TRUE)

ACRO_clust_medium_df <- data.frame(
  Site1 = acro_site_names[acro_medium_lower_idx[, 1]],
  Site2 = acro_site_names[acro_medium_lower_idx[, 2]],
  ACRO_sync = acro_medium_sync_mat[lower.tri(acro_medium_sync_mat)],
  timescale = "medium"
)

acro_long_sync_mat <- as.matrix(acro_all_clust_long) # synmat
acro_site_names <- row.names(acro_df_short_clean$cdat) #synmat

rownames(acro_long_sync_mat) <- acro_site_names
colnames(acro_long_sync_mat) <- acro_site_names

acro_long_lower_idx <- which(lower.tri(acro_long_sync_mat), arr.ind = TRUE)

ACRO_clust_long_df <- data.frame(
  Site1 = acro_site_names[acro_long_lower_idx[, 1]],
  Site2 = acro_site_names[acro_long_lower_idx[, 2]],
  ACRO_sync = acro_long_sync_mat[lower.tri(acro_long_sync_mat)],
  timescale = "long"
)



ACRO_clust_short_df <- sort_sites(ACRO_clust_short_df)
ACRO_clust_medium_df <- sort_sites(ACRO_clust_medium_df)
ACRO_clust_long_df <- sort_sites(ACRO_clust_long_df)

ACRO_mats<-rbind(ACRO_clust_long_df,ACRO_clust_medium_df,ACRO_clust_short_df)


# Now merge ACRO with predictors
ACRO_predictor_sync_mats<-merge(ACRO_mats,predictor_sync_mats) %>%
  mutate(site_pair = paste(pmin(Site1, Site2), pmax(Site1, Site2), sep = " - "))



# Then run multiple linear regression to see how predictors
ACRO_sync_MLR <- lm(ACRO_sync ~ (alg_sync + DTR_sync + DHD_sync) * timescale - timescale,
                    data = ACRO_predictor_sync_mats)

ACRO_sync_rand_MLR <- lmer(ACRO_sync ~ (alg_sync + DTR_sync + DHD_sync) * timescale - timescale + (1|site_pair),
                           data = ACRO_predictor_sync_mats)


ACRO_sync_MLR
summary(ACRO_sync_MLR)
anova(ACRO_sync_MLR)
ACRO_sync_tidy <- broom::tidy(ACRO_sync_MLR, conf.int = TRUE)
ACRO_sync_tidy


get_trends <- function(model, var, label) {
  em <- emtrends(model, ~ timescale, var = var)
  sm <- summary(em, infer = c(TRUE, TRUE)) %>% as_tibble()

  sm %>%
    transmute(timescale,
              estimate = !!sym(paste0(var, ".trend")),
              std.error = SE,
              conf.low = lower.CL,
              conf.high = upper.CL,
              p.value = p.value,
              sig = p.value < 0.05,
              Predictor = label)
}

# collect slopes for all three predictors
ACRO_MLR_slopes_all <- bind_rows(
  get_trends(ACRO_sync_rand_MLR, "alg_sync", "Macroalgae"),
  get_trends(ACRO_sync_rand_MLR, "DTR_sync", "DTR"),
  get_trends(ACRO_sync_rand_MLR, "DHD_sync", "DHD")
) %>%
  mutate(timescale = factor(timescale, levels = c("short","medium","long"))) %>%
  mutate(timescale_written = case_when(
    timescale == "short"  ~ "2-5y",
    timescale == "medium" ~ "5-10y",
    timescale == "long"   ~ "2-10y"
  ))

ACRO_MLR_slopes_all$timescale_written = factor(ACRO_MLR_slopes_all$timescale_written,
                                               levels=c("2-5y","5-10y","2-10y"))





acro_r2_val <- summary(ACRO_sync_MLR)$r.squared
acro_r2_lab <- paste0("r² = ", round(acro_r2_val, 2))

acro_lmer_r2_vals <- r.squaredGLMM(ACRO_sync_rand_MLR)

acro_r2_lab <- bquote(
  atop(
    R[m]^2 == .(round(acro_lmer_r2_vals[1], 2)),
    R[c]^2 == .(round(acro_lmer_r2_vals[2], 2))
  )
)
acro_r2_lab <- sprintf(
  "atop(R[m]^2==%.2f, R[c]^2==%.2f)",
  round(acro_lmer_r2_vals[1], 2),
  round(acro_lmer_r2_vals[2], 2)
)

# Plot standardized coefficients
acro_MLR_plot<-ggplot(ACRO_MLR_slopes_all, aes(x = estimate, y = timescale_written, color = Predictor)) +
  geom_rect(aes(xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf),
            fill = "lightblue", alpha = 0.2, color = NA) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "mistyrose", alpha = 0.2, color = NA) +
  annotate("text", x = -0.23, y = 0.5, label = "Desynchronized",
           hjust = 0.5, size = 3.5, fontface = "italic") +
  annotate("text", x = 0.4, y = 0.5, label = "Moran effect",
           hjust = 0.5, size = 3.5, fontface = "italic") +
  geom_point(aes(fill = ifelse(sig, Predictor, NA_character_), group = Predictor),
             size = 3, shape = 21, position = position_dodge(width = 0.5), na.rm = TRUE) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high,group = Predictor),
                 height = 0.01, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Macroalgae" = "darkgreen",
                                "DTR" = "blue",
                                "DHD" = "red")) +
  scale_fill_manual(values = c("Macroalgae" = "darkgreen",
                               "DTR" = "blue",
                               "DHD" = "red"),
                    na.value = "white") +
  labs(x = "Std. Effect (β) on pairwise pop. synchrony", y = "Timescale",
       title = "(D) *Acropora*",
       color = "Predictor") +
  guides(fill = "none") +
  theme_bw() +
  coord_cartesian(xlim = c(-0.4,0.7)) +
  theme(plot.title = element_markdown()) +
  annotate("text",
           x = -0.38,
           y = 3.3,
           label = acro_r2_lab,
           parse = T,
           hjust = 0,
           size = 3)


combined_sync_MLR<-poc_MLR_plot + por_MLR_plot + mont_MLR_plot + acro_MLR_plot + plot_layout(guides = "collect",axes = "collect",axis_titles = "collect")



ggsave('./figures/Figure_4.jpeg',
       combined_sync_MLR,
       dpi = 500,
       height = 7,
       width = 8)

# END #
