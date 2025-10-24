# ============================================================
# brms Models - Assessing level noise 
# 
# Description:
# This script processes reef survey data collected via images and transects for the model.
# Data are used into a model to distangle bias from noise sources.
# Model diagnostics are performed. 
# Author: Julie Vercelloni
# Modified data: 15/10/2025
# ============================================================

#### Setup ----

rm(list = ls())
setwd(here::here())
source("R/functions.R")
source("R/packages.R")

#### Read lookup tables ----

lookup_col <- read.csv("data/lookup_colors.csv")

lookup_table_main <- read.csv("data/lookup_table_classes_main.csv") %>%
  rename(class = image_class) %>%
  rename(Comm_benthic = transect_class)

#### Exploring transect data ----

## Raw transect data
dat_raw_wide <- read.csv("data/Ningaloo_macroalgae_10m_transect_raw_data.csv")

dat_raw_transect <- dat_raw_wide %>%
  gather(key = Comm_benthic, Proportion, 3:42) %>%
  mutate(Proportion = Proportion / 10)


## Final transect dataset
dat_transect_final <- dat_raw_transect %>%
  filter(Proportion > 0.06) %>%
  left_join(
    lookup_table_main %>%
      filter(Comm_benthic %in% unique(dat_raw_transect$Comm_benthic)) %>%
      dplyr::select(-class) %>%
      distinct()
  ) %>%
  group_by(Site, Transect, merge) %>%
  rename(Comm = merge) %>%
  mutate(Observer = paste("Diver T", Transect, sep = "")) %>%
  mutate(Cover = sum(Proportion)) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  dplyr::select(Observer, Site, Comm, Cover) %>%
  mutate(Source = "Transect")

## Add "Other" category to transect data
dat_transect_final <- dat_transect_final %>%
  group_by(Observer, Site, Source) %>%
  group_split() %>%
  map_dfr(~ add_row(
    .x, 
    Observer = first(.x$Observer), 
    Site = first(.x$Site), 
    Comm = "Other", 
    Cover = 1 - sum(.x$Cover),
    Source = first(.x$Source)
  ))

dat_transect_final %>%
  group_by(Site, Observer) %>%
  summarize(sum_cov = sum(Cover))


#### Image-based data ----

## Raw image data
dat_raw_image <- read.csv("data/classification_raw.csv")

## Reformat
dat_image <- dat_raw_image %>%
  dplyr::select(starts_with("display")) %>%
  cbind(dat_raw_image %>% dplyr::select(image_name)) %>%
  pivot_longer(!image_name, names_to = "observer", values_to = "class") %>%
  mutate(observer = stringr::str_split(observer, "_") %>% map_chr(3)) %>%
  mutate(survey   = stringr::str_split(image_name, "/") %>% map_chr(2)) %>%
  mutate(site     = stringr::str_split(image_name, "/") %>% map_chr(3)) %>%
  mutate(image    = stringr::str_split(image_name, "/") %>% map_chr(4)) %>%
  group_by(image, observer) %>%
  mutate(point = row_number()) %>%
  dplyr::select(-image_name) %>%
  ungroup() %>%
  group_by(observer) %>%
  mutate(obs_id = cur_group_id()) %>%
  ungroup() %>%
  group_by(class) %>%
  mutate(class_id = cur_group_id())


## Aggregate by site
dat_image_site <- dat_image %>%
  group_by(observer, class, site) %>%
  summarise(COUNT = n()) %>%
  ungroup() %>%
  group_by(observer, site) %>%
  mutate(TOTAL = sum(COUNT)) %>%
  ungroup() %>%
  mutate(Mean_prop = COUNT / TOTAL)


## Remove "Exclude" class and recalculate proportions
new_total <- dat_image_site %>%
  filter(str_detect(class, "Exclude")) %>%
  group_by(observer, site, TOTAL) %>%
  summarise(COUNT_rm = sum(COUNT)) %>%
  mutate(TOTAL_new = TOTAL - COUNT_rm) %>%
  dplyr::select(observer, site, TOTAL_new)

dat_image_site_reprop <- dat_image_site %>%
  inner_join(new_total) %>%
  mutate(PROP_new = round(COUNT / TOTAL_new, 3)) %>%
  filter(!str_detect(class, "Exclude"))


## Final image dataset
dat_image_final_raw <- dat_image_site_reprop %>%
  filter(PROP_new > 0.02)

dat_image_final_raw2 <- dat_image_final_raw %>%
  inner_join(
    lookup_table_main %>%
      filter(class %in% unique(dat_image_final_raw$class)) %>%
      dplyr::select(-Comm_benthic) %>%
      distinct()
  ) %>%
  ungroup() %>%
  group_by(site, observer, merge) %>%
  summarise(Cover = mean(PROP_new)) %>%
  rename(Site = site, Comm = merge, Observer = observer) %>%
  dplyr::select(Observer, Site, Comm, Cover) %>%
  mutate(Source = "Image")

dat_image_final_raw2 %>%
  group_by(Site, Observer) %>%
  summarize(sum_cov = sum(Cover))


## Add "Other" category and unify observer IDs
dat_image_final <- dat_image_final_raw2 %>%
  group_by(Observer, Site, Source) %>%
  group_split() %>%
  map_dfr(~ add_row(
    .x,
    Observer = first(.x$Observer),
    Site     = first(.x$Site),
    Comm     = "Other",
    Cover    = 1 - sum(.x$Cover),
    Source   = first(.x$Source)
  )) %>%
  mutate(Observer = case_when(
    Observer == "cal"   ~ "Observer 3",
    Observer == "chris" ~ "Observer 2",
    Observer == "jamie" ~ "Observer 1",
    Observer == "mary"  ~ "Observer 4",
    Observer == "shaun" ~ "Observer 5",
    TRUE                ~ "NA"
  ))

dat_image_final %>%
  group_by(Site, Observer) %>%
  summarize(sum_cov = sum(Cover))


#### Combine image and transect datasets ----

dat_all <- rbind(dat_image_final, dat_transect_final)

#### Model 1: Effect of source ----

# Keep classes present in both method and at more than one site 

to_keep <- dat_all %>% group_by(Comm, Source) %>% count() %>% group_by(Comm) %>% count() %>% filter(n == 2)
response_vars <- unique(to_keep$Comm) 

dat_all <- dat_all %>%
 filter(Comm %in% response_vars)

to_keep2 <- dat_all %>% group_by(Comm, Site) %>% count() %>% group_by(Comm) %>% count() %>% filter(n> 1)
response_vars <- unique(to_keep2$Comm)  

#### Filter for Comm of interest 

dat_all <- dat_all %>%
 filter(Comm %in% response_vars) %>%
  mutate(Cover = replace(Cover, Cover == 0, 0.001))  %>%
  mutate(across(c(Observer, Site, Source), as.factor)) %>%
  mutate(Observer = str_remove(Observer, " T\\d+")) %>%
  data.frame() %>%
  mutate(Source = relevel(Source, ref = "Transect"))

######### Model formula 

model_formula <- bf(Cover ~ Source + 
         (1 | Site) + (1 | Observer)) 

##################### Run model for the main groups 
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
rstan_options(threads_per_chain = 1)

plan(multisession, workers = 4) 

fits <- response_vars %>%
  set_names() %>%
  future_map(~ brm(
    formula = model_formula,
    data = dat_all %>% filter(Comm == .x),
    family = "Beta",
   iter = 1e4, warmup = 5e3, cores = 4, chains = 4,
   sample_prior = "yes", seed = 10, save_pars = save_pars(all = TRUE),
   control = list(max_treedepth = 20, adapt_delta = 0.999),
   backend = "cmdstanr"
  ))

saveRDS(fits, file = "data/model_fits.RDS")
future::plan(sequential) 

# Read model outputs - if model saved 

#fits <- readRDS("data/model_fits.RDS")

# Conditional effects

names(fits) <- response_vars

nd_a <- data.frame(
  Source = c("Transect", "Image")
)

palette = c("#FBA72A", "#D3D4D8")

conditional_plots <- imap(fits, ~ {

     title_plot <- dplyr::case_when(
     .y == "Dictyopteris" ~ "<i>Dictyopteris</i>",
     .y == "Lobophora"  ~ "<i>Lobophora</i>",
     .y == "Sargassum"  ~ "<i>Sargassum</i>",
    TRUE ~ .y 
   )


  brms::posterior_epred(.x, newdata = nd_a, re_formula = NA) |>
    as.data.frame() |>
    rename(Transect = V1, Image = V2) |>
    pivot_longer(everything(), names_to = "Source") |>
    ggplot(aes(x = value*100, fill = Source)) +
      stat_halfeye(.width = 0.75, adjust = 2, alpha = 0.6, color = NA) +
      labs(x = "Cover (%)", y = "Density", fill = "", title = title_plot) +
      scale_fill_manual(values = palette) +
      theme_pubr() +
      theme(legend.position = "top",
            legend.position.inside = c(0.8, 0.8),
            axis.title = element_text(size = 12),
            plot.title = ggtext::element_markdown(size = 12))
})

p_lay <- "
AB
CD
"
 combined_cond_plot <- (conditional_plots[[1]] + conditional_plots[[2]] + conditional_plots[[4]] + conditional_plots[[5]] ) +
  plot_layout(design = p_lay, guides = "collect", axis_titles = "collect") &
  theme(plot.tag.position  = c(.055, 1), legend.position = "bottom", legend.text = element_text(size = 14)) 
 combined_cond_plot
 
ggsave(combined_cond_plot, file = "../R/figures_paper/main/method_bias.png", width = 8, height = 8)

# Random effects

random_df_list <- imap(fits, ~ {
spread_draws(.x, r_Observer[Observer, ]) %>%
  filter(Observer != "Diver") %>% 
  compare_levels(r_Observer, by = Observer) %>%
  ungroup() %>%
  mutate(condition = reorder(Observer, r_Observer), 
        class = .y) 
})

random_df <- bind_rows(random_df_list)

ran_1 <- fits[[1]] %>%
  spread_draws(r_Observer[Observer, ]) %>%
  filter(Observer != "Diver") %>% 
  median_qi(condition_mean = r_Observer, .width = .75) %>%
  mutate(Observer = str_replace(Observer, "Observer\\.", "Annotator ")) %>%
  ggplot(aes(y = reorder(Observer, condition_mean), 
             x = condition_mean, 
             xmin = .lower, xmax = .upper)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_pointinterval(size = 1, fatten_point = 2, color = "#2C3E50") +
  labs(
    x = "Effect Size (median ± 95% CI)",
    y = "",
    title = bquote(italic(.(response_vars[1])))
  ) +
  theme_pubr() +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "none"
  )  +
    scale_x_continuous(
      breaks = seq(-0.55, 0.40, by = 0.2),
      limits = c(-0.65, 0.40)
    )

ran_2 <- fits[[2]] %>%
  spread_draws(r_Observer[Observer, ]) %>%
  filter(Observer != "Diver") %>% 
  median_qi(condition_mean = r_Observer, .width = .75) %>%
  mutate(Observer = str_replace(Observer, "Observer\\.", "Annotator ")) %>%
  ggplot(aes(y = reorder(Observer, condition_mean), 
             x = condition_mean, 
             xmin = .lower, xmax = .upper)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_pointinterval(size = 1, fatten_point = 2, color = "#2C3E50") +
  labs(
    x = "Effect Size (median ± 95% CI)",
    y = "",
    title = bquote(italic(.(response_vars[2])))
  ) +
  theme_pubr() +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "none"
  )  +
    scale_x_continuous(
      breaks = seq(-0.55, 0.40, by = 0.2),
      limits = c(-0.65, 0.40)
    )

ran_4 <- fits[[4]] %>%
  spread_draws(r_Observer[Observer, ]) %>%
  filter(Observer != "Diver") %>% 
  median_qi(condition_mean = r_Observer, .width = .75) %>%
  mutate(Observer = str_replace(Observer, "Observer\\.", "Annotator ")) %>%
  ggplot(aes(y = reorder(Observer, condition_mean), 
             x = condition_mean, 
             xmin = .lower, xmax = .upper)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_pointinterval(size = 1, fatten_point = 2, color = "#2C3E50") +
  labs(
    x = "Effect Size (median ± 95% CI)",
    y = "",
    title = bquote(italic(.(response_vars[4])))
  ) +
  theme_pubr() +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "none"
  )  +
    scale_x_continuous(
      breaks = seq(-0.55, 0.40, by = 0.2),
      limits = c(-0.65, 0.40)
    )

ran_5 <- fits[[5]] %>%
  spread_draws(r_Observer[Observer, ]) %>%
  filter(Observer != "Diver") %>% 
  median_qi(condition_mean = r_Observer, .width = .75) %>%
  mutate(Observer = str_replace(Observer, "Observer\\.", "Annotator ")) %>%
  ggplot(aes(y = reorder(Observer, condition_mean), 
             x = condition_mean, 
             xmin = .lower, xmax = .upper)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_pointinterval(size = 1, fatten_point = 2, color = "#2C3E50") +
  labs(
    x = "Effect Size (median ± 95% CI)",
    y = "",
    title = response_vars[5]
  ) +
  theme_pubr() +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "none"
  )  +
    scale_x_continuous(
      breaks = seq(-0.55, 0.40, by = 0.2),
      limits = c(-0.65, 0.40)
    )

combined_ran_plot <- (ran_1 + ran_2 + ran_4 + ran_5) +
  plot_layout(design = p_lay, guides = "collect", axis_titles = "collect") 

#ggsave(combined_ran_plot, file = "level_noise_method.png", width = 8, height = 6) 

# Model residuals

residuals <- imap(fits, ~ {
 check_brms(.x, integer = FALSE)
}
)

for (i in 1:length(response_vars)) {
  png(paste0("residuals_method_", response_vars[i], ".png"),
      width = 12, height = 12, units = "cm", res = 300)
  plotQQunif(simulationOutput = residuals[[i]], 
             testDispersion = FALSE,
             testUniformity = FALSE,
             testOutliers = FALSE,
             main = response_vars[i])
  dev.off()
}

# Model diagnostics

diagnostic_plots <- imap(fits, ~ {
  p_pp <- brms::pp_check(.x, type = "dens_overlay", ndraws = 300)
  p_check <- brms::pp_check(.x, type = "scatter_avg", ndraws = 300)
  # Combine plots 
  (p_pp | p_check) + 
    plot_annotation(
      title = .y, tag_levels = "a", tag_suffix = ')', 
      theme = theme(plot.title = element_text(hjust = 0))  # Center title
    )
})

combined_plot <- wrap_elements(diagnostic_plots[[1]]) + wrap_elements(diagnostic_plots[[2]]) +
wrap_elements(diagnostic_plots[[4]]) + wrap_elements(diagnostic_plots[[5]])

#ggsave(combined_plot, file = "model_method_diagnostics.png", width = 12, height = 10) 

#### Model 1.2: Null model ----

model_formula <- bf(Cover ~ 
         (1 | Site) + (1 | Observer)) 

##################### Run model for the main groups 
 options(mc.cores = parallel::detectCores())
 rstan_options(auto_write = TRUE)
 rstan_options(threads_per_chain = 1)

 plan(multisession, workers = 4) 

 fits_null <- response_vars %>%
   set_names() %>%
   future_map(~ brm(
     formula = model_formula,
     data = dat_all %>% filter(Comm == .x),
     family = "Beta",
    iter = 1e4, warmup = 5e3, cores = 4, chains = 4,
    sample_prior = "yes", seed = 10, save_pars = save_pars(all = TRUE),
    control = list(max_treedepth = 20, adapt_delta = 0.999),
    backend = "cmdstanr"
   ))

saveRDS(fits_null, file = "data/model_fits_null.RDS")
future::plan(sequential) 


#### Model 2: Effect of observer experience ----
dat_all <- dat_image_final

observer_field_exp <- c("Observer 2", "Observer 5", "Observer 3")

dat_all <- dat_all %>%
 mutate(exp_type = case_when(
    Observer %in% observer_field_exp ~ "field",
    TRUE ~ "none"
  )) 

# Keep classes present at all sites 
to_keep <- dat_all %>% group_by(Comm, Site) %>% count() %>% group_by(Comm) %>% count() %>% filter(n> 1)
response_vars <- unique(to_keep$Comm)  

#### Filter for Comm of interest 

dat_all <- dat_all %>%
 filter(Comm %in% response_vars) %>%
  mutate(Cover = replace(Cover, Cover == 0, 0.001))  %>%
  mutate(across(c(Observer, Site, exp_type), as.factor)) %>%
  mutate(exp_type = relevel(exp_type, ref = "field")) %>%
  droplevels() %>%
  data.frame() 

######### Model formula 

model_formula <- bf(Cover ~ exp_type + 
         (1 | Site) + (1 | Observer)) 

##################### Run model for the main groups 
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
rstan_options(threads_per_chain = 1)

plan(multisession, workers = 4) 

fits <- response_vars %>%
   set_names() %>%
   future_map(~ brm(
   formula = model_formula,
     data = dat_all %>% filter(Comm == .x),
     family = "Beta",
    iter = 1e4, warmup = 5e3, cores = 4, chains = 4,
    sample_prior = "yes", seed = 10, save_pars = save_pars(all = TRUE),
    control = list(max_treedepth = 20, adapt_delta = 0.999),
    backend = "cmdstanr"
   ))

saveRDS(fits, file = "data/model_fits_exp.RDS")
future::plan(sequential) 

# Read model outputs (if pre-saved)

#fits <- readRDS("data/model_fits_exp.RDS")

# Conditional effects

names(fits) <- response_vars

nd_a <- data.frame(
  exp_type = c("field", "none")
)

palette <- c("#E55C00", "#7D7F85")

conditional_plots <- imap(fits, ~ {

     title_plot <- dplyr::case_when(
     .y == "Dictyopteris" ~ "<i>Dictyopteris</i>",
     .y == "Lobophora"  ~ "<i>Lobophora</i>",
     .y == "Sargassum"  ~ "<i>Sargassum</i>",
    TRUE ~ .y 
   )


  brms::posterior_epred(.x, newdata = nd_a, re_formula = NA) |>
    as.data.frame() |>
    rename(`Field experience` = V1, None = V2) |>
    pivot_longer(everything(), names_to = "exp_type") |>
    ggplot(aes(x = value*100, fill = exp_type)) +
      stat_halfeye(.width = 0.75, adjust = 2, alpha = 0.6, color = NA) +
      labs(x = "Cover (%)", y = "Density", fill = "", title = title_plot) +
      scale_fill_manual(values = palette) +
      theme_pubr() +
      theme(legend.position = "top",
            legend.position.inside = c(0.8, 0.8),
            axis.title = element_text(size = 12),
            plot.title = ggtext::element_markdown(size = 12))
})

p_lay <- "
AB
CD
"
combined_cond_plot <- (conditional_plots[[1]] + conditional_plots[[2]] + conditional_plots[[4]] + conditional_plots[[5]] ) +
  plot_layout(design = p_lay, guides = "collect", axis_titles = "collect") &
  theme(plot.tag.position  = c(.055, 1), legend.position = "bottom", legend.text = element_text(size = 14)) 
 combined_cond_plot
 
ggsave(combined_cond_plot, file = "pattern_noise_field.png", width = 8, height = 8) 

# Random effects

random_df_list <- imap(fits, ~ {
spread_draws(.x, r_Observer[Observer, ]) %>%
  compare_levels(r_Observer, by = Observer) %>%
  ungroup() %>%
  mutate(condition = reorder(Observer, r_Observer), 
        class = .y) 
})

random_df <- bind_rows(random_df_list)

ran_1 <- fits[[1]] %>%
  spread_draws(r_Observer[Observer, ]) %>%
  median_qi(condition_mean = r_Observer, .width = .75) %>%
  mutate(Observer = str_replace(Observer, "Observer\\.", "Annotator ")) %>%
  ggplot(aes(y = reorder(Observer, condition_mean), 
             x = condition_mean, 
             xmin = .lower, xmax = .upper)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_pointinterval(size = 1, fatten_point = 2, color = "#2C3E50") +
  labs(
    x = "Effect Size (median ± 95% CI)",
    y = "",
    title = bquote(italic(.(response_vars[1])))
  ) +
  theme_pubr() +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "none"
  )  +
    scale_x_continuous(
      breaks = seq(-0.90, 0.70, by = 0.4),
      limits = c(-0.90, 0.70)
    )

ran_2 <- fits[[2]] %>%
  spread_draws(r_Observer[Observer, ]) %>%
  median_qi(condition_mean = r_Observer, .width = .75) %>%
  mutate(Observer = str_replace(Observer, "Observer\\.", "Annotator ")) %>%
  ggplot(aes(y = reorder(Observer, condition_mean), 
             x = condition_mean, 
             xmin = .lower, xmax = .upper)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_pointinterval(size = 1, fatten_point = 2, color = "#2C3E50") +
  labs(
    x = "Effect Size (median ± 95% CI)",
    y = "",
    title = bquote(italic(.(response_vars[2])))
  ) +
  theme_pubr() +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "none"
  )   +
    scale_x_continuous(
      breaks = seq(-0.90, 0.70, by = 0.4),
      limits = c(-0.90, 0.70)
    )

ran_4 <- fits[[4]] %>%
  spread_draws(r_Observer[Observer, ]) %>%
  median_qi(condition_mean = r_Observer, .width = .75) %>%
  mutate(Observer = str_replace(Observer, "Observer\\.", "Annotator ")) %>%
  ggplot(aes(y = reorder(Observer, condition_mean), 
             x = condition_mean, 
             xmin = .lower, xmax = .upper)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_pointinterval(size = 1, fatten_point = 2, color = "#2C3E50") +
  labs(
    x = "Effect Size (median ± 95% CI)",
    y = "",
    title = bquote(italic(.(response_vars[4])))
  ) +
  theme_pubr() +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "none"
  )   +
    scale_x_continuous(
      breaks = seq(-0.90, 0.70, by = 0.4),
      limits = c(-0.90, 0.70)
    )

ran_5 <- fits[[5]] %>%
  spread_draws(r_Observer[Observer, ]) %>%
  median_qi(condition_mean = r_Observer, .width = .75) %>%
  mutate(Observer = str_replace(Observer, "Observer\\.", "Annotator ")) %>%
  ggplot(aes(y = reorder(Observer, condition_mean), 
             x = condition_mean, 
             xmin = .lower, xmax = .upper)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_pointinterval(size = 1, fatten_point = 2, color = "#2C3E50") +
  labs(
    x = "Effect Size (median ± 95% CI)",
    y = "",
    title = response_vars[5]
  ) +
  theme_pubr() +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "none"
  )  +
    scale_x_continuous(
      breaks = seq(-0.90, 0.70, by = 0.4),
      limits = c(-0.90, 0.70)
    )

combined_ran_plot <- (ran_1 + ran_2 + ran_4 + ran_5) +
  plot_layout(design = p_lay, guides = "collect", axis_titles = "collect") 

ggsave(combined_ran_plot, file = "level_noise_exp.png", width = 8, height = 6)

# Model residuals

residuals <- imap(fits, ~ {
 check_brms(.x, integer = FALSE)
}
)

for (i in 1:length(response_vars)) {
  png(paste0("residuals_exp_", response_vars[i], ".png"), 
      width = 12, height = 12, units = "cm", res = 300)
  plotQQunif(simulationOutput = residuals[[i]], 
             testDispersion = FALSE,
             testUniformity = FALSE,
             testOutliers = FALSE,
             main = response_vars[i])
  dev.off()
}

# Model diagnostics

diagnostic_plots <- imap(fits, ~ {
  p_pp <- brms::pp_check(.x, type = "dens_overlay", ndraws = 300)
  p_check <- brms::pp_check(.x, type = "scatter_avg", ndraws = 300)
  # Combine plots 
  (p_pp | p_check) + 
    plot_annotation(
      title = .y, tag_levels = "a", tag_suffix = ')', 
      theme = theme(plot.title = element_text(hjust = 0))  # Center title
    )
})

combined_plot <- wrap_elements(diagnostic_plots[[1]]) + wrap_elements(diagnostic_plots[[2]]) +
wrap_elements(diagnostic_plots[[4]]) + wrap_elements(diagnostic_plots[[5]])

ggsave(combined_plot, file = "model_method_diagnostics_exp.png", width = 12, height = 10) 


#### Model 2.2: Null model ----

model_formula <- bf(Cover ~ 
         (1 | Site) + (1 | Observer)) 

##################### Run model for the main groups 
 options(mc.cores = parallel::detectCores())
 rstan_options(auto_write = TRUE)
 rstan_options(threads_per_chain = 1)

 plan(multisession, workers = 4) 

fits_null <- response_vars %>%
   set_names() %>%
   future_map(~ brm(
     formula = model_formula,
     data = dat_all %>% filter(Comm == .x),
     family = "Beta",
    iter = 1e4, warmup = 5e3, cores = 4, chains = 4,
    sample_prior = "yes", seed = 10, save_pars = save_pars(all = TRUE),
    control = list(max_treedepth = 20, adapt_delta = 0.999),
    backend = "cmdstanr"
   ))

saveRDS(fits_null, file = "data/model_fits_null_exp.RDS")
future::plan(sequential) 