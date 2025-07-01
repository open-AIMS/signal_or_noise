# ============================================================
# GLVVM Models - Assessing Observer Differences in Benthic Cover
# 
# Description:
# This script fits Generalized Linear Latent Variable Models (GLLVM)
# to reef survey data collected via images and transects. It:
#  - Cleans and merges image-based and transect-based observations
#  - Builds models to estimate the effect of observation method (image/transect)
#    and observer experience on benthic community structure
#  - Validates models and visualizes coefficient and fit diagnostics
#
# Author: Julie Vercelloni
# Modified data: 01/07/2025
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
dat_all$Observer <- as.factor(dat_all$Observer)
dat_all$Comm     <- as.factor(as.character(dat_all$Comm))


#### Plot community composition ----

lookup_col_temp <- lookup_col %>%
  filter(Comm %in% unique(dat_all$Comm))

fillScale <- scale_fill_manual(name = "Class", values = lookup_col_temp$colors)

ggplot(dat_all, aes(y = Cover * 100, x = Site, fill = Comm)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  fillScale +
  xlab("") +
  ylab("Cover (%)") +
  facet_wrap(~Observer, ncol = 3) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 11, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 11),
    strip.background = element_rect(fill = 'white')
  ) +
  guides(fill = guide_legend(nrow = 2))


#### Model 1: Effect of source ----

dat_clust <- dat_all %>%
  pivot_wider(names_from = Comm, values_from = Cover, values_fill = 0.001) %>%
  mutate(across(c(Observer, Site, Source), as.factor)) %>%
  mutate(Observer = str_remove(Observer, " T\\d+")) %>%
  data.frame() %>%
  mutate(Source = relevel(Source, ref = "Transect"))

y <- as.matrix(dat_clust[4:ncol(dat_clust)])
y[y == 0] <- 0.001
class(y) <- "matrix"

X <- dat_clust[1:3] %>%
  group_by(Observer) %>%
  mutate(obs_id = cur_group_id()) %>%
  ungroup() %>%
  group_by(Site) %>%
  mutate(site_id = cur_group_id()) %>%
  ungroup()

GroupStruc <- data.frame(Group = 1:nrow(y), X$site_id, X$obs_id) %>%
  rename(site_id = X.site_id, obs_id = X.obs_id)

shapeForm <- rep(1, ncol(y))

m_0 <- gllvm(y, family = "beta")

m_1 <- gllvm(y, X, family = "beta", formula = ~ Source,
             num.lv = 2, studyDesign = GroupStruc, 
             row.eff = ~(1 | site_id) + (1 | obs_id))

m_ <- m_1

p <- coef.plot(m_, "Xcoef.SourceImage")
ggsave(p, file = "appendix/figures/coefplotm1.png", width = 8, height = 6)

AIC(m_0, m_)

plot_resid <- plot_residuals(m_, X, X$Source, "residm1")

df_plot <- goodness_of_fit(m_, y)
pfit <- plot_fit(df_plot)
ggsave(pfit, file = "appendix/figures/fitm1.png", width = 8, height = 6)


#### Model 2: Effect of observer experience ----

observer_field_exp <- c("Observer 2", "Observer 5", "Observer 3")

dat_clust <- dat_all %>%
  filter(!Source == "Transect") %>%
  mutate(exp_type = case_when(
    Observer %in% observer_field_exp ~ "field",
    TRUE ~ "none"
  )) %>%
  pivot_wider(names_from = Comm, values_from = Cover, values_fill = 0.001) %>%
  mutate(across(c(Observer, Site, exp_type), as.factor)) %>%
  mutate(exp_type = relevel(exp_type, ref = "field")) %>%
  droplevels()

y <- as.matrix(dat_clust[5:ncol(dat_clust)])
y[y == 0] <- 0.001
class(y) <- "matrix"

X <- dat_clust[1:4] %>%
  group_by(Observer) %>%
  mutate(obs_id = cur_group_id()) %>%
  ungroup() %>%
  group_by(Site) %>%
  mutate(site_id = cur_group_id()) %>%
  ungroup()

GroupStruc <- data.frame(Group = 1:nrow(y), X$site_id, X$obs_id) %>%
  rename(site_id = X.site_id, obs_id = X.obs_id)

m_02 <- gllvm(y, family = "beta")

m_2 <- gllvm(y, X, family = "beta", formula = ~ exp_type,
             num.lv = 2, studyDesign = GroupStruc,
             row.eff = ~(1 | site_id) + (1 | obs_id))

m_ <- m_2

p <- coef.plot(m_, "Xcoef.exp_typenone")
ggsave(p, file = "appendix/figures/coefplotm2.png", width = 8, height = 6)

AIC(m_02, m_)

plot_resid <- plot_residuals(m_, X, X$exp_type, "residm2")

df_plot <- goodness_of_fit(m_, y)
pfit <- plot_fit(df_plot)
ggsave(pfit, file = "appendix/figures/fitm2.png", width = 8, height = 6)
