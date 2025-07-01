# ============================================================
# Noise with and without annotator calibration
#
# Description:
# This script processes reef survey data collected via images and transects.
# It performs data cleaning and merging of multiple observer datasets,
# computes inter-observer agreement statistics (Fleiss’ Kappa),
# compares agreement before and after calibration,
# and visualizes agreement levels across benthic classes.
#
#
# Author: Julie Vercelloni
# Modified data: 04/06/2025
# ============================================================

#### Setup ----
rm(list = ls())
setwd(here::here())
source("R/functions.R")
source("R/packages.R")

#### Load and Prepare Pre-Calibration Data ----
dat_pre_cal <- read.csv("data/cal-precal_csv.csv") %>%
  rename(identifier_cal = identifier,
         display_name_cal = display_name) %>%
  dplyr::select(x, y, identifier_cal, display_name_cal, image_name)

dat_pre_chris <- read.csv("data/chris-precal_csv.csv") %>%
  rename(identifier_chris = identifier,
         display_name_chris = display_name) %>%
  dplyr::select(x, y, identifier_chris, display_name_chris, image_name)

dat_pre_jamie <- read.csv("data/jamie-precal_csv.csv") %>%
  rename(identifier_jamie = identifier,
         display_name_jamie = display_name) %>%
  dplyr::select(x, y, identifier_jamie, display_name_jamie, image_name)

dat_pre_all <- reduce(list(dat_pre_cal, dat_pre_chris, dat_pre_jamie), dplyr::left_join)

#### Compute Fleiss’ Kappa for Overall Pre-Calibration Agreement ----
agg_details_overall_pre <- kappam.fleiss(dat_pre_all %>%
  dplyr::select(starts_with("display")),
  detail = FALSE)

agg_details_pre <- kappam.fleiss(dat_pre_all %>%
  dplyr::select(starts_with("display")),
  detail = TRUE)

details_agg <- agg_details_pre$detail %>%
  data.frame()

#### Load and Prepare Post-Calibration Data ----
dat_raw <- read.csv("data/classification_raw.csv")

agg_details_overall <- kappam.fleiss(dat_raw %>%
  dplyr::select(starts_with("display")),
  detail = FALSE)

agg_details <- kappam.fleiss(dat_raw %>%
  dplyr::select(starts_with("display")),
  detail = TRUE)

#### Create and Prepare Summary Tables ----
post_dat <- agg_details$detail %>%
  data.frame() %>%
  pivot_wider(names_from = Var2, values_from = Freq) %>%
  filter(p.value <= 0.05) %>%
  arrange(desc(Kappa)) %>%
  mutate(Type = "post_calibration")

pre_dat <- agg_details_pre$detail %>%
  data.frame() %>%
  pivot_wider(names_from = Var2, values_from = Freq) %>%
  mutate(Type = "pre_calibration") %>%
  mutate(Kappa = -Kappa) %>%
  mutate(Kappa = case_when(
    (Kappa <= 0 & Kappa > -.003) ~ -0.1,
    TRUE ~ Kappa
  )) %>%
  filter(Var1 %in% post_dat$Var1)

#### Load Lookup Tables and Join ----
pre_lookup <- read.csv("data/pre_calib_labelset_unity.csv") %>%
  dplyr::select(-X)

post_lookup <- read.csv("data/post_calib_labelset_unity.csv") %>%
  dplyr::select(-X)

post_dat <- post_dat %>%
  inner_join(post_lookup) %>%
  mutate(Comm = str_trim(sapply(str_split(Var1, "[|]"), function(x) x[length(x)])))

pre_dat <- pre_dat %>%
  inner_join(pre_lookup) %>%
  mutate(Comm = str_trim(sapply(str_split(Var1, "[|]"), function(x) x[length(x)]))) %>%
  mutate(Plot = case_when(
    UnityLabels == Comm ~ Comm,
    TRUE ~ paste(UnityLabels, Comm, sep = "_")
  )) %>%
  mutate(across(c(Var1, Comm, Plot), factor))

post_dat <- post_dat %>%
  filter(Var1 %in% pre_dat$Var1) %>%
  droplevels() %>%
  mutate(Plot = case_when(
    UnityLabels == Comm ~ Comm,
    TRUE ~ paste(UnityLabels, Comm, sep = "_")
  )) %>%
  arrange(desc(Kappa)) %>%
  mutate(across(c(Var1, Comm, Plot), factor)) %>%
  mutate(Plot = fct_inorder(Plot))

pre_dat <- pre_dat %>%
  mutate(Plot = fct_reorder(Plot, levels(post_dat$Plot)))

dat_full <- rbind(post_dat, pre_dat) %>%
  mutate(Int = case_when(
    Kappa >= .75 & Kappa <= -0.75 ~ "excellent agreement",
    Kappa <= -0.75 ~ "excellent agreement",
    (Kappa > 0.40 & Kappa < 0.75) ~ "good agreement",
    (Kappa < -0.40 & Kappa > -0.75) ~ "good agreement",
    (Kappa <= .4 & Kappa > 0) ~ "poor agreement",
    (Kappa >= -.4 & Kappa < 0) ~ "poor agreement",
    Kappa == 0 ~ "poor agreement",
    TRUE ~ "Other"
  )) %>%
  filter(!UnityLabels == "Fishes")

#### Plot Agreement Barplot ----
p_bar <- ggplot(dat_full, aes(x = Plot, y = Kappa, group = Type, fill = Int), col = "black") +
  geom_bar(col = "black", stat = "identity", position = "identity", alpha = .4) +
  theme_bw() +
  ylim(-.8, .8) +
  ylab("Kappa score") +
  xlab("Class") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0.75, linetype = "dashed", col = "gray80") +
  geom_hline(yintercept = -0.75, linetype = "dashed", col = "gray80") +
  geom_hline(yintercept = 0.4, linetype = "dashed", col = "gray80") +
  geom_hline(yintercept = -0.4, linetype = "dashed", col = "gray80") +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 11),
    strip.background = element_rect(fill = "white"),
    panel.grid = element_blank()
  ) +
  geom_text(label = "with calibration", aes(x = 1.9, y = .80), size = 4.5) +
  geom_text(label = "without calibration", aes(x = 1.9, y = -.79), size = 4.5) +
  geom_text(label = "0.75", aes(x = 19.2, y = -.79), size = 3.7, fontface = "italic") +
  geom_text(label = "0.75", aes(x = 19.2, y = .71), size = 3.7, fontface = "italic") +
  geom_text(label = "0.40", aes(x = 19.2, y = .36), size = 3.7, fontface = "italic") +
  geom_text(label = "0.40", aes(x = 19.2, y = -.44), size = 3.7, fontface = "italic") +
  scale_fill_manual("", values = viridisLite::viridis(8)[c(3, 5, 1)])

p_bar