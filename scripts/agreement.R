# ============================================================
# Noise with and without annotator calibration
#
# Description:
# This script processes reef survey data collected via images and transects.
# It performs data cleaning and merging of multiple observer datasets,
# computes inter-observer agreement statistics (Fleiss’ Kappa),
# compares agreement before and after calibration,
# and visualizes agreement levels across benthic classes.
# Author: Julie Vercelloni
# Modified data: 04/10/2025
# ============================================================

#### Setup ----
rm(list = ls())
setwd(here::here())

source("R/functions.R")
source("R/packages.R")

#### Load and Prepare Pre-Calibration Data ----
dat_pre_3 <- read.csv("data/3-precal_csv.csv") %>%
  rename(identifier_3 = identifier,
         display_name_3 = display_name) %>%
  dplyr::select(x, y, identifier_3, display_name_3, image_name)

dat_pre_2 <- read.csv("data/2-precal_csv.csv") %>%
  rename(identifier_2 = identifier,
         display_name_2 = display_name) %>%
  dplyr::select(x, y, identifier_2, display_name_2, image_name)

dat_pre_1 <- read.csv("data/1-precal_csv.csv") %>%
  rename(identifier_1 = identifier,
         display_name_1 = display_name) %>%
  dplyr::select(x, y, identifier_1, display_name_1, image_name)

dat_pre_all <- reduce(list(dat_pre_3, dat_pre_2, dat_pre_1), dplyr::left_join)

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
p_bar <- ggplot(dat_full, aes(x = Plot, y = Kappa, group = Type, fill = Int)) +
  geom_bar(col = "black", stat = "identity", position = "identity", alpha = .4) +
  theme_bw() +
  ylim(-.8, .8) +
  ylab("Kappa score") +
  xlab("Label") +
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
  geom_text(label = "calibrated", aes(x = 1.9, y = .80), size = 4.5) +
  geom_text(label = "uncalibrated", aes(x = 1.9, y = -.79), size = 4.5) +
  geom_text(label = "0.75", aes(x = 19.2, y = -.79), size = 3.7, fontface = "italic") +
  geom_text(label = "0.75", aes(x = 19.2, y = .71), size = 3.7, fontface = "italic") +
  geom_text(label = "0.40", aes(x = 19.2, y = .36), size = 3.7, fontface = "italic") +
  geom_text(label = "0.40", aes(x = 19.2, y = -.44), size = 3.7, fontface = "italic") +
  scale_fill_manual("", values = viridisLite::viridis(8)[c(3, 5, 1)]) +
  scale_x_discrete(
    labels = c(
      Substrate = "Substrate",
      Lobophora = expression(italic("Lobophora")),  
      Sargassum = expression(italic("Sargassum")),
      Turbinaria = expression(italic("Turbinaria")),
      Sargassopsis  = expression(italic("Sargassopsis ")),
      Padina = expression(italic("Padina")),
      Laurencia = expression(italic("Laurencia")),
      Dictyopteris = expression(italic("Dictyopteris")),
      Hormophysa = expression(italic("Hormophysa")),
      Dictyota = expression(italic("Dictyota")),
      `Hard coral_Corymbose` = "Hard coral_Corymbose",
      `Hard coral_Massive` = "Hard coral_Massive",
      `Cant Tell ` = "Cant Tell ",
      `Hard coral_Tabulate` = "Hard coral_Tabulate",
      `Hard coral_Encrusting` = "Hard coral_Encrusting",
      `Other_Cnidaria` = "Other_Cnidaria",
       Acanthophora = expression(italic("Acanthophora")),
      `Other_Sampling Infrastructure` = "Other_Sampling Infrastructure",
      `Hard coral_Sub-massive` = "Hard coral_Sub-massive"
    )
  )

p_bar

ggsave(p_bar, file = "agreement_prevspost.png", width = 10.5, height = 7) 
