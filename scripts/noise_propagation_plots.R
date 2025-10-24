# ============================================================
# Propagation of noise into machine learning models
# Description:
# This script analyzes classification agreement across datasets pre- and post-calibration. 
# It generates heatmaps of F1 scores, visualizes entropy distributions by agreement counts,
# and compares percent cover estimates from individual and ensemble models,
# as well as diver-based data, for different benthic classes.
# The script also merges data from multiple sources and creates comparative visualizations.
# Author: Julie Vercelloni
# Modified data: 04/10/2025
# ============================================================

#### Setup ----
rm(list = ls())
setwd(here::here())
source("R/functions.R")
source("R/packages.R")

#### Read and Prepare F1 Score Data ----
f1_pre <- read.csv("data/f1_scores_heatmap_pre_calibration.csv") %>% mutate(type = "uncalibrated")
f1_post <- read.csv("data/f1_scores_heatmap_post_calibration.csv") %>% mutate(type = "calibrated")

f1_all <- rbind(f1_pre, f1_post) %>%
  mutate(across(c(type), as.factor))%>%
  mutate(type = relevel(type, ref = "uncalibrated")) %>%
  mutate(
    Dataset.Name   = str_replace_all(Dataset.Name, "Observer", "Annotator"),
    Dataset.Name.2 = str_replace_all(Dataset.Name.2, "Observer", "Annotator")
  )

#### Plot F1 Score Heatmap ----
p_heat <- ggplot(f1_all, aes(x = Dataset.Name, y = Dataset.Name.2, fill = F1.Score)) +
  geom_tile(color = "white") +
  geom_text(aes(
    label = round(F1.Score, 2),
    color = F1.Score > 0.99  
  ), size = 4) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") +
  scale_fill_gradient("F1", low = "#f7fbff", high = "#08306b") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.text = element_text(size = 12),
    legend.position = "none"
  ) +
  facet_wrap(~type, scales = "free") + 
  xlab("Validation") +
  ylab("Prediction")


#### Read and Prepare Entropy Data ----
ent_pre <- read.csv("data/entropies_precal.csv") %>% mutate(type = "uncalibrated")
ent_post <- read.csv("data/entropies_postcal.csv") %>% mutate(type = "calibrated")

ent_all <- rbind(ent_pre, ent_post) %>% 
  mutate(across(c(type), as.factor))%>%
  mutate(type = relevel(type, ref = "uncalibrated")) 


#### Plot Entropy Boxplots ----
p_ent <- 
ggplot(ent_all, aes(x = as.factor(Agreement.Counts), y = Entropies)) + 
  geom_boxplot(
    aes(fill = as.factor(Agreement.Counts)), 
    outlier.shape = 16,                     
    outlier.size = 2,
    alpha=.4                          
  ) +
  facet_wrap(~type) + 
  labs(
    x = "", 
    y = "Entropies",                      
  ) +
  scale_fill_viridis_d("Number of Annotators in Agreement") +  
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_blank(), 
    axis.title = element_text(size = 14),            
    axis.text = element_text(size = 12),                
    panel.grid.major = element_line(color = "gray90"),   
    panel.grid.minor = element_line(color = "gray95"),
    axis.text.x = element_blank(),
    legend.position = "bottom"
  ) 

######################### Prepare Percent Cover Data Pre and Post Calibration for Individual Models ----

cov_pre <- read.csv("data/predictions_precal.csv") 
cov_post <- read.csv("data/predictions_postcal.csv") 

cov_pre2 <- cov_pre %>%
            dplyr::select(image_name,prediction, Dataset.Name) %>%
            mutate(image = stringr::str_split(image_name, "/") %>% map_chr(., 7)) %>%
            group_by(image, Dataset.Name) %>%
            mutate(point = row_number()) %>%
            dplyr::select(., -image_name) %>%
            ungroup() %>%
            group_by(Dataset.Name) %>%
            mutate(obs_id = cur_group_id())%>%
            ungroup() %>%
            rename(class = prediction) %>%
            group_by(class) %>%
            mutate(class_id = cur_group_id()) 

cov_pre2 %>% group_by(image, Dataset.Name) %>% count()

cov_post2 <- cov_post %>% 
            dplyr::select(image_name,prediction, Dataset.Name) %>%
            mutate(image = stringr::str_split(image_name, "/") %>% map_chr(., 4)) %>%
            group_by(image, Dataset.Name) %>%
            mutate(point = row_number()) %>%
            dplyr::select(., -image_name) %>%
            ungroup() %>%
            group_by(Dataset.Name) %>%
            mutate(obs_id = cur_group_id())%>%
            ungroup() %>%
            rename(class = prediction) %>%
            group_by(class) %>%
            mutate(class_id = cur_group_id()) 

cov_post2 %>% group_by(image, Dataset.Name) %>% count()

#### Read Lookup Tables for Calibration Labels ----
pre_lookup <- read.csv("data/pre_calib_labelset_unity.csv") %>% 
  dplyr::select(-X) %>%
  rename(class = Var1)

post_lookup <- read.csv("data/post_calib_labelset_unity.csv") %>% 
  dplyr::select(-X) %>%
  rename(class = Var1)

#### Prepare Joined Tables ----
cov_pre3 <- cov_pre2 %>%
  inner_join(pre_lookup) %>%
  mutate(Comm = str_trim(sapply(str_split(class, "[|]"), function(x) x[length(x)]))) %>%
  mutate(Plot = case_when(UnityLabels == Comm ~ Comm,
                          TRUE ~ paste(UnityLabels, Comm, sep = "_"))) %>%
  mutate(across(c(Comm,Plot), factor)) 

cov_post3 <- cov_post2 %>%
  inner_join(post_lookup) %>%
  mutate(Comm = str_trim(sapply(str_split(class, "[|]"), function(x) x[length(x)]))) %>%
  mutate(Plot = case_when(UnityLabels == Comm ~ Comm,
                          TRUE ~ paste(UnityLabels, Comm, sep = "_"))) %>%
  mutate(across(c(Comm,Plot), factor)) 

#### Aggregate Percent Cover by Class for Individual Models ----
dat_pre <- cov_pre3 %>%
  filter(!Dataset.Name == "Ensemble") %>%
  group_by(Dataset.Name, UnityLabels) %>%
  summarise(COUNT = n()) %>%
  ungroup() %>%
  group_by(Dataset.Name) %>%
  mutate(TOTAL=sum(COUNT)) 


dat_post <- cov_post3 %>%
  filter(!Dataset.Name == "Ensemble") %>%
  group_by(Dataset.Name, UnityLabels) %>%
  summarise(COUNT = n()) %>%
  ungroup() %>%
  group_by(Dataset.Name) %>%
  mutate(TOTAL=sum(COUNT)) 

#### Summarize Percent Cover Statistics Pre and Post Calibration ----
dat_pre_class <- dat_pre %>% 
  group_by(UnityLabels) %>%
  summarize(mean_COUNT = mean(COUNT), sd_COUNT = sd(COUNT), mean_TOTAL = mean(TOTAL), n = n()) %>%
  mutate(mean = mean_COUNT / mean_TOTAL,
         sd = sd_COUNT / mean_TOTAL,
         SE = sd/ sqrt(n)) %>% 
  arrange(desc(mean)) %>%
  slice(1:5) %>% 
  mutate(type = "pre calibration") %>%
  mutate(source = "Individual") %>%
  dplyr::select(UnityLabels, mean, SE, type, source)


dat_post_class <- dat_post %>% 
  group_by(UnityLabels) %>%
  summarize(mean_COUNT = mean(COUNT), sd_COUNT = sd(COUNT), mean_TOTAL = mean(TOTAL), n = n()) %>%
  mutate(mean = mean_COUNT / mean_TOTAL,
         sd = sd_COUNT / mean_TOTAL,
         SE = sd/ sqrt(n)) %>% 
  arrange(desc(mean)) %>%
  slice(1:5) %>% 
  mutate(type = "post calibration") %>%
  mutate(source = "Individual") %>%
  dplyr::select(UnityLabels, mean, SE, type, source)


#### Aggregate Percent Cover for Ensemble Models ----
dat_pre_ens <- cov_pre3 %>%
  filter(Dataset.Name == "Ensemble") %>%
  group_by(Dataset.Name, UnityLabels) %>%
  summarise(COUNT = n()) %>%
  ungroup() %>%
  group_by(Dataset.Name) %>%
  mutate(TOTAL=sum(COUNT)) 


dat_post_ens <- cov_post3 %>%
  filter(Dataset.Name == "Ensemble") %>%
  group_by(Dataset.Name, UnityLabels) %>%
  summarise(COUNT = n()) %>%
  ungroup() %>%
  group_by(Dataset.Name) %>%
  mutate(TOTAL=sum(COUNT)) 

#### Summarize Ensemble Model Statistics ----
dat_pre_ens_class <- dat_pre_ens %>% 
  group_by(UnityLabels) %>%
  summarize(mean_COUNT = mean(COUNT), sd_COUNT = sd(COUNT), mean_TOTAL = mean(TOTAL), n = n()) %>%
  mutate(mean = mean_COUNT / mean_TOTAL,
         sd = sd_COUNT / mean_TOTAL,
         SE = sd/ sqrt(n)) %>% 
  arrange(desc(mean)) %>%
  slice(1:5) %>% 
  mutate(type = "pre calibration") %>%
  mutate(source = "Ensemble") %>%
  dplyr::select(UnityLabels, mean, SE, type, source) %>%
  filter(UnityLabels %in% dat_pre_class$UnityLabels)


dat_post_ens_class <- dat_post_ens %>% 
  group_by(UnityLabels) %>%
  summarize(mean_COUNT = mean(COUNT), sd_COUNT = sd(COUNT), mean_TOTAL = mean(TOTAL), n = n()) %>%
  mutate(mean = mean_COUNT / mean_TOTAL,
         sd = sd_COUNT / mean_TOTAL,
         SE = sd/ sqrt(n)) %>% 
  arrange(desc(mean)) %>%
  slice(1:5) %>% 
  mutate(type = "post calibration") %>%
  mutate(source = "Ensemble") %>%
  dplyr::select(UnityLabels, mean, SE, type, source) %>%
  filter(UnityLabels %in% dat_post_class$UnityLabels)

######################### Add Percent Cover from Diver-Based Data ----

dat_raw_wide <- read.csv("data/Ningaloo_macroalgae_10m_transect_raw_data.csv") 
lookup_calib <- read.csv("data/lookup_prepost_calib.csv")

dat_diver_pre<- dat_raw_wide %>%
  gather(key = Comm_benthic, Proportion, 3:42) %>%
  group_by(Comm_benthic)  %>%
  mutate(Proportion = Proportion / 10) %>%
  inner_join(lookup_calib) %>%
  group_by(UnityLabels, Site, Transect) %>%
  summarize(sum_prop = sum(Proportion)) %>%
  ungroup() %>%
  group_by(UnityLabels) %>%  
  summarize(mean = round(mean(sum_prop),3), sd = sd(sum_prop), n = n()) %>%
  mutate(SE = sd / sqrt(n)) %>% 
  dplyr::select(-n) %>%
  mutate(type = "pre calibration") %>%
  mutate(source = "Diver") %>%
  dplyr::select(UnityLabels, mean, SE, type, source) %>%
  filter(UnityLabels %in% dat_pre_class$UnityLabels)


dat_diver_post <- dat_raw_wide %>%
  gather(key = Comm_benthic, Proportion, 3:42) %>%
  group_by(Comm_benthic)  %>%
  mutate(Proportion = Proportion / 10) %>%
  inner_join(lookup_calib) %>%
  group_by(UnityLabels, Site, Transect) %>%
  summarize(sum_prop = sum(Proportion)) %>%
  ungroup() %>%
  group_by(UnityLabels) %>%  
  summarize(mean = round(mean(sum_prop),3), sd = sd(sum_prop), n = n()) %>%
  mutate(SE = sd / sqrt(n)) %>% 
  dplyr::select(-n) %>%
  mutate(type = "post calibration") %>%
  mutate(source = "Diver") %>%
  dplyr::select(UnityLabels, mean, SE, type, source) %>%
  filter(UnityLabels %in% dat_post_class$UnityLabels)

### Merge All Percent Cover Datasets ----

dat_class <- rbind(dat_pre_class, dat_post_class,
                    dat_pre_ens_class, dat_post_ens_class,
                    dat_diver_pre, dat_diver_post) %>% 
  mutate(across(c(type, UnityLabels), as.factor))%>%
  mutate(type = relevel(type, ref = "pre calibration")) 

dat_class <- dat_class %>%
  mutate(UnityLabels = case_when(
    UnityLabels %in% c("Sargassum", "Lobophora", "Sargassopsis") ~ paste0("*", UnityLabels, "*"),
    TRUE ~ UnityLabels
  ))

#### Plot Percent Cover Comparison Across Models and Calibration States ----
p1 <- ggplot(dat_class, aes(x = fct_reorder(UnityLabels, mean, .desc = TRUE), 
                      y = mean*100, col = source)) + 
  geom_point(size = 3, position=position_dodge(width=0.5)) +  # Slightly larger, clean color
  geom_errorbar(aes(ymin = mean*100 - 1.96 * SE*100, 
                    ymax = mean*100 + 1.96 * SE*100), 
                width = 0.2,
        position=position_dodge(width=0.5)) +  
  facet_wrap(~type, scales = "free_x") +
  xlab("Label") +
  ylab("Cover (%)") +
   theme_minimal(base_size = 14) + 
  theme(
    axis.text.x = element_markdown(angle = 45, hjust = 1),
    strip.text = element_blank(),  
    axis.title = element_text(size = 14),                
    axis.text = element_text(size = 12),                
    panel.grid.major = element_line(color = "gray90"),    
    panel.grid.minor = element_line(color = "gray95"),
    legend.position = "bottom"
  )  +
  scale_color_manual(values = c("#F1BB7B", "#FD6467", "#5B1A18") 
  )

#### Combine Plots for Presentation ----
p_fig <- p_heat / p_ent / p1 +
  plot_layout(heights = c(2, 1, 1)) +  
  plot_annotation(tag_levels = 'A')

ggsave(p_fig, file = "noise_propagation.png", width = 8, height = 10) 
