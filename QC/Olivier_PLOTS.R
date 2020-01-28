pdf("olivier_sha_plots.pdf", onefile = T)
for (i in 1:(length(kalivas_lga_measures)/2)){
  g_cohort <-lga_allsubjects_tograph %>% 
    dplyr::filter(active_lever_excel == active_lever_raw,
                  inactive_lever_excel == inactive_lever_raw,
                  infusions_excel == infusions_raw) %>% 
    dplyr::filter(resolution != "FLAG_EXPERIMENT"|is.na(resolution) ) %>% 
    mutate(session = as.numeric(session) %>% as.factor) %>% 
    ggplot(aes(x = session, group = session)) + 
    geom_boxplot(aes_string(y = kalivas_lga_measures[i])) + 
    facet_grid(~ cohort_number) +
    labs(title = paste0(gsub("_(?!excel)", " ", gsub("_excel", "", kalivas_lga_measures[i]), perl = T), "_Verified_Data_U01_Kalivas", "\n", "By cohort"),
         y = gsub("_(?!excel)", " ", gsub("_excel", "", kalivas_lga_measures[i]), perl = T)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  g_room <-lga_allsubjects_tograph %>% 
    dplyr::filter(active_lever_excel == active_lever_raw,
                  inactive_lever_excel == inactive_lever_raw,
                  infusions_excel == infusions_raw) %>% 
    mutate(session = as.numeric(session) %>% as.factor) %>%
    dplyr::filter(resolution != "FLAG_EXPERIMENT"|is.na(resolution)  ) %>% 
    ggplot(aes(x = session, group = session)) + 
    geom_boxplot(aes_string(y = kalivas_lga_measures[i])) + 
    facet_grid(~ self_administration_room) +
    labs(title = paste0(gsub("_(?!excel)", " ", gsub("_excel", "", kalivas_lga_measures[i]), perl = T), "_Verified_Data_U01_Kalivas", "\n", "By self admin room"),
         y = gsub("_(?!excel)", " ", gsub("_excel", "", kalivas_lga_measures[i]), perl = T)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  g_box <-lga_allsubjects_tograph %>% 
    dplyr::filter(active_lever_excel == active_lever_raw,
                  inactive_lever_excel == inactive_lever_raw,
                  infusions_excel == infusions_raw) %>% 
    mutate(session = as.numeric(session) %>% as.factor) %>%
    dplyr::filter(resolution != "FLAG_EXPERIMENT"|is.na(resolution)  ) %>% 
    ggplot(aes(x = session, group = session)) + 
    geom_boxplot(aes_string(y = kalivas_lga_measures[i])) + 
    facet_grid(~ self_administration_box) +
    labs(title = paste0(gsub("_(?!excel)", " ", gsub("_excel", "", kalivas_lga_measures[i]), perl = T), "_Verified_Data_U01_Kalivas", "\n", "By self admin box"),
         y = gsub("_(?!excel)", " ", gsub("_excel", "", kalivas_lga_measures[i]), perl = T)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  g_box_x <-lga_allsubjects_tograph %>% 
    dplyr::filter(active_lever_excel == active_lever_raw,
                  inactive_lever_excel == inactive_lever_raw,
                  infusions_excel == infusions_raw) %>% 
    mutate(session = as.numeric(session) %>% as.factor,
           self_administration_box = as.numeric(self_administration_box) %>% as.factor) %>%
    dplyr::filter(resolution != "FLAG_EXPERIMENT"|is.na(resolution)  ) %>% 
    ggplot(aes(x = self_administration_box, group = self_administration_box)) + 
    geom_boxplot(aes_string(y = kalivas_lga_measures[i])) + 
    # facet_grid(~ self_administration_box) +
    labs(title = paste0(gsub("_(?!excel)", " ", gsub("_excel", "", kalivas_lga_measures[i]), perl = T), "_Verified_Data_U01_Kalivas", "\n", "By self admin box"),
         y = gsub("_(?!excel)", " ", gsub("_excel", "", kalivas_lga_measures[i]), perl = T)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  g_rein <-lga_allsubjects_tograph %>% 
    dplyr::filter(active_lever_excel == active_lever_raw,
                  inactive_lever_excel == inactive_lever_raw,
                  infusions_excel == infusions_raw) %>% 
    mutate(session = as.numeric(session) %>% as.factor) %>% 
    dplyr::filter(resolution != "FLAG_EXPERIMENT"|is.na(resolution)) %>% 
    ggplot(aes(x = session, group = session, fill = reinforcement_schedule)) + 
    geom_boxplot(aes_string(y = kalivas_lga_measures[i])) + 
    # facet_grid(~ reinforcement_schedule) +
    labs(title = paste0(gsub("_(?!excel)", " ", gsub("_excel", "", kalivas_lga_measures[i]), perl = T), "_Verified_Data_U01_Kalivas", "\n", "By reinforcement schedule"),
         y = gsub("_(?!excel)", " ", gsub("_excel", "", kalivas_lga_measures[i]), perl = T)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  g_saline <-lga_allsubjects_tograph %>% 
    dplyr::filter(active_lever_excel == active_lever_raw,
                  inactive_lever_excel == inactive_lever_raw,
                  infusions_excel == infusions_raw) %>% 
    dplyr::filter(resolution != "FLAG_EXPERIMENT"|is.na(resolution) ) %>% 
    mutate(session = as.numeric(session) %>% as.factor) %>% 
    ggplot(aes(x = session, group = session)) + 
    geom_boxplot(aes_string(y = kalivas_lga_measures[i])) + 
    facet_grid(~ heroin_or_saline) +
    labs(title = paste0(gsub("_(?!excel)", " ", gsub("_excel", "", kalivas_lga_measures[i]), perl = T), "_Verified_Data_U01_Kalivas", "\n", "By heroin/saline"),
         y = gsub("_(?!excel)", " ", gsub("_excel", "", kalivas_lga_measures[i]), perl = T)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  # lga_allsubjects_tograph %>% select(cohort_number, self_administration_room) %>% table %>% kable 
  # lga_allsubjects_tograph %>% select(cohort_number, reinforcement_schedule) %>% table %>% kable 
  # lga_allsubjects_tograph %>% select(sex, self_administration_box) %>% table %>% kable 
  # lga_allsubjects_tograph %>% select(cohort_number, self_administration_box) %>% table %>% kable 
  # 
  
  print(g_cohort)
  print(g_room)
  print(g_box)
  print(g_box_x)
  print(g_rein)
  print(g_saline)
  
}

dev.off()
