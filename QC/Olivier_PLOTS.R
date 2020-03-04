library(knitr)
library(ggforce)

## QC PLOTS FOR DATA VERIFIED BY EXCEL AND FOR RAW

## see objects from the raw vs excel qc 
### long access 

## uniform variables: session_length_hours, reinforcer, bolus_volume_ml, dose_ug_kg_infusion, time_out_seconds
setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/QC")




################################
########## SHA #################
################################

pdf("olivier_shaverified.pdf", onefile = T)
for (i in 1:(length(olivier_sha_measures))){
  g_cohort <- rewards_sha_tograph %>% 
    dplyr::filter(rewards_raw == rewards_excel) %>% 
    # dplyr::filter(resolution != "FLAG_EXPERIMENT"|is.na(resolution) ) %>% 
    mutate(exp = as.factor(exp)) %>% 
    ggplot(aes(x = exp, group = exp)) + 
    geom_boxplot(aes(y = rewards_raw), outlier.size=0.1) + 
    ggforce::facet_wrap_paginate(~ cohort, ncol = 2, nrow = 2, page = i, strip.position="top") +
    # facet_grid(~ cohort) +
    labs(title = paste0("SHA_Rewards", "_Verified_Data_U01_Olivier", "\n", "By cohort"),
         y = "SHA_Rewards") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) 
  
  g_sex <- rewards_sha_tograph %>%
    dplyr::filter(rewards_raw == rewards_excel) %>%
    # dplyr::filter(resolution != "FLAG_EXPERIMENT"|is.na(resolution) ) %>%
    mutate(exp = as.factor(exp)) %>%
    ggplot(aes(x = exp, y = rewards_raw)) +
    geom_boxplot(aes(fill = sex),outlier.size=0.1) +
    ggforce::facet_wrap_paginate(~ cohort, ncol = 2, nrow = 2, page = i, strip.position="top") +
    # facet_grid(~ cohort) +
    labs(title = paste0("SHA_Rewards", "_Verified_Data_U01_Olivier", "\n", "By cohort and sex"),
         y = "SHA_Rewards") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
  
  g_dir <- rewards_sha_tograph %>%
    dplyr::filter(rewards_raw == rewards_excel) %>%
    # dplyr::filter(resolution != "FLAG_EXPERIMENT"|is.na(resolution) ) %>%
    mutate(exp = as.factor(exp)) %>%
    ggplot(aes(x = exp, y = rewards_raw)) +
    geom_boxplot(aes(fill = directory),outlier.size=0.1) +
    ggforce::facet_wrap_paginate(~ cohort, ncol = 2, nrow = 2, page = i, strip.position="top", scales="free_y") +
    # facet_grid(~ cohort) +
    labs(title = paste0("SHA_Rewards", "_Verified_Data_U01_Olivier", "\n", "By cohort and directory"),
         y = "SHA_Rewards") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
  
  g_dir_sex <- rewards_sha_tograph %>%
    dplyr::filter(rewards_raw == rewards_excel) %>%
    mutate(dir_sex = paste0(directory, "_", sex)) %>% 
    # dplyr::filter(resolution != "FLAG_EXPERIMENT"|is.na(resolution) ) %>%
    mutate(exp = as.factor(exp)) %>%
    ggplot(aes(x = exp, y = rewards_raw, group = dir_sex)) +
    geom_point(aes(fill = dir_sex), stat='summary', fun.y=median) + 
    stat_summary(fun.y=median, geom="line", aes(color = dir_sex)) +
    # ggforce::facet_wrap_paginate(~ box, ncol = 2, nrow = 2, page = i, strip.position="top", scales="free_y") +
    # facet_grid(~ cohort) +
    labs(title = paste0("SHA_Rewards", "_Verified_Data_U01_Olivier", "\n", "By directory and sex"),
         y = "SHA_Rewards") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
  
  g_age <- rewards_sha_tograph %>% 
    dplyr::filter(rewards_raw == rewards_excel) %>% 
    # dplyr::filter(resolution != "FLAG_EXPERIMENT"|is.na(resolution) ) %>% 
    mutate(exp = as.factor(exp)) %>% 
    ggplot(aes(x = exp, group = exp)) + 
    geom_boxplot(aes(y = exp_age), outlier.size=0.1) + 
    ggforce::facet_wrap_paginate(~ cohort, ncol = 2, nrow = 2, page = i, strip.position="top", scales="free_y") +
    # facet_grid(~ cohort) +
    labs(title = paste0("SHA_Age", "_Verified_Data_U01_Olivier", "\n", "By cohort"),
         y = "SHA_Age") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) 
  
  g_age_cohort <- rewards_sha_tograph %>%
    dplyr::filter(rewards_raw == rewards_excel) %>%
    # dplyr::filter(resolution != "FLAG_EXPERIMENT"|is.na(resolution) ) %>%
    mutate(exp = as.factor(exp)) %>%
    ggplot(aes(x = exp, y = exp_age, group = cohort)) +
    geom_point(aes(fill = cohort), stat='summary', fun.y=median) + 
    stat_summary(fun.y=median, geom="line", aes(color = cohort)) +
    # ggforce::facet_wrap_paginate(~ box, ncol = 2, nrow = 2, page = i, strip.position="top", scales="free_y") +
    # facet_grid(~ cohort) +
    labs(title = paste0("SHA_Age_Median", "_Verified_Data_U01_Olivier", "\n", "By cohort"),
         y = "SHA_Age") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
  
  
  g_age_ind <- rewards_sha_tograph %>% 
    dplyr::filter(rewards_raw == rewards_excel) %>% 
    # dplyr::filter(resolution != "FLAG_EXPERIMENT"|is.na(resolution) ) %>% 
    mutate(exp = as.factor(exp)) %>% 
    ggplot(aes(x = exp, group = labanimalid, color = cohort)) + 
    geom_line(aes(y = exp_age)) + 
    geom_point(aes(x = exp, y = exp_age)) +
    ggforce::facet_wrap_paginate(~ cohort, ncol = 2, nrow = 2, page = i, strip.position="top", scales="free_y") +
    # facet_grid(~ cohort) +
    labs(title = paste0("SHA_Age", "_Verified_Data_U01_Olivier", "\n", "By cohort"),
         y = "SHA_Age") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) 
  
  # g_rein <-rewards_sha_tograph %>% 
  #   dplyr::filter(active_lever_excel == active_lever_raw,
  #                 inactive_lever_excel == inactive_lever_raw,
  #                 infusions_excel == infusions_raw) %>% 
  #   mutate(session = as.numeric(session) %>% as.factor) %>% 
  #   dplyr::filter(resolution != "FLAG_EXPERIMENT"|is.na(resolution)) %>% 
  #   ggplot(aes(x = session, group = session, fill = reinforcement_schedule)) + 
  #   geom_boxplot(aes_string(y = olivier_sha_measures[i])) + 
  #   # facet_grid(~ reinforcement_schedule) +
  #   labs(title = paste0(gsub("_(?!excel)", " ", gsub("_excel", "", olivier_sha_measures[i]), perl = T), "_Verified_Data_U01_olivier", "\n", "By reinforcement schedule"),
  #        y = gsub("_(?!excel)", " ", gsub("_excel", "", olivier_sha_measures[i]), perl = T)) + 
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  # 
  # g_saline <-rewards_sha_tograph %>% 
  #   dplyr::filter(active_lever_excel == active_lever_raw,
  #                 inactive_lever_excel == inactive_lever_raw,
  #                 infusions_excel == infusions_raw) %>% 
  #   dplyr::filter(resolution != "FLAG_EXPERIMENT"|is.na(resolution) ) %>% 
  #   mutate(session = as.numeric(session) %>% as.factor) %>% 
  #   ggplot(aes(x = session, group = session)) + 
  #   geom_boxplot(aes_string(y = olivier_sha_measures[i])) + 
  #   facet_grid(~ heroin_or_saline) +
  #   labs(title = paste0(gsub("_(?!excel)", " ", gsub("_excel", "", olivier_sha_measures[i]), perl = T), "_Verified_Data_U01_olivier", "\n", "By heroin/saline"),
  #        y = gsub("_(?!excel)", " ", gsub("_excel", "", olivier_sha_measures[i]), perl = T)) + 
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  # 
  # # rewards_sha_tograph %>% select(cohort, self_administration_room) %>% table %>% kable 
  # # rewards_sha_tograph %>% select(cohort, reinforcement_schedule) %>% table %>% kable 
  # # rewards_sha_tograph %>% select(sex, self_administration_box) %>% table %>% kable 
  # # rewards_sha_tograph %>% select(cohort, self_administration_box) %>% table %>% kable 
  # # 
  
  print(g_cohort)
  print(g_sex)
  print(g_dir)
  print(g_dir_sex)
  print(g_age)
  print(g_age_cohort)
  print(g_age_ind)
  # print(g_box)
  # print(g_box_x)
  # print(g_rein)
  # print(g_saline)
  # 
}

dev.off()










################################
########## LGA #################
################################

pdf("olivier_lgaverified.pdf", onefile = T)
for (i in 1:(length(olivier_lga_measures))){
  g_cohort <- rewards_lga_tograph %>% 
    dplyr::filter(rewards_raw == rewards_excel) %>% 
    # dplyr::filter(resolution != "FLAG_EXPERIMENT"|is.na(resolution) ) %>% 
    mutate(exp = as.factor(exp)) %>% 
    ggplot(aes(x = exp, group = exp)) + 
    geom_boxplot(aes(y = rewards_raw), outlier.size=0.1) + 
    ggforce::facet_wrap_paginate(~ cohort, ncol = 2, nrow = 2, page = i, strip.position="top", scales="free_y") +
    # facet_grid(~ cohort) +
    labs(title = paste0("LGA_Rewards", "_Verified_Data_U01_Olivier", "\n", "By cohort"),
         y = "LGA_Rewards") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) 

  g_sex <- rewards_lga_tograph %>%
    dplyr::filter(rewards_raw == rewards_excel) %>%
    # dplyr::filter(resolution != "FLAG_EXPERIMENT"|is.na(resolution) ) %>%
    mutate(exp = as.factor(exp)) %>%
    ggplot(aes(x = exp, y = rewards_raw)) +
    geom_boxplot(aes(fill = sex),outlier.size=0.1) +
    ggforce::facet_wrap_paginate(~ cohort, ncol = 2, nrow = 2, page = i, strip.position="top") +
    # facet_grid(~ cohort) +
    labs(title = paste0("LGA_Rewards", "_Verified_Data_U01_Olivier", "\n", "By cohort and sex"),
         y = "LGA_Rewards") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
  
  
  g_sex_upto14 <- rewards_lga_tograph  %>% subset(!grepl("1[5-9]|2[0-9]", exp)) %>%
    dplyr::filter(rewards_raw == rewards_excel) %>%
    # dplyr::filter(resolution != "FLAG_EXPERIMENT"|is.na(resolution) ) %>%
    mutate(exp = as.factor(exp)) %>%
    ggplot(aes(x = exp, y = rewards_raw)) +
    geom_boxplot(aes(fill = sex),outlier.size=0.1) +
    ggforce::facet_wrap_paginate(~ cohort, ncol = 2, nrow = 2, page = i, strip.position="top") +
    # facet_grid(~ cohort) +
    labs(title = paste0("LGA_Rewards (Up to 14 sessions)", "_Verified_Data_U01_Olivier", "\n", "By cohort and sex"),
         y = "LGA_Rewards") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
  
  
  
  g_dir <- rewards_lga_tograph %>%
    dplyr::filter(rewards_raw == rewards_excel) %>%
    # dplyr::filter(resolution != "FLAG_EXPERIMENT"|is.na(resolution) ) %>%
    mutate(exp = as.factor(exp)) %>%
    ggplot(aes(x = exp, y = rewards_raw)) +
    geom_boxplot(aes(fill = directory),outlier.size=0.1) +
    ggforce::facet_wrap_paginate(~ cohort, ncol = 2, nrow = 2, page = i, strip.position="top", scales="free_y") +
    # facet_grid(~ cohort) +
    labs(title = paste0("LGA_Rewards", "_Verified_Data_U01_Olivier", "\n", "By cohort and directory"),
         y = "LGA_Rewards") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

  g_dir_sex <- rewards_lga_tograph %>%
    dplyr::filter(rewards_raw == rewards_excel) %>%
    mutate(dir_sex = paste0(directory, "_", sex)) %>% 
    # dplyr::filter(resolution != "FLAG_EXPERIMENT"|is.na(resolution) ) %>%
    mutate(exp = as.factor(exp)) %>%
    ggplot(aes(x = exp, y = rewards_raw, group = dir_sex)) +
    geom_point(aes(fill = dir_sex), stat='summary', fun.y=median) + 
    stat_summary(fun.y=median, geom="line", aes(color = dir_sex)) +
    # ggforce::facet_wrap_paginate(~ box, ncol = 2, nrow = 2, page = i, strip.position="top", scales="free_y") +
    # facet_grid(~ cohort) +
    labs(title = paste0("LGA_Rewards", "_Verified_Data_U01_Olivier", "\n", "By directory and sex"),
         y = "LGA_Rewards") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

  g_age <- rewards_lga_tograph %>% 
    dplyr::filter(rewards_raw == rewards_excel) %>% 
    # dplyr::filter(resolution != "FLAG_EXPERIMENT"|is.na(resolution) ) %>% 
    mutate(exp = as.factor(exp)) %>% 
    ggplot(aes(x = exp, group = exp)) + 
    geom_boxplot(aes(y = exp_age), outlier.size=0.1) + 
    ggforce::facet_wrap_paginate(~ cohort, ncol = 2, nrow = 2, page = i, strip.position="top", scales="free_y") +
    # facet_grid(~ cohort) +
    labs(title = paste0("LGA_Age", "_Verified_Data_U01_Olivier", "\n", "By cohort"),
         y = "LGA_Age") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) 
  
  g_age_cohort <- rewards_lga_tograph %>%
    dplyr::filter(rewards_raw == rewards_excel) %>%
    # dplyr::filter(resolution != "FLAG_EXPERIMENT"|is.na(resolution) ) %>%
    mutate(exp = as.factor(exp)) %>%
    ggplot(aes(x = exp, y = exp_age, group = cohort)) +
    geom_point(aes(fill = cohort), stat='summary', fun.y=median) + 
    stat_summary(fun.y=median, geom="line", aes(color = cohort)) +
    # ggforce::facet_wrap_paginate(~ box, ncol = 2, nrow = 2, page = i, strip.position="top", scales="free_y") +
    # facet_grid(~ cohort) +
    labs(title = paste0("LGA_Age", "_Verified_Data_U01_Olivier", "\n", "By cohort"),
         y = "LGA_Age") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

  # g_rein <-rewards_lga_tograph %>% 
  #   dplyr::filter(active_lever_excel == active_lever_raw,
  #                 inactive_lever_excel == inactive_lever_raw,
  #                 infusions_excel == infusions_raw) %>% 
  #   mutate(session = as.numeric(session) %>% as.factor) %>% 
  #   dplyr::filter(resolution != "FLAG_EXPERIMENT"|is.na(resolution)) %>% 
  #   ggplot(aes(x = session, group = session, fill = reinforcement_schedule)) + 
  #   geom_boxplot(aes_string(y = olivier_lga_measures[i])) + 
  #   # facet_grid(~ reinforcement_schedule) +
  #   labs(title = paste0(gsub("_(?!excel)", " ", gsub("_excel", "", olivier_lga_measures[i]), perl = T), "_Verified_Data_U01_olivier", "\n", "By reinforcement schedule"),
  #        y = gsub("_(?!excel)", " ", gsub("_excel", "", olivier_lga_measures[i]), perl = T)) + 
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  # 
  # g_saline <-rewards_lga_tograph %>% 
  #   dplyr::filter(active_lever_excel == active_lever_raw,
  #                 inactive_lever_excel == inactive_lever_raw,
  #                 infusions_excel == infusions_raw) %>% 
  #   dplyr::filter(resolution != "FLAG_EXPERIMENT"|is.na(resolution) ) %>% 
  #   mutate(session = as.numeric(session) %>% as.factor) %>% 
  #   ggplot(aes(x = session, group = session)) + 
  #   geom_boxplot(aes_string(y = olivier_lga_measures[i])) + 
  #   facet_grid(~ heroin_or_saline) +
  #   labs(title = paste0(gsub("_(?!excel)", " ", gsub("_excel", "", olivier_lga_measures[i]), perl = T), "_Verified_Data_U01_olivier", "\n", "By heroin/saline"),
  #        y = gsub("_(?!excel)", " ", gsub("_excel", "", olivier_lga_measures[i]), perl = T)) + 
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  # 
  # # rewards_lga_tograph %>% select(cohort, self_administration_room) %>% table %>% kable 
  # # rewards_lga_tograph %>% select(cohort, reinforcement_schedule) %>% table %>% kable 
  # # rewards_lga_tograph %>% select(sex, self_administration_box) %>% table %>% kable 
  # # rewards_lga_tograph %>% select(cohort, self_administration_box) %>% table %>% kable 
  # # 
  
  print(g_cohort)
  print(g_sex)
  print(g_dir)
  print(g_dir_sex)
  print(g_age)
  print(g_age_cohort)
  print(g_sex_upto14)
  # print(g_box)
  # print(g_box_x)
  # print(g_rein)
  # print(g_saline)
  # 
}

dev.off()
