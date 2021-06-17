## create csv 

## wfu metadata 
setwd("~/Desktop/Database/csv files/u01_olivier_george_cocaine")
write.csv(rats_allcohorts_weights, "rats_allcohorts_weights_c01_11.csv", row.names = F)

write.csv(oliviercocaine_14_wfu_metadata, "~/Desktop/Database/csv files/u01_olivier_george_cocaine/oliviercocaine_14_wfu_metadata.csv", row.names = F)


# raw data
write.csv(pr_raw_df, "~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/U01/Olivier_George_U01DA043799 (Cocaine)/excel_and_csv_files/cocaine_pr_raw_c01_11_oldnewdirs.csv", row.names = F)


## nida gwas traits

cocaine_xl_sha_df_qc %>% 
  select(cohort, labanimalid, rfid, exp, room, box, rewards_raw, rewards_QC) %>% 
  rename("rewards" = "rewards_raw") %>% 
  subset(rewards_QC == "pass") %>% 
  spread(exp, rewards) %>% 
  # mutate(rewards = replace(rewards, rewards_QC != "pass", NA)) %>% 
  select(-rewards_QC) %>% 
  select(-matches("sha0[1-7]")) %>% 
  mutate(mean_sha_08_10 = rowMeans(select(., starts_with("sha")), na.rm = T)) %>% 
  rename("mean_sha_08_10_box" = "box") %>% 
  select(-matches("^sha")) %>% 
  
  full_join(
    cocaine_xl_lga_df_qc %>% 
      select(cohort, labanimalid, rfid, exp, room, box, rewards_raw, rewards_QC) %>% 
      rename("rewards" = "rewards_raw") %>% 
      subset(rewards_QC == "pass") %>% 
      spread(exp, rewards) %>% 
      select(-rewards_QC) %>% 
      mutate_at(vars(matches("lga1[1-4]$")), list(esc = ~.-lga01)) %>% 
      mutate(mean_lga_esc_11_14 = rowMeans(select(., ends_with("_esc")), na.rm = TRUE)) %>% 
      rename("mean_lga_esc_11_14_box" = "box") %>% 
      select(-matches("^lga")),
   by = c("cohort", "labanimalid", "rfid", "room")
  ) %>%
  mutate(sex = str_extract(labanimalid, "[MF]")) %>%
  bind_rows(read.csv("~/Downloads/db_cocaine_selfadmin_phenotypes.csv", colClasses = "character") %>% 
              mutate(mean_sha_08_10 = as.numeric(mean_sha_08_10),
                     mean_lga_esc_11_14 = as.numeric(mean_lga_esc_11_14))) %>% 
  select(cohort, labanimalid, rfid, sex, room, mean_sha_08_10, mean_lga_esc_11_14, mean_sha_08_10_box, mean_lga_esc_11_14_box) %>% 
  write.csv("~/Desktop/Database/csv files/u01_olivier_george_cocaine/gwas_cocaine_phenotypes_w_c10_11.csv", row.names = F)

read.csv("~/Downloads/db_cocaine_selfadmin_phenotypes.csv", colClasses = "character") %>% 
  distinct(cohort) 


