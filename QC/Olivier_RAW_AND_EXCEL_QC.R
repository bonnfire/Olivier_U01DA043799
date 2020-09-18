# Add more information to selfadmin data and check validity

library(dplyr)
library(data.table)
library(tidyverse)
library(readxl)

# install.packages('splitstackshape')
# install.packages('janitor')
library(splitstackshape)
library(janitor)
library(stringr)


# extract the GWAS database Giordano excel sheet
setwd("~/Dropbox (Palmer Lab)/Olivier_George_U01")

u01.importxlsx <- function(xlname){
  df <- lapply(excel_sheets(path = xlname), read_excel, path = xlname)
  names(df) <- excel_sheets(xlname)
  return(df)
}

## TO DO:  
# import the xl, extract the first table, extract the colorcoded table, assign the logic to create new variables in selfadmin_df
########################
# GWAS New Excel - DrB #
########################
GWAS_database_xl <- u01.importxlsx("GWAS database Giordano NEW.xlsx")
selfadmindatabasexl <- GWAS_database_xl[[1]]
selfadmindatabasexl_subset <- selfadmindatabasexl %>% 
  select(`RAT #`, `Acquis ShA`:`SHOCK`) %>%
  rename_all(tolower) # XX REMOVE SPACES FROM VARNAMES 
names(selfadmindatabasexl_subset) <- gsub(selfadmindatabasexl_subset, )

selfadmin_df_gwascopy <- selfadmin_df %>% 
  mutate(avgsha7_10 = rowMeans(.[, sha07:sha10]),
         avglga11_14 = rowMeans(.[, lga11:lga14]),
         avg_zlga11_14 = 
         ) %>% # rowMeans(a[c('high', 'low')], na.rm=TRUE) XX PICK UP HERE -- NOT SURE IF NA SHOULD BE REMOVED FROM MEAN CALCULATIONS
  mutate(`acquis_sha` = ifelse(avgsha7_10 < 4, "No", "Yes"),
         `escalation` = ifelse(avglga11_14 > 4, "Vulnerable", "Resistant")) %>% 
  group_by(sex) %>% 
  mutate(sex_mean = _ - mean(), 
         sex_sd = sd()) %>% 
  ungroup() %>%
  mutate(add_ind = (maxpr - sex_mean)/sex_df,
         add_ind_cat = ifelse(addind>0.49|addind==0.49, "High AI", "Low AI"))
# `add ind` =IF(BW2>=0.49,"High AI","Low AI")
# `fr` =IF(BT2>=0.49,"High FR","Low FR")
# `pr` =IF(BU2>=0.49,"High PR","Low PR")
# `shock`=IF(BV3>=0.49,"Res Shock","Vuln Shock")

rowMeans(a[c('high', 'low')], na.rm=TRUE)


# check ids from the raw file, find the four 14 digit id's (should be 15)


## DONE: (see documentation in u01_spleenceca_shipment.R) 
# check if spleen and ceca are in the naive
# extract naive ids from wfu 
# extract ids from the spleen and ceca shipment files 
# check if matches 

setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/20190829_WFU_U01_ShippingMaster")
# extract the naive cases from the wfu data in u01_qc.R


##
# missing data from self admin(for gwas)

selfadmin_df %>% group_by(cohort) %>% count()
WFU_Olivier_co_test_df %>% group_by(cohort) %>% count()
setdiff(WFU_Olivier_co_test_df$rfid, selfadmin_df$rfid) %>% length # 198 animals
setdiff(WFU_Olivier_co_test_df$rfid, selfadmin_df$rfid)[which(setdiff(WFU_Olivier_co_test_df$rfid, selfadmin_df$rfid) %>% nchar!=15)] # same four are the ones that don't have 15 characters

# sending spleens but no excel sheet data 
anti_join(olivier_spleen_list_df, selfadmin_df, by = "rfid") 

# data from spleen shipments
setdiff(WFU_Olivier_co_test_df$rfid,olivier_spleen_list_df$rfid)
olivier_spleen_list_df %>% group_by(cohort) %>% count()

## checking for comments 
specificcomments_list_df <- lapply(list, `[[`, 3) %>% rbindlist(fill = T)
dead_list_df <- lapply(list, `[[`, 'dead') %>% rbindlist(fill = T)
# compare here: 
missingdataspleenextractions <- anti_join(olivier_spleen_list_df[which(olivier_spleen_list_df$experiment=="Cocaine"),], selfadmin_df, by = "rfid") %>% dplyr::filter(!is.na(rfid)) %>% select(labanimalid)
specificcomments_list_df %>% dplyr::filter(labanimalid %in% missingdataspleenextractions)










##################
## SHA ########### 
##################
setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/QC")

### QCing raw data

sha_rewards_new %>% 
  mutate_at(vars(exp), as.factor) %>% 
  mutate(sex = str_extract(labanimalid, "[M|F]")) %>% 
  ggplot(aes(exp, rewards, group = labanimalid, color = sex)) + 
  geom_line() +
  facet_grid( ~ cohort) + 
  labs(title = "SHA Rewards New (Raw Only) Directories, For Each Rat") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

sha_rewards_old %>% 
  mutate_at(vars(exp), as.factor) %>% 
  mutate(sex = str_extract(labanimalid, "[M|F]")) %>% 
  ggplot(aes(exp, rewards, group = labanimalid, color = sex)) + 
  geom_line() +
  facet_grid( ~ cohort) + 
  labs(title = "SHA Rewards Old (Raw Only) Directories, For Each Rat") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  

### QCing raw vs excel data


# sha NEW and OLD combine 
## 5/20 XX add valid to sha_rewards_new, remove when not needed
rewards <- rbindlist(
  list(
    "new_sha" = sha_rewards_new_valid,
    "old_sha" = sha_rewards_old,
    "new_lga" = lga_rewards_new_valid,
    "old_lga" = lga_rewards_old,
    "new_pr" = pr_rewards_new_valid,
    "old_pr" = pr_rewards_old,
    "new_shocks" = shock_rewards_new_valid
    
  ),
  idcol = "directory",
  fill = T
)



######### JOIN TO WFU DATABASE 
  
# add notes about missingness (file or dead)

Olivier_Cocaine_df <- WFU_OlivierCocaine_test_df %>%
  select(cohort, rfid, sex, comment, dob) %>%
  # rename("wfu_labanimalid" = "labanimalid") %>%
  mutate(cohort = paste0("C", cohort)) %>%
  dplyr::filter(grepl("^\\d", rfid)) %>% #811 (ignore the blanks and annotations in the excel)
  left_join(., rat_info_allcohort_xl_df[, c("rat", "rfid")], by = "rfid") %>% # add labanimalid number ## 06/08 use rat_info_allcohort_xl_df rather than allcohorts2
  rename("labanimalid" = "rat") %>% 
  mutate(labanimalid = gsub("HS", "", labanimalid)) %>%
  left_join(., ratinfo_list_deaths_processed %>% select(-c("naive", "datedropped")) %>% subset(grepl("surgery", reasoning, ignore.case = T)) %>% subset(!(rfid == "933000320047576"&reasoning=="Died, during surgery")), by = c("rfid", "cohort")) %>% # 811 # deaths/compromises before any experiments XX ASK TEAM AND FIX THE TWO RFID'S 
  left_join(., ratinfo_list_replacements_processed %>% subset(grepl("^RENUMBERED", comment, ignore.case = T)) %>% select(cohort, originalrat, replacement), by = c("tailmark"="originalrat", "cohort")) %>% # replacements, when the animal dies labanimalid changes XX WAITING FOR THEM TO CONFIRM MISSING RFID
  left_join(., ratinfo_list_replacements_processed %>% subset(grepl("Not Renumbered", comment, ignore.case = T)) %>% mutate(comment_replace = paste("Replacing", originalrat, "But", comment)) %>% select(cohort, rfidreplacement, comment_replace), by = c("rfid"="rfidreplacement", "cohort")) %>% 
  mutate(labanimalid = coalesce(labanimalid, replacement),
         tailmark = ifelse(!is.na(tailmark), paste(tailmark, "originally but replaced"), tailmark),
         comment = ifelse(!is.na(reasoning)&is.na(comment), reasoning,
                          ifelse(!is.na(reasoning)&!is.na(comment), paste0(comment, "; ", reasoning), comment))) %>%
  select(-c("replacement", "reasoning")) %>% # replacements, when the animal is the replacement labanimalid changes XX WAITING FOR THEM TO CONFIRM MISSING RFID
  left_join(., computernotes_coc %>% subset(!grepl("cohort_notes", exp)) %>% select(cohort, exp, computernote), by = "cohort") %>% # 21525 (explains missing files for every session, every rat)
  rename("computernote_exp" = "computernote") %>% 
  left_join(., computernotes_coc %>% subset(grepl("cohort_notes", exp)) %>% select(cohort, computernote), by = "cohort") %>%
  left_join(., rewards, by = c("labanimalid", "cohort", "exp")) %>% # 21526 ## ADDING THE RAW REWARDS DATA # M155 WFU_OlivierOxycodone_test_df %>%rename("wfu_labanimalid" = "labanimalid") %>%mutate(cohort = paste0("C", cohort)) %>%dplyr::filter(grepl("^\\d", rfid)) %>% #811 (ignore the blanks and annotations in the excel)left_join(., olivieroxy_excel[, c("labanimalid", "rfid")], by = "rfid") %>% # add labanimalid numberleft_join(., computernotes_oxy, by = "cohort") %>% # 21525 (explains missing files for every session, every rat)subset(cohort == "C01") %>% add_count(labanimalid) %>% rename("n_cnotes_count" = "n") %>% left_join(., rewards, by = c("labanimalid", "cohort", "exp")) %>% add_count(labanimalid) %>% rename("n_crewards_count" = "n") %>% subset(n_cnotes_count != n_crewards_count) %>% View() 
  left_join(., ratinfo_list_deaths_processed %>% select(-c("tailmark", "naive")), by = c("rfid", "cohort")) %>% # 21608 # deaths/compromises # look back on 933000120138753 and 933000320047461 # bc we are trying to use data for as many days as possible, so hteh deatsh table may have repeats
  mutate_at(vars(contains("date")), lubridate::ymd) %>%
  group_by(labanimalid) %>%
  mutate(
    flag = case_when(
      grepl("Died", reasoning) & date >= datedropped ~ "DEAD_EXCLUDE", ## if the animal has died, remove all data on the data and after
      !grepl("Died", reasoning) & date == datedropped ~ "COMP_EXCLUDE" ## if the animal was compromised, only flag that day
    )
  ) %>%
  ungroup() %>%
  mutate(exp_age = difftime(as.POSIXct(date), as.POSIXct(dob), units = "days") %>% as.numeric(.)) %>% 
  select(cohort, rfid, labanimalid, sex, exp, rewards, date, time, filename, tailmark, computernote_exp, computernote, everything())

## exclude some columns that aren't needed
Olivier_Cocaine_df_sql <- Olivier_Cocaine_df %>% 
  select(cohort, rfid, labanimalid, exp, rewards, exp_age, experiment_duration, comment, reasoning, flag) %>% 
  mutate(labanimalid = replace(labanimalid, experiment_duration == 360.5000&exp == "LGA16"&labanimalid=="F516", "F507")) ## TEMPORARY FIX, BUT NEED TO GO BACK AND FIX 

## XX Olivier_Cocaine_df_sql %>% dim should be the same as Olivier_Cocaine_df_sql %>% distinct() %>% dim

ratinfo_list_replacements_processed %>% subset(grepl("^RENUMBERED", comment, ignore.case = T)) 
ori F724 933000320046616, replaced with F737 (933000320046143)


%>% 
  
  left_join(., computernotes_coc, by = "cohort") %>% # 27313 (explains missing files for every session, every rat)
  left_join(., ratinfo_list_replacements_processed %>% subset(grepl("^RENUMBERED", comment, ignore.case = T)) %>% select(cohort, originalrat, replacement), by = c("tailmark"="originalrat", "cohort")) %>% # replacements, when the animal dies labanimalid changes XX WAITING FOR THEM TO CONFIRM MISSING RFID
  left_join(., rewards, by = c("labanimalid", "cohort", "exp")) %>% # 15527 ## ADDING THE RAW REWARDS DATA
  # left_join(.,
  #   allcohorts2 %>% select(labanimalid, rfid, matches("^sha")) %>% distinct() %>%
  #     gather(exp, rewards_excel, sha01:sha10) %>% mutate(exp = toupper(exp)),
  #   by = c("labanimalid", "rfid", "exp")
  # ) %>% ## 5/20 not sure why added this, perhaps to add excel rewards data
  # rename("rewards_raw" = "rewards",
  #        "exp_date" = "date",
  #        "exp_time" = "time") %>% # 15527
  left_join(., ratinfo_list_deaths_processed, by = c("rfid", "cohort")) %>% # deaths/compromises
  mutate_at(vars(contains("date")), lubridate::ymd) %>%
  group_by(labanimalid) %>%
  mutate(
    flag = case_when(
      grepl("Died", reasoning) & date >= datedropped ~ "DEAD_EXCLUDE", ## if the animal has died, remove all data on the data and after
      !grepl("Died", reasoning) & date == datedropped ~ "COMP_EXCLUDE" ## if the animal was compromised, only flag that day
    )
  ) %>%
  ungroup()

Olivier_Cocaine_df %>% select(cohort, rfid, exp, rewards, datedropped, flag) %>% subset(!is.na(flag))



## 08/03/2020
genotyped_ids_1 <- read.csv("genotyped_cocaine_ids_n357.csv")
genotyped_ids_1 <- genotyped_ids_1 %>% 
  mutate_all(as.character) %>% 
  select(-X) # change data type and remove row number column
## WFU_OlivierCocaine_test_df %>% left_join(genotyped_ids_1, by = "rfid") %>% subset(!is.na(project_name)) %>% dim matches the number we expect = 357

## Olivier_Cocaine_df %>% left_join(genotyped_ids_1, by = "rfid") %>% subset(!is.na(project_name)) %>% distinct(rfid) %>% dim 


calc_sa_phenotype <- function(x){
  x <- Olivier_Cocaine_df %>% subset(cohort%in%x) %>% 
    mutate(labanimalid = str_extract(toupper(labanimalid), "[MF]\\d+")) %>% 
    select(cohort, rfid, labanimalid, sex, exp, rewards) %>% 
    distinct() %>% 
    spread(exp, rewards) %>%
    select(matches("cohort|rfid|labanimalid|sex|LGA|PR0[123]|SHA\\d+|SHOCK03")) %>%
    # select(matches("cohort|rfid|labanimalid|sex|LGA(01|1[234])|PR0[23]|SHA\\d+|SHOCK03")) %>%
    mutate(PR_max = pmax(PR02, PR03, na.rm = F)) %>% # exclude animal if the PR02 03 is not complete data
    group_by(sex) %>% 
    mutate(LGA01_mean = mean(LGA01, na.rm = T),
           LGA01_sd = sd(LGA01, na.rm = T),
           PR_max_mean = mean(PR_max, na.rm = T),
           PR_max_sd = sd(PR_max, na.rm = T),
           SHOCK_mean = mean(SHOCK03, na.rm = T),
           SHOCK_sd = sd(SHOCK03, na.rm = T)) %>% 
    ungroup() %>% 
    mutate_at(vars(matches("LGA1[234]")), list(esc = ~.-LGA01_mean)) %>%
    mutate_at(vars(ends_with("_esc")), ~./LGA01_sd) %>%
    mutate(PR_index = (PR_max - PR_max_mean)/PR_max_sd,
           SHOCK_index = (SHOCK03 - SHOCK_mean)/SHOCK_sd) %>%
    # rowwise() %>% 
    mutate(ind_esc_mean = rowMeans(select(., matches("LGA1[234]_esc")), na.rm = TRUE)) %>% 
    group_by(sex) %>% 
    mutate(esc_mean = mean(ind_esc_mean, na.rm = T),
           esc_sd = sd(LGA01, na.rm = T)) %>%  
    ungroup() %>% 
    mutate(esc_index = (ind_esc_mean - esc_mean)/esc_sd) %>% 
    mutate(addiction_index = rowMeans(select(., ends_with("index")), na.rm = TRUE)) %>% 
    mutate(SHA_last3_mean = rowMeans(select(., matches("SHA(0[89]|10)")), na.rm = F)) 
  # %>% 
    # select(matches("cohort|rfid|labanimalid|sex|SHA_last3_mean|PR_index|esc_index|SHOCK_index|addiction_index")) 
  
  # "How long does it take the animal to hit 5 rewards in SHA?"
  # y <- Olivier_Cocaine_df %>% subset(cohort == x|cohort%in% x) %>% 
  #   mutate(labanimalid = str_extract(toupper(labanimalid), "[MF]\\d+")) %>% 
  #   select(cohort, rfid, labanimalid, sex, exp, rewards) %>% 
  #   distinct() %>% 
  #   subset(grepl("SHA", exp)) %>% 
  #   group_by(rfid) %>% 
  #   mutate(SHA_daysto5 = cumsum(case_when(rewards >= 5 ~ 1, TRUE ~ 0))) %>% 
  #   dplyr::filter(rewards >= 5|SHA_daysto5 == 0) %>% slice(1) %>% 
  #   ungroup() %>% 
  #   select(SHA_daysto5)
  #   
  # final_phenotypes <- cbind(x, y)
  # 
  # return(final_phenotypes)
  return(x)
}


setwd("~/Desktop/Database/csv files/u01_olivier_george_cocaine")

# cohort 1
cohort1_sa_phenotype_df <- calc_sa_phenotype("C01") 
cohort1_sa_phenotype_df %>% subset(is.na(labanimalid))
# cohort1_sa_phenotype_df %<>% mutate(labanimalid = replace(labanimalid, rfid == "933000120124702", "F17")) # have a dead table, instead of fixing here
cohort1_sa_phenotype_df <- cohort1_sa_phenotype_df[gtools::mixedorder(cohort1_sa_phenotype_df$labanimalid),]
# check for no duplicates
cohort1_sa_phenotype_df %>% get_dupes(rfid)
write.csv(cohort1_sa_phenotype_df, file = "cohort1_sa_phenotype.csv", row.names = F)  

# cohort 2
cohort2_sa_phenotype_df <- calc_sa_phenotype("C02")
cohort2_sa_phenotype_df %>% subset(is.na(labanimalid))
cohort2_sa_phenotype_df <- cohort2_sa_phenotype_df[gtools::mixedorder(cohort2_sa_phenotype_df$labanimalid),]
# check for no duplicates
cohort2_sa_phenotype_df %>% get_dupes(rfid)
write.csv(cohort2_sa_phenotype_df, file = "cohort2_sa_phenotype.csv", row.names = F)  

# cohort 3
cohort3_sa_phenotype_df <- calc_sa_phenotype("C03")
cohort3_sa_phenotype_df %>% subset(is.na(labanimalid)|!grepl("^[MF]\\d+$", labanimalid))
cohort3_sa_phenotype_df <- cohort3_sa_phenotype_df[gtools::mixedorder(cohort3_sa_phenotype_df$labanimalid),]
# check for no duplicates
cohort3_sa_phenotype_df %>% get_dupes(rfid)
write.csv(cohort3_sa_phenotype_df, file = "cohort3_sa_phenotype.csv", row.names = F) 

# cohort 4
cohort4_sa_phenotype_df <- calc_sa_phenotype("C04")
cohort4_sa_phenotype_df %>% subset(is.na(labanimalid)|!grepl("^[MF]\\d+$", labanimalid))
cohort4_sa_phenotype_df <- cohort4_sa_phenotype_df[gtools::mixedorder(cohort4_sa_phenotype_df$labanimalid),]
# check for no duplicates
cohort4_sa_phenotype_df %>% get_dupes(rfid)
write.csv(cohort4_sa_phenotype_df, file = "cohort4_sa_phenotype.csv", row.names = F) 

# cohort 5
cohort5_sa_phenotype_df <- calc_sa_phenotype("C05")
cohort5_sa_phenotype_df %>% subset(is.na(labanimalid)|!grepl("^[MF]\\d+$", labanimalid))
cohort5_sa_phenotype_df <- cohort5_sa_phenotype_df[gtools::mixedorder(cohort5_sa_phenotype_df$labanimalid),]
# check for no duplicates
cohort5_sa_phenotype_df %>% get_dupes(rfid)
write.csv(cohort5_sa_phenotype_df, file = "cohort5_sa_phenotype.csv", row.names = F) 

# cohort 6 # ABORTED
cohort6_sa_phenotype_df <- calc_sa_phenotype("C06")
cohort6_sa_phenotype_df %>% subset(is.na(labanimalid)|!grepl("^[MF]\\d+$", labanimalid))
cohort6_sa_phenotype_df <- cohort6_sa_phenotype_df[gtools::mixedorder(cohort6_sa_phenotype_df$labanimalid),]
# check for no duplicates
cohort6_sa_phenotype_df %>% get_dupes(rfid)
write.csv(cohort6_sa_phenotype_df, file = "cohort6_sa_phenotype.csv", row.names = F) 

# cohort 7
cohort7_sa_phenotype_df <- calc_sa_phenotype("C07")
cohort7_sa_phenotype_df %>% subset(is.na(labanimalid)|!grepl("^[MF]\\d+$", labanimalid))
cohort7_sa_phenotype_df <- cohort7_sa_phenotype_df[gtools::mixedorder(cohort7_sa_phenotype_df$labanimalid),]
# check for no duplicates
cohort7_sa_phenotype_df %>% get_dupes(rfid)
write.csv(cohort7_sa_phenotype_df, file = "cohort7_sa_phenotype.csv", row.names = F) 

# cohort 8
cohort8_sa_phenotype_df <- calc_sa_phenotype("C08")
cohort8_sa_phenotype_df %>% subset(is.na(labanimalid)|!grepl("^[MF]\\d+$", labanimalid))
cohort8_sa_phenotype_df <- cohort8_sa_phenotype_df[gtools::mixedorder(cohort8_sa_phenotype_df$labanimalid),]
# check for no duplicates
cohort8_sa_phenotype_df %>% get_dupes(rfid)
write.csv(cohort8_sa_phenotype_df, file = "cohort8_sa_phenotype.csv", row.names = F) 



## Olivier_Cocaine_df doesn't have cohort 9 and on yet 
# cohort 9
cohort9_sa_phenotype_df <- calc_sa_phenotype("C09")
cohort9_sa_phenotype_df %>% subset(is.na(labanimalid)|!grepl("^[MF]\\d+$", labanimalid))
cohort9_sa_phenotype_df <- cohort9_sa_phenotype_df[gtools::mixedorder(cohort9_sa_phenotype_df$labanimalid),]
# check for no duplicates
cohort9_sa_phenotype_df %>% get_dupes(rfid)
write.csv(cohort9_sa_phenotype_df, file = "cohort9_sa_phenotype.csv", row.names = F) 







## 08/05/2020 use this to plot all cohorts
Olivier_Cocaine_C01_09 <-  Olivier_Cocaine_df %>% 
  mutate(labanimalid = str_extract(toupper(labanimalid), "[MF]\\d+")) %>% 
  select(cohort, rfid, labanimalid, sex, box, exp, rewards) %>% 
  distinct() %>% 
  spread(exp, rewards) %>%
  select(matches("cohort|rfid|labanimalid|sex|box|LGA(01|1[234])|PR0[23]|SHA\\d+|SHOCK03")) %>%
  mutate(PR_max = pmax(PR02, PR03, na.rm = F)) %>% # exclude animal if the PR02 03 is not complete data
  group_by(sex) %>% 
  mutate(LGA01_mean = mean(LGA01, na.rm = T),
         LGA01_sd = sd(LGA01, na.rm = T),
         PR_max_mean = mean(PR_max, na.rm = T),
         PR_max_sd = sd(PR_max, na.rm = T),
         SHOCK_mean = mean(SHOCK03, na.rm = T),
         SHOCK_sd = sd(SHOCK03, na.rm = T)) %>% 
  ungroup() %>% 
  mutate_at(vars(matches("LGA1[234]")), list(esc = ~.-LGA01_mean)) %>%
  mutate_at(vars(ends_with("_esc")), ~./LGA01_sd) %>%
  mutate(PR_index = (PR_max - PR_max_mean)/PR_max_sd,
         SHOCK_index = (SHOCK03 - SHOCK_mean)/SHOCK_sd) %>%
  # rowwise() %>% 
  mutate(ind_esc_mean = rowMeans(select(., matches("LGA1[234]_esc")), na.rm = TRUE)) %>% 
  group_by(sex) %>% 
  mutate(esc_mean = mean(ind_esc_mean, na.rm = T),
         esc_sd = sd(LGA01, na.rm = T)) %>%  
  ungroup() %>% 
  mutate(esc_index = (ind_esc_mean - esc_mean)/esc_sd) %>% 
  mutate(addiction_index = rowMeans(select(., ends_with("index")), na.rm = TRUE)) %>% 
  mutate(SHA_last3_mean = rowMeans(select(., matches("SHA(0[89]|10)")), na.rm = F)) %>% 
  select(matches("cohort|rfid|labanimalid|sex|box|SHA_last3_mean|PR_index|esc_index|SHOCK_index|addiction_index")) %>% 
  left_join(Olivier_Cocaine_df %>% # "How long does it take the animal to hit 5 rewards in SHA?"
      mutate(labanimalid = str_extract(toupper(labanimalid), "[MF]\\d+")) %>% 
      select(cohort, rfid, labanimalid, sex, exp, rewards) %>% 
      distinct() %>% 
      subset(grepl("SHA", exp)) %>% 
      group_by(rfid) %>% 
      mutate(SHA_daysto5 = cumsum(case_when(rewards >= 5 ~ 1, TRUE ~ 0))) %>% 
      dplyr::filter(rewards >= 5|SHA_daysto5 == 0) %>% slice(1) %>% 
      ungroup() %>% 
      select(rfid, SHA_daysto5), by = "rfid") %>%  # left_join instead of cbind, bc cbind has cohort12, from wfu but no data
  clean_names()

## troubleshoot cohort 4
Olivier_Cocaine_df %>% 
  subset(cohort == "C04") %>% 
  mutate(labanimalid = str_extract(toupper(labanimalid), "[MF]\\d+")) %>% 
  select(cohort, rfid, labanimalid, sex, exp, rewards) %>% 
  distinct() %>% 
  spread(exp, rewards) %>%
  select(matches("cohort|rfid|labanimalid|sex|LGA(01|1[234])")) %>%
  group_by(sex) %>% 
  mutate(LGA01_mean = mean(LGA01, na.rm = T),
         LGA01_sd = sd(LGA01, na.rm = T)) %>% 
  ungroup() %>% 
  mutate_at(vars(matches("LGA1[234]")), list(esc = ~.-LGA01_mean)) %>%
  mutate_at(vars(ends_with("_esc")), ~./LGA01_sd) %>%
  mutate(PR_index = (PR_max - PR_max_mean)/PR_max_sd,
         SHOCK_index = (SHOCK03 - SHOCK_mean)/SHOCK_sd) %>%
  # rowwise() %>% 
  mutate(ind_esc_mean = rowMeans(select(., matches("LGA1[234]_esc")), na.rm = TRUE)) %>% 
  group_by(sex) %>% 
  mutate(esc_mean = mean(ind_esc_mean, na.rm = T),
         esc_sd = sd(LGA01, na.rm = T)) %>%  
  ungroup() %>% 
  mutate(esc_index = (ind_esc_mean - esc_mean)/esc_sd) %>% 
  mutate(addiction_index = rowMeans(select(., ends_with("index")), na.rm = TRUE)) %>% 
  mutate(SHA_last3_mean = rowMeans(select(., matches("SHA(0[89]|10)")), na.rm = F)) %>% 
  select(matches("cohort|rfid|labanimalid|sex|SHA_last3_mean|PR_index|esc_index|SHOCK_index|addiction_index"))


## Extract the Excel data for the indices
# 08/07/2020
setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE")

Olivier_phenotype_xl <- openxlsx::read.xlsx("Cocaine_GWAS_ADDIND.xlsx") %>% 
  clean_names() %>% 
  select_if(~sum(!is.na(.)) > 0) %>% 
  select(-matches("x\\d+$")) %>%  # manually removing, after looking at column contents on 08/07/2020
  rename("labanimalid" = "rat")

Olivier_Cocaine_C01_09 %>% 
  left_join(Olivier_phenotype_xl %>% 
              select(-cohort), by = "rfid") %>% 
  subset(!is.na(labanimalid.y)) %>% 
  subset(cohort == "C03")

Olivier_Cocaine_C01_09 %>% 
  left_join(Olivier_phenotype_xl %>% 
              select(-cohort), by = "rfid") %>% 
  subset(cohort == "C04")




### EXTRACT THE EXCEL INDICES (AND INTERMEDIATE VALUES FOR COHORTS 01-07)
setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE")
cocaine_gwas_xl <- u01.importxlsx("Addiction indices for C01-C07 Cocaine.xlsx") %>% 
  lapply(function(x){
    if(grepl("F701", x[1,1])){ # if cohort 7 table
      names(x)[1] <- "rats"
    }
    x <- x %>% clean_names
    i <- 1
    if(names(x)[1] == "x1"){
      names(x) <- x[i,]
      x <- x %>% clean_names
      x <- x[-i, ]
      i = i + 1 # go down the rows until the test expression is false
    }
    
    x <- x %>%
      select(starts_with("rat"), starts_with("fr"), starts_with("z"), starts_with("add"))    
    return(x)
  }) 

# extract the columns for each cohort, rather than coding the rename, and extract the columns at the end of the summary end of the sheet # change to the first occurence for c03
cocaine_gwas_xl$C01<- cocaine_gwas_xl$C01[, c("rats", "fr_zscore_2", "z_score_pr_2", "z_score_shock_2", "addiction_index")]
cocaine_gwas_xl$C02<- cocaine_gwas_xl$C02[, c("rats", "fr_zscore", "z_score_pr_2", "z_score_shock_2", "addiction_index")]
cocaine_gwas_xl$C03<- cocaine_gwas_xl$C03[, c("rats_1", "z_intake_32", "z_pr_37", "z_shock_43", "addiction_index")]
cocaine_gwas_xl$C04<- cocaine_gwas_xl$C04[, c("rat", "z_intake_2", "z_pr_2", "z_shock_2", "addiction_index")]
cocaine_gwas_xl$C05<- cocaine_gwas_xl$C05[, c("rats_1", "z_intake_48", "z_pr_49", "z_shock_50", "addiction_index")]
cocaine_gwas_xl$C07<- cocaine_gwas_xl$C07[, c("rats", "z_esc_45", "zpr_46", "z_shock_47", "add_ind")]


cocaine_gwas_xl_df <- cocaine_gwas_xl %>% 
  lapply(function(x){
    # names(x) <- c("labanimalid", "z_esc", "z_pr", "z_shock", "add_ind")
    names(x) <- c("labanimalid", "escalation_zscore", "pr_zscore", "shock_zscore", "addiction_index")
    return(x)
  }) %>% 
  rbindlist(idcol = "cohort", fill = T) %>%  
  mutate(labanimalid = str_extract(labanimalid, "[MF]\\d+")) %>% 
  mutate_at(vars(-matches("cohort|labanimal")), as.numeric) %>% 
  subset(!is.na(labanimalid)) %>% 
  left_join(rat_info_allcohort_xl_df[, c("rfid", "labanimalid")], by = "labanimalid") %>%  # extract rfid from mapping files
  left_join(Olivier_Cocaine_C01_09[, c("rfid", "labanimalid", "box")])


  

#### CREATE RAW VS EXCEL SHEETS

setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE")
cocaine_intermediates_xl_lga <- u01.importxlsx("Addiction indices for C01-C07 Cocaine.xlsx") %>% 
  lapply(function(x){
    if(grepl("F701", x[1,1])){ # if cohort 7 table
      names(x)[1] <- "rats"
    }
    x <- x %>% clean_names
    i <- 1
    if(names(x)[1] == "x1"){
      names(x) <- x[i,]
      x <- x %>% clean_names
      x <- x[-i, ]
      i = i + 1 # go down the rows until the test expression is false
    }
    
    x <- x %>%
      select(starts_with("rat"), starts_with("lg"), starts_with("day"))
    return(x)
  }) 
cocaine_intermediates_xl_lga$C01<- cocaine_intermediates_xl_lga$C01 %>% select(matches("rats$"), matches("^day_\\d+$"), -matches("^day_\\d+_\\d+"))
cocaine_intermediates_xl_lga$C02<- cocaine_intermediates_xl_lga$C02 %>% select(matches("rats$"), matches("^day_\\d+$"), -matches("^day_\\d+_\\d+"))
cocaine_intermediates_xl_lga$C03<- cocaine_intermediates_xl_lga$C03 %>% select(matches("rats_1$"), matches("^day_\\d_\\d$"), matches("^day_([9]|1\\d)_[1]\\d")) 
cocaine_intermediates_xl_lga$C04<- cocaine_intermediates_xl_lga$C04 %>% select(matches("rat$"), matches("^lg_a\\d+$"), -matches("^lg_a\\d+_2"))
cocaine_intermediates_xl_lga$C05<- cocaine_intermediates_xl_lga$C05 %>% select(matches("rats_1$"), matches("^day_\\d_\\d$"), matches("^day_([9]|1\\d)_[1]\\d")) 
cocaine_intermediates_xl_lga$C07<- cocaine_intermediates_xl_lga$C07 %>% select(matches("rats$"), matches("^lg_a\\d_\\d$"), matches("^lg_a([9]|1\\d)_[1]\\d")) 

# join to raw, wait for fixes, and then join to corrected values
cocaine_intermediates_xl_lga_df <- cocaine_intermediates_xl_lga %>% 
  lapply(function(x){
    names(x) <- c("labanimalid", paste0("lga_", str_pad(1:14, "2", "left", "0"))) 
    return(x)
  }) %>% 
  rbindlist(idcol = "cohort", fill = T) %>%  
  mutate(labanimalid = str_extract(labanimalid, "[MF]\\d+")) %>% 
  mutate_at(vars(-matches("cohort|labanimal")), as.numeric) %>% 
  subset(!is.na(labanimalid)) %>% 
  left_join(rat_info_allcohort_xl_df[, c("rfid", "labanimalid")], by = "labanimalid") %>% 
  gather("exp", "rewards_xl", -cohort, -labanimalid, -rfid) %>% 
  left_join(Olivier_Cocaine_df %>% 
              mutate(labanimalid = str_extract(toupper(labanimalid), "[MF]\\d+")) %>% 
              select(cohort, rfid, labanimalid, sex, exp, rewards) %>% distinct() %>% 
              spread(exp, rewards) %>% select(cohort, rfid, labanimalid, sex, matches("LGA(0[1-9]|1[0-4])")) %>% 
              gather("exp", "rewards_raw", -cohort, -rfid, -labanimalid, -sex) %>% 
              mutate(exp = gsub("LGA", "lga_", exp)),
            by = c("cohort", "labanimalid", "rfid", "exp")) %>% 
  mutate(rewards_QC_diff = rewards_xl - rewards_raw,
         rewards_QC = ifelse(rewards_QC_diff == 0, "pass", "fail"))

### QC'ING AND SHOWING THE ONES THAT FAIL QC
# copy for Palmer Lab, long

cocaine_qc_long <- cocaine_intermediates_xl_df %>% subset(rewards_QC == "fail") %>% 
  left_join(Olivier_Cocaine_df %>% 
              mutate(exp = gsub("LGA", "lga_", exp)) %>% 
              select(cohort, labanimalid, rfid, exp, sex, filename), 
            by = c("cohort", "labanimalid", "rfid", "exp", "sex"))

# copy for Olivier lab, wide 

cocaine_qc_wide <- cocaine_intermediates_xl_df %>% 
  subset(rewards_QC == "fail") %>% 
    left_join(Olivier_Cocaine_df %>% 
              mutate(exp = gsub("LGA", "lga_", exp)) %>% 
              select(cohort, labanimalid, rfid, exp, sex, filename), 
            by = c("cohort", "labanimalid", "rfid", "exp", "sex")) %>% 
  spread(exp, rewards_xl) %>% 
  mutate(labanimalid_num = parse_number(labanimalid)) %>% 
  arrange(cohort, sex, labanimalid_num) %>% select(-labanimalid_num)

library(openxlsx)

## Split data apart by a grouping variable;
##   makes a named list of tables
cocaine_qc_wide_bycohort <- split(cocaine_qc_wide, cocaine_qc_wide$cohort)

## Create a blank workbook
wb <- createWorkbook()

## Loop through the list of split tables as well as their names
##   and add each one as a sheet to the workbook
Map(function(data, name){
  
  addWorksheet(wb, name)
  writeData(wb, name, data)
  
}, cocaine_qc_wide_bycohort, names(cocaine_qc_wide_bycohort))


## Save workbook to working directory
saveWorkbook(wb, file = "cocaine_qc_bycohort.xlsx", overwrite = TRUE)

## fix the data using the lab's response 08/20/2020
setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE")
cocaine_qc_decisions_lga <- u01.importxlsx("decision_cocaine_qc_bycohort (updated).xlsx")[-1] %>% 
  rbindlist()
cocaine_qc_decisions_lga %>% select(source) %>% table()

### find the QC NA cases
# cocaine_qc_wide_na <- cocaine_intermediates_xl_df %>% 
#   subset(is.na(rewards_QC)) %>% 
#   left_join(Olivier_Cocaine_df %>% 
#               mutate(exp = gsub("LGA", "lga_", exp)) %>% 
#               select(cohort, labanimalid, rfid, exp, sex, filename), 
#             by = c("cohort", "labanimalid", "rfid", "exp", "sex")) %>% 
#   spread(exp, rewards_xl) %>% 
#   mutate(labanimalid_num = parse_number(labanimalid)) %>% 
#   arrange(cohort, sex, labanimalid_num) %>% select(-labanimalid_num)
# 
# library(openxlsx)
# 
# ## Split data apart by a grouping variable;
# ##   makes a named list of tables
# cocaine_qc_wide_bycohort <- split(cocaine_qc_wide, cocaine_qc_wide$cohort)
# 
# ## Create a blank workbook
# wb <- createWorkbook()
# 
# ## Loop through the list of split tables as well as their names
# ##   and add each one as a sheet to the workbook
# Map(function(data, name){
#   
#   addWorksheet(wb, name)
#   writeData(wb, name, data)
#   
# }, cocaine_qc_wide_bycohort, names(cocaine_qc_wide_bycohort))
# 
# 
# ## Save workbook to working directory
# saveWorkbook(wb, file = "cocaine_qc_bycohort.xlsx", overwrite = TRUE)


### join and correct the raw values
cocaine_intermediates_xl_lga_corrected <- cocaine_intermediates_xl_lga_df %>% 
  left_join(cocaine_qc_decisions_lga %>% 
              select(cohort, labanimalid, rfid, source, starts_with("lga")) %>% 
              gather("exp", "rewards", -cohort, -labanimalid, -rfid, -source) %>% 
              subset(!is.na(rewards)), 
            by = c("cohort", "labanimalid", "rfid", "exp")) %>% 
  naniar::replace_with_na_all(condition = ~.x %in% c("NA", "NaN")) %>% 
  mutate(rewards = ifelse(is.na(rewards_QC)&!is.na(rewards_xl), rewards_xl,  
                          ifelse(rewards_QC == "pass"|source == "raw"|is.na(rewards_xl), rewards_raw, rewards_xl)),
         sex = str_extract(labanimalid, "[MF]")) %>% ## make sure that there are no gaps 
         # source = ifelse(is.na(source), "raw", source)) %>% 
  select(cohort, rfid, labanimalid, sex, exp, rewards) # remove source from data



# cohorts 01-07 (LGA)
cohort01_07_cocaine_lga <- cocaine_intermediates_xl_lga_corrected %>% 
  mutate(exp = paste0(exp, "_rewards")) %>% 
  spread(exp, rewards) %>% 
  mutate_at(vars(matches("lga_(0[2-9]|1[0-4])_rewards$")), list(esc = ~.-lga_01_rewards)) %>% 
  mutate(esc11_14_mean = rowMeans(select(., matches("1[1-4]_rewards_esc")), na.rm = TRUE)) %>% 
  left_join(subjects_exp_age %>% select(matches("cohort|rfid|labanimalid|^lga")), c("cohort", "rfid", "labanimalid")) %>%  # join the age data xx *new line* # join the active and inactive lever presses
  select(cohort, rfid, sex, labanimalid, everything())
# quick qc before upload 
cohort01_07_cocaine_lga %>% get_dupes(rfid)
cohort01_07_cocaine_lga %>% ggplot() + geom_density(aes(x = esc11_14_mean)) + facet_grid(rows = "cohort")
# writing onto desktop folder bc Dropbox does not work
setwd("~/Desktop/Database/csv files/u01_olivier_george_cocaine/")
write.csv(cohort01_07_cocaine_lga, file = "cohort01_07_lga_phenotypes.csv", row.names = F)


# cohort01-08 (SHA)
cohort01_08_cocaine_sha_qc <- rbind(sha_rewards_new_valid %>% subset(grepl("SHA(0[89]|10)", exp)) %>% select(cohort, labanimalid, exp, rewards, filename),
      sha_rewards_old %>% subset(grepl("SHA(0[89]|10)", exp)) %>% mutate(filename = gsub("(.*/){4}", "", filename)) %>% select(cohort, labanimalid, exp, rewards, filename)) %>% 
  mutate(labanimalid = replace(labanimalid, labanimalid == "F414"&exp == "SHA09", "F424"),
         rewards = replace(rewards, labanimalid == "F511"&exp =="SHA08", NA)) %>% # notes from rat info and cocaine self admin excel 
  rowwise() %>% 
  mutate(labanimalid = replace(labanimalid, cohort=="C05"&grepl("F1", labanimalid), gsub("F1", "F5", labanimalid))) %>% 
  left_join(rat_info_allcohort_xl_df[c("labanimalid", "rfid")], by = "labanimalid") %>%
  # mutate(rewards = replace(rewards, rfid %in% unique(compromised_rats[grep("died", compromised_rats$death_comment, ignore.case = T),]$rfid), NA)) %>%   # XX make values NA if dead after date
  rename("rewards_raw" =  "rewards") %>% 
  left_join(olivierxl_df %>% select(rfid, labanimalid, matches("^sha(0[89]|10)$")) %>% gather("exp", "rewards_xl", -rfid, -labanimalid) %>% mutate(exp = toupper(exp)), by = c("labanimalid", "rfid", "exp")) %>%
  mutate(sex = str_extract(labanimalid, "[MF]"),
         rewards_QC_diff = rewards_xl - rewards_raw,
         rewards_QC = ifelse(rewards_QC_diff == 0, "pass", "fail")) %>% 
  ungroup()

setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE")
cohort01_08_cocaine_sha_qc %>% subset(rewards_QC == "fail") %>% 
  spread(exp, rewards_xl) %>% 
  mutate(labanimalid_num = parse_number(labanimalid)) %>% 
  arrange(cohort, sex, labanimalid_num) %>% select(-labanimalid_num) %>% 
  openxlsx::write.xlsx(file = "cocaine_qc_sha.xlsx")

## create the corrected objects (SHA)  
setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE")
cocaine_qc_decisions_sha <- u01.importxlsx("decision_cocaine_qc_sha.xlsx")[[1]]
  
cocaine_intermediates_xl_sha_corrected <- cohort01_08_cocaine_sha_qc %>% 
  left_join(cocaine_qc_decisions_sha %>% 
              select(cohort, labanimalid, rfid, source, starts_with("SHA")) %>% 
              gather("exp", "rewards", -cohort, -labanimalid, -rfid, -source) %>% ## maybe add F52[345]
              subset(!is.na(rewards)), 
            by = c("cohort", "labanimalid", "rfid", "exp")) %>% 
  naniar::replace_with_na_all(condition = ~.x %in% c("NA", "NaN")) %>% 
  mutate(rewards = ifelse(is.na(rewards_QC)&!is.na(rewards_xl), rewards_xl,  
                          ifelse(rewards_QC == "pass"|source == "raw"|is.na(rewards_xl), rewards_raw, rewards_xl))) %>% ## make sure that there are no gaps 
  # source = ifelse(is.na(source), "raw", source)) %>% 
  subset(!(is.na(rewards)&grepl("F52[345]", labanimalid))) %>% 
  select(cohort, rfid, labanimalid, sex, exp, rewards) # remove source from data


# mean sha
cohort01_07_cocaine_sha <- cocaine_intermediates_xl_sha_corrected %>% 
  mutate(exp = tolower(exp)) %>% 
  spread(exp, rewards) %>% 
  mutate(mean_sha_last3 = rowMeans(select(., starts_with("sha")), na.rm = T)) 

  


# exclude the un-qc'ed id's for the database 09/14/2020
cohort01_08_cocaine_sha <- cohort01_08_cocaine_sha_qc %>% 
  subset(rfid %in% unique(cohort01_08_cocaine_sha_qc[which(cohort01_08_cocaine_sha_qc$rewards_QC == "pass"),]$rfid)) %>% 
  # group_by(rfid, cohort, labanimalid, sex) %>% 
  # summarize(mean_sha_last3 = mean(rewards_raw, is.na = T)) ## XX come back to this -- fix because it is taking the group by incorrectly


# qc raw vs excel ## extract multiple cohorts at once for the phenotypes 
calc_sa_phenotype(c("C01", "C02")) %>% select(cohort, rfid, labanimalid, sex, matches("^PR0[123]$"), SHOCK03, matches("SHA_last3_mean")) %>% head(3)


# cohort01-08 (PR)
cohort01_08_cocaine_pr_qc <- rbind(pr_rewards_new_valid %>% subset(grepl("PR0[123]", exp)) %>% select(cohort, labanimalid, exp, rewards, filename),
                                    pr_rewards_old %>% subset(grepl("PR0[123]", exp)) %>% mutate(filename = gsub("(.*/){4}", "", filename)) %>% select(cohort, labanimalid, exp, rewards, filename)) %>% 
  left_join(rat_info_allcohort_xl_df[c("labanimalid", "rfid")], by = "labanimalid") %>% # give rfid
  rename("rewards_raw" =  "rewards") %>%
  subset(!is.na(rfid)) %>% 
  right_join(olivierxl_df %>% select(cohort, rfid, labanimalid, matches("^pr0[123]$")) %>% gather("exp", "rewards_xl", -rfid, -labanimalid, -cohort) %>% mutate(exp = toupper(exp)), by = c("cohort", "labanimalid", "rfid", "exp")) %>%
  mutate(sex = str_extract(labanimalid, "[MF]"),
         rewards_QC_diff = rewards_xl - rewards_raw,
         rewards_QC = ifelse(rewards_QC_diff == 0, "pass", "fail"))

setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE")
cohort01_08_cocaine_pr_qc %>% subset(rewards_QC == "fail") %>% 
  spread(exp, rewards_xl) %>% 
  mutate(labanimalid_num = parse_number(labanimalid)) %>% 
  arrange(cohort, sex, labanimalid_num) %>% select(-labanimalid_num) %>% 
  openxlsx::write.xlsx(file = "cocaine_qc_pr.xlsx")

## create the corrected objects (PR)  
setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE")
cocaine_qc_decisions_pr <- u01.importxlsx("decision_cocaine_qc_pr.xlsx")[[1]]

cocaine_intermediates_xl_pr_corrected <- cohort01_08_cocaine_pr_qc %>% 
  left_join(cocaine_qc_decisions_pr %>% 
              select(cohort, labanimalid, rfid, source, starts_with("PR")) %>% 
              gather("exp", "rewards", -cohort, -labanimalid, -rfid, -source) %>% ## maybe add F52[345]
              subset(!is.na(rewards)), 
            by = c("cohort", "labanimalid", "rfid", "exp")) %>% 
  naniar::replace_with_na_all(condition = ~.x %in% c("NA", "NaN")) %>% 
  mutate(rewards = ifelse(is.na(rewards_QC)&!is.na(rewards_xl), rewards_xl,  
                          ifelse(rewards_QC == "pass"|source == "raw"|is.na(rewards_xl), rewards_raw, rewards_xl))) %>% ## make sure that there are no gaps 
  # source = ifelse(is.na(source), "raw", source)) %>% 
  subset(!is.na(rfid)) %>%
  mutate(exp = tolower(exp)) %>% 
  select(cohort, rfid, labanimalid, sex, exp, rewards) %>%  # remove source from data
  spread(exp, rewards)






# cohort01-08 (SHOCKS)
cohort01_08_cocaine_shock_qc <- shock_rewards_new_valid %>% subset(grepl("SHOCK0[3]", exp)) %>% select(cohort, labanimalid, exp, rewards, filename) %>% # no need for rbind bc there are no old files for shock rewards 
  left_join(rat_info_allcohort_xl_df[c("labanimalid", "rfid")], by = "labanimalid") %>%
  rename("rewards_raw" =  "rewards") %>% 
  left_join(olivierxl_df %>% select(rfid, labanimalid, matches("^shock0[3]$")) %>% gather("exp", "rewards_xl", -rfid, -labanimalid) %>% mutate(exp = toupper(exp)), by = c("labanimalid", "rfid", "exp")) %>%
  mutate(sex = str_extract(labanimalid, "[MF]"),
         rewards_QC_diff = rewards_xl - rewards_raw,
         rewards_QC = ifelse(rewards_QC_diff == 0, "pass", "fail"))

setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE")
cohort01_08_cocaine_shock_qc %>% subset(rewards_QC == "fail") %>% 
  spread(exp, rewards_xl) %>% 
  mutate(labanimalid_num = parse_number(labanimalid)) %>% 
  arrange(cohort, sex, labanimalid_num) %>% select(-labanimalid_num) %>% 
  openxlsx::write.xlsx(file = "cocaine_qc_shock.xlsx")

## create the corrected objects (SHOCK)  
setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE")
cocaine_qc_decisions_shock <- u01.importxlsx("decision_cocaine_qc_shock.xlsx")[[1]]

cocaine_intermediates_xl_shock_corrected <- cohort01_08_cocaine_shock_qc %>% 
  left_join(cocaine_qc_decisions_shock %>% 
              select(cohort, labanimalid, rfid, source, starts_with("SHOCK")) %>% 
              gather("exp", "rewards", -cohort, -labanimalid, -rfid, -source) %>% ## maybe add F52[345]
              subset(!is.na(rewards)), 
            by = c("cohort", "labanimalid", "rfid", "exp")) %>% 
  naniar::replace_with_na_all(condition = ~.x %in% c("NA", "NaN")) %>% 
  mutate(rewards = ifelse(is.na(rewards_QC)&!is.na(rewards_xl), rewards_xl,  
                          ifelse(rewards_QC == "pass"|source == "raw"|is.na(rewards_xl), rewards_raw, rewards_xl))) %>% ## make sure that there are no gaps 
  # source = ifelse(is.na(source), "raw", source)) %>% 
  subset(!is.na(rfid)) %>%
  mutate(exp = tolower(exp)) %>% 
  select(cohort, rfid, labanimalid, sex, exp, rewards) %>%  # remove source from data
  spread(exp, rewards)

cocaine_ints_c01_08_merge_pr_shock <- full_join(cocaine_intermediates_xl_pr_corrected, cocaine_intermediates_xl_shock_corrected, by = c("cohort", "rfid", "labanimalid", "sex"))
# 
cohort01_07_phenotypes_merge <- full_join(cohort01_07_cocaine_lga, cohort01_07_cocaine_sha, by = c("cohort", "rfid", "labanimalid", "sex"))

cocaine_phenotypes_merge <- full_join(cohort01_07_phenotypes_merge, cocaine_ints_c01_08_merge_pr_shock, by = c("cohort", "rfid", "sex", "labanimalid"))
cocaine_phenotypes_merge %>% write.csv(file = "cocaine_phenotypes_n223_all_vars.csv", row.names = F) # current csv is written for the 223, not 365
cocaine_phenotypes_merge %>% select(cohort, rfid, sex, labanimalid, esc11_14_mean, mean_sha_last3, starts_with("pr"), shock03) %>% 
  left_join(subjects_exp_age %>% select(cohort, labanimalid, rfid, lga_11_age, sha_08_age, starts_with("pr_"), shock_03_age), 
            by = c("cohort", "labanimalid", "rfid")) %>%  # get age 
  mutate(lga_11_age = replace(lga_11_age, is.na(esc11_14_mean), NA),
         sha_08_age = replace(sha_08_age, is.na(mean_sha_last3), NA),
         pr_01_age = replace(pr_01_age, is.na(pr01), NA),
         pr_02_age = replace(pr_02_age, is.na(pr02), NA),
         pr_03_age = replace(pr_03_age, is.na(pr03), NA),
         shock_03_age = replace(shock_03_age, is.na(shock03), NA)) %>% 
  left_join(box_metadata_wide, by = c("cohort", "labanimalid", "rfid")) %>% 
  mutate(lga_11_box = replace(lga_11_box, is.na(esc11_14_mean), NA),
         sha_08_box = replace(sha_08_box, is.na(mean_sha_last3), NA),
         pr_01_box = replace(pr_01_box, is.na(pr01), NA),
         pr_02_box = replace(pr_02_box, is.na(pr02), NA),
         pr_03_box = replace(pr_03_box, is.na(pr03), NA),
         shock_03_box = replace(shock_03_box, is.na(shock03), NA)) %>% 
  mutate(sex = str_extract(labanimalid, "[MF]")) %>% 
  write.csv(file = "cocaine_phenotypes_n223_stripped_vars.csv", row.names = F)



# will return to cohorts for the expert opinion values 
## once the values have been corrected, now calculate the indices 
cohort1_lga <- cocaine_intermediates_xl_lga_corrected %>% 
  subset(cohort == "C01") %>% 
  spread(exp, rewards) %>% 
  # mutate(lga01_mean = mean(lga_01, na.rm = T), lga01_sd = sd(lga_01, na.rm = T)) %>% 
  # ungroup() %>% 
  mutate_at(vars(matches("lga_1[1-4]$")), list(esc = ~.-lga_01)) %>% 
  mutate(esc11_14_mean = rowMeans(select(., ends_with("_esc")), na.rm = TRUE)) %>% 
  select(cohort, rfid, sex, labanimalid, matches("_esc$"), matches("_mean$"))
  # mutate_at(vars(matches("lga_\\d+$")), list(esc = ~.-lga01_mean)) %>% ## XX SAVE CODE FOR GENERATING THE EXPERT OPINION DATA
  # mutate_at(vars(ends_with("_esc")), ~./lga01_sd) 
# check for row id dupes before calculating means 
# cohort1_lga %>% get_dupes(rfid)
# cohort1_lga <- cohort1_lga %>% 
#   mutate(ind_esc_mean = rowMeans(select(., ends_with("_esc")), na.rm = TRUE)) %>% 
#   group_by(sex) %>% 
#   mutate(fr_zscore = scale(ind_esc_mean)) %>% 
#   ungroup() %>% 
#   mutate_all(as.character)
# plot(density(as.numeric(cohort1_lga$fr_zscore)))

# setwd("~/Dropbox (Palmer Lab)/PalmerLab_Datasets/u01_george_oliviercocaine/database/C01")
# write.csv(cohort1_lga, file = "cohort1_lga.csv", row.names = F)
# writing onto desktop folder bc Dropbox does not work
setwd("~/Desktop/Database/csv files/u01_olivier_george_cocaine/")
write.csv(cohort1_lga, file = "cohort1_lga.csv", row.names = F)



## COHORT 2  
cohort2_lga <- cocaine_intermediates_xl_df_corrected %>% 
  subset(cohort == "C02") %>% 
  spread(exp, rewards) %>% 
  group_by(sex) %>% 
  mutate(lga01_mean = mean(lga_01, na.rm = T), lga01_sd = sd(lga_01, na.rm = T)) %>% 
  ungroup() %>% 
  mutate_at(vars(matches("lga_\\d+$")), list(esc = ~.-lga01_mean)) %>% 
  mutate_at(vars(ends_with("_esc")), ~./lga01_sd) 
# check for row id dupes before calculating means 
cohort2_lga %>% get_dupes(rfid)

cohort2_lga <- cohort2_lga %>% 
  mutate(ind_esc_mean = rowMeans(select(., ends_with("_esc")), na.rm = TRUE)) %>% 
  group_by(sex) %>% 
  mutate(fr_zscore = scale(ind_esc_mean)) %>% 
  ungroup() %>% 
  mutate_all(as.character)
plot(density(as.numeric(cohort2_lga$fr_zscore)))

# setwd("~/Dropbox (Palmer Lab)/PalmerLab_Datasets/u01_george_oliviercocaine/database/C01")
# write.csv(cohort1_lga, file = "cohort1_lga.csv", row.names = F)
# writing onto desktop folder bc Dropbox does not work
setwd("~/Desktop/Database/csv files/u01_olivier_george_cocaine/")
write.csv(cohort2_lga, file = "cohort2_lga.csv", row.names = F)


## COHORT 3  
cohort3_lga <- cocaine_intermediates_xl_df_corrected %>% 
  subset(cohort == "C03") %>% 
  spread(exp, rewards) %>% 
  group_by(sex) %>% 
  mutate(lga01_mean = mean(lga_01, na.rm = T), lga01_sd = sd(lga_01, na.rm = T)) %>% 
  ungroup() %>% 
  mutate_at(vars(matches("lga_\\d+$")), list(esc = ~.-lga01_mean)) %>% 
  mutate_at(vars(ends_with("_esc")), ~./lga01_sd) 
# check for row id dupes before calculating means 
cohort3_lga %>% get_dupes(rfid)
cohort3_lga %>% naniar::vis_miss()

cohort3_lga <- cohort3_lga %>% 
  mutate(ind_esc_mean = rowMeans(select(., ends_with("_esc")), na.rm = TRUE)) %>% 
  group_by(sex) %>% 
  mutate(fr_zscore = scale(ind_esc_mean)) %>% 
  ungroup() %>% 
  mutate_all(as.character)
plot(density(as.numeric(cohort3_lga$fr_zscore)))

# setwd("~/Dropbox (Palmer Lab)/PalmerLab_Datasets/u01_george_oliviercocaine/database/C01")
# write.csv(cohort1_lga, file = "cohort1_lga.csv", row.names = F)
# writing onto desktop folder bc Dropbox does not work
setwd("~/Desktop/Database/csv files/u01_olivier_george_cocaine/")
write.csv(cohort3_lga, file = "cohort3_lga.csv", row.names = F)



## COHORT 4  
## XX uniform variables??? 
cohort4_lga <- cocaine_intermediates_xl_df_corrected %>% 
  subset(cohort == "C04") %>% 
  spread(exp, rewards) %>% 
  group_by(sex) %>% 
  mutate(lga01_mean = mean(lga_01, na.rm = T), lga01_sd = sd(lga_01, na.rm = T)) %>% 
  ungroup() %>% 
  mutate_at(vars(matches("lga_\\d+$")), list(esc = ~.-lga01_mean)) %>% 
  mutate_at(vars(ends_with("_esc")), ~./lga01_sd)
# check for row id dupes before calculating means 
cohort4_lga %>% get_dupes(rfid)
cohort4_lga %>% naniar::vis_miss()

cohort4_lga <- cohort4_lga %>% 
  mutate(ind_esc_mean = rowMeans(select(., matches("1[1-4]_esc$")), na.rm = TRUE)) %>%  ## ONLY USE LAST FOUR LGA VALUES 
  group_by(sex) %>% 
  mutate(fr_zscore = scale(ind_esc_mean)) %>% 
  ungroup() %>% 
  mutate_all(as.character)
plot(density(as.numeric(cohort4_lga$fr_zscore)))

# setwd("~/Dropbox (Palmer Lab)/PalmerLab_Datasets/u01_george_oliviercocaine/database/C01")
# write.csv(cohort1_lga, file = "cohort1_lga.csv", row.names = F)
# writing onto desktop folder bc Dropbox does not work
setwd("~/Desktop/Database/csv files/u01_olivier_george_cocaine/")
write.csv(cohort4_lga, file = "cohort4_lga.csv", row.names = F)


## COHORT 5  
## XX uniform variables??? 
cohort5_lga <- cocaine_intermediates_xl_df_corrected %>% 
  subset(cohort == "C05") %>% 
  spread(exp, rewards) %>% 
  group_by(sex) %>% 
  mutate(lga01_mean = mean(lga_01, na.rm = T), lga01_sd = sd(lga_01, na.rm = T)) %>% 
  ungroup() %>% 
  mutate_at(vars(matches("lga_\\d+$")), list(esc = ~.-lga01_mean)) %>% 
  mutate_at(vars(ends_with("_esc")), ~./lga01_sd) 
# check for row id dupes before calculating means 
cohort5_lga %>% get_dupes(rfid)
cohort5_lga %>% naniar::vis_miss()

cohort5_lga <- cohort5_lga %>% 
  mutate(ind_esc_mean = rowMeans(select(., matches("1[1-5]_esc$")), na.rm = TRUE)) %>% ## ONLY USE LAST THREE LGA VALUES 
  group_by(sex) %>% 
  mutate(fr_zscore = scale(ind_esc_mean)) %>% 
  ungroup() %>% 
  mutate_all(as.character)
plot(density(as.numeric(cohort5_lga$fr_zscore)))

# setwd("~/Dropbox (Palmer Lab)/PalmerLab_Datasets/u01_george_oliviercocaine/database/C01")
# write.csv(cohort1_lga, file = "cohort1_lga.csv", row.names = F)
# writing onto desktop folder bc Dropbox does not work
setwd("~/Desktop/Database/csv files/u01_olivier_george_cocaine/")
write.csv(cohort5_lga, file = "cohort5_lga.csv", row.names = F)


# =================================================
## PR
setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE")
cocaine_intermediates_xl <- u01.importxlsx("Addiction indices for C01-C07 Cocaine.xlsx") %>% 
  lapply(function(x){
    if(grepl("F701", x[1,1])){ # if cohort 7 table
      names(x)[1] <- "rats"
    }
    x <- x %>% clean_names
    i <- 1
    if(names(x)[1] == "x1"){
      names(x) <- x[i,]
      x <- x %>% clean_names
      x <- x[-i, ]
      i = i + 1 # go down the rows until the test expression is false
    }
    
    x <- x %>%
      select(starts_with("rat"), starts_with("pr"), starts_with("day"))
    return(x)
  }) 
cocaine_intermediates_xl$C01<- cocaine_intermediates_xl$C01 %>% select(matches("rats$"), matches("^day_\\d+$"), -matches("^day_\\d+_\\d+"))
cocaine_intermediates_xl$C02<- cocaine_intermediates_xl$C02 %>% select(matches("rats$"), matches("^day_\\d+$"), -matches("^day_\\d+_\\d+"))
cocaine_intermediates_xl$C03<- cocaine_intermediates_xl$C03 %>% select(matches("rats_1$"), matches("^day_\\d_\\d$"), matches("^day_([9]|1\\d)_[1]\\d")) 
cocaine_intermediates_xl$C04<- cocaine_intermediates_xl$C04 %>% select(matches("rat$"), matches("^lg_a\\d+$"), -matches("^lg_a\\d+_2"))
cocaine_intermediates_xl$C05<- cocaine_intermediates_xl$C05 %>% select(matches("rats_1$"), matches("^day_\\d_\\d$"), matches("^day_([9]|1\\d)_[1]\\d")) 
cocaine_intermediates_xl$C07<- cocaine_intermediates_xl$C07 %>% select(matches("rats$"), matches("^lg_a\\d_\\d$"), matches("^lg_a([9]|1\\d)_[1]\\d")) 


cocaine_intermediates_xl_df <- cocaine_intermediates_xl %>% 
  lapply(function(x){
    names(x) <- c("labanimalid", paste0("lga_", str_pad(1:14, "2", "left", "0"))) 
    return(x)
  }) %>% 
  rbindlist(idcol = "cohort", fill = T) %>%  
  mutate(labanimalid = str_extract(labanimalid, "[MF]\\d+")) %>% 
  mutate_at(vars(-matches("cohort|labanimal")), as.numeric) %>% 
  subset(!is.na(labanimalid)) %>% 
  left_join(rat_info_allcohort_xl_df[, c("rfid", "labanimalid")], by = "labanimalid") %>% 
  gather("exp", "rewards_xl", -cohort, -labanimalid, -rfid) %>% 
  left_join(Olivier_Cocaine_df %>% 
              mutate(labanimalid = str_extract(toupper(labanimalid), "[MF]\\d+")) %>% 
              select(cohort, rfid, labanimalid, sex, exp, rewards) %>% distinct() %>% 
              spread(exp, rewards) %>% select(cohort, rfid, labanimalid, sex, matches("LGA(0[1-9]|1[0-4])")) %>% 
              gather("exp", "rewards_raw", -cohort, -rfid, -labanimalid, -sex) %>% 
              mutate(exp = gsub("LGA", "lga_", exp)),
            by = c("cohort", "labanimalid", "rfid", "exp")) %>% 
  mutate(rewards_QC_diff = rewards_xl - rewards_raw,
         rewards_QC = ifelse(rewards_QC_diff == 0, "pass", "fail"))





##### 
## XX ask cocaine team why this, especially for the C09+  
WFU_OlivierCocaine_test_df %>%
  rename("wfu_labanimalid" = "labanimalid") %>%
  mutate(cohort = paste0("C", cohort)) %>%
  dplyr::filter(grepl("^\\d", rfid)) %>% #811 (ignore the blanks and annotations in the excel)
  left_join(., allcohorts2[, c("labanimalid", "rfid")], by = "rfid") %>% subset(is.na(labanimalid))


# tried and troubleshoot TEST %>% subset(labanimalid %in% c("M464", "F516", "M360"))

## Error in `==.default`(date, dplyr::first(na.omit(datedropped))) : 
# comparison (1) is possible only for atomic and list types

left_join(., ratinfo_list_deaths_processed, by = c("rfid", "cohort"))
  
  


sha_rewards_new <- sha_rewards_new %>%
  left_join(., allcohorts2[, c("labanimalid", "rfid")], by = "labanimalid") %>% # add rfid number
  left_join(., ratinfo_list_deaths_processed, by = c("rfid", "cohort")) %>%
  mutate_at(vars(contains("date")), lubridate::ymd) %>%
  group_by(labanimalid) %>%
  mutate(
    comment = case_when(
      !grepl("Died", dplyr::first(na.omit(reasoning)), ignore.case = T) &
        date == dplyr::first(na.omit(datedropped))
      ~ "COMP_EXCLUDE",
      grepl("Died", dplyr::first(na.omit(reasoning)), ignore.case = T) &
        date >= dplyr::first(na.omit(datedropped)) ~ "DEAD_EXCLUDE"
    )
  ) %>%
  ungroup() # %>%
# subset(!is.na(comment))

rewards_sha_tograph <- rewards %>% 
  subset(grepl("SHA", exp)) %>% 
  merge(., allcohorts2 %>% select(labanimalid, rfid, matches("^sha")) %>% distinct() %>% 
                    gather(exp, rewards_excel, sha01:sha10) %>% mutate(exp = toupper(exp)),
                  by = c("labanimalid", "exp")) %>% 
  rename("rewards_raw"= "rewards") %>% 
  mutate(sex = str_extract(labanimalid, "\\D")) %>% 
  left_join(., WFU_OlivierCocaine_test_df[, c("rfid", "dob")], by = "rfid") %>% 
  mutate_at(vars(one_of("date", "dob")), lubridate::ymd) %>% 
  mutate(exp_age = as.integer(difftime(date, dob, unit = "days")))

olivier_sha_measures <- grep("rewards", names(rewards_sha_tograph), value = T) 
rewards_sha_tograph <- rewards_sha_tograph %>% 
  mutate_at(olivier_sha_measures, as.numeric)

# create plots 
setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/QC")

pdf("olivier_sha.pdf", onefile = T)
for (i in 1:(length(olivier_sha_measures)/2)){
  g <-  rewards_sha_tograph %>% 
    mutate(directory = ifelse(grepl("new", directory), "New", "Old")) %>% 
    ggplot(aes_string(x = olivier_sha_measures[i], y = olivier_sha_measures[i+1])) + 
    geom_point(aes(color = directory)) + 
    labs(title = paste0("Excel vs raw data comparison of ", gsub("_(excel|raw)", "", olivier_sha_measures[i]), "\n")) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = 20)) 
  
  # g_cohort <-  ggplot(rewards_sha_tograph, aes_string(x = olivier_sha_measures[i], y = olivier_sha_measures[i+3])) + 
  #   geom_point(aes(color = cohort_number)) + 
  #   facet_grid(~ cohort_number)
  #   labs(title = paste0(olivier_sha_measures[i], "_Raw_VS_Excel_U01_Kalivas", "\n")) + 
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  print(g)
  # print(g_cohort)
}

dev.off()
rewards_sha_tograph %>% dim
rewards_sha_tograph %>% subset(rewards_raw != rewards_excel) %>% dim
rewards_sha_tograph %>% subset(rewards_raw == rewards_excel) %>% dim
rewards_sha_tograph %>% subset(is.na(rewards_raw)|is.na(rewards_excel)) %>% dim
rewards_sha_tograph %>% subset(rewards_raw != rewards_excel) %>% 
  select(labanimalid, exp, filename, rewards_raw, rewards_excel) %>% 
  openxlsx::write.xlsx("sha_compare.xlsx")




# clean the raw files
rewards_sha_df_graph[duplicated(rewards_sha_df_graph[,c("labanimalid", "file_exp")]), ] #should only be one data point for each animal but we are seeing multiple

rewards_sha_df_graph <- rewards_sha_df %>% 
  dplyr::filter(bin == "total") %>% 
  merge(., rewards_sha_df %>% 
          dplyr::filter(bin != "total") %>% 
          group_by(labanimalid, file_exp) %>% 
          summarize(tot_counts_calc = sum(counts))) # merge df of sum of existing numbers to calculated sums

rewards_sha_df_graph %>% dplyr::filter(tot_counts_calc != counts) %>% 
  select(labanimalid, file_cohort, counts, tot_counts_calc, filename, file_exp)# show cases that these values don't match

##################
## LGA ########### 
##################
setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/QC")

# lga NEW and OLD combine 
lga_rewards <- rbindlist(list("lga_rewards_new" = lga_rewards_new, 
                              "lga_rewards_old" = lga_rewards_old), idcol = "directory", fill = T) %>% select(-filename1)
rewards_lga_tograph <- lga_rewards %>% merge(., allcohorts2 %>% select(labanimalid, rfid, matches("^lga")) %>% distinct() %>% 
                                               gather(exp, rewards_excel, lga01:lga23) %>% mutate(exp = toupper(exp)),
                                             by = c("labanimalid", "exp")) %>% 
  rename("rewards_raw"= "rewards") %>% 
  mutate(sex = str_extract(labanimalid, "\\D")) %>% 
  left_join(., WFU_OlivierCocaine_test_df[, c("rfid", "dob")], by = "rfid") %>% 
  mutate_at(vars(one_of("date", "dob")), lubridate::ymd) %>% 
  mutate(exp_age = as.integer(difftime(date, dob, unit = "days")))

olivier_lga_measures <- grep("rewards", names(rewards_lga_tograph), value = T) 
rewards_lga_tograph <- rewards_lga_tograph %>% 
  mutate_at(olivier_lga_measures, as.numeric)

# NOT IN DATABASE 
# > rewards_lga_tograph %>% subset(is.na(exp_age)) %>% select(labanimalid) %>% unique
labanimalid
2250        F724 (died and replace)
2907        F830 933000320047252
2923        F837 933000320047268


# create plots 
pdf("olivier_lga.pdf", onefile = T)
for (i in 1:(length(olivier_lga_measures)/2)){
  g <-  ggplot(rewards_lga_tograph, aes_string(x = olivier_lga_measures[i], y = olivier_lga_measures[i+1])) + 
    geom_point(aes(color = directory), size = 0.1) + 
    facet_grid(~ cohort) +
    labs(title = paste0(olivier_lga_measures[i], "_Raw_VS_Excel_U01_Olivier", "\n")) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  # g_cohort <-  ggplot(rewards_lga_tograph, aes_string(x = olivier_lga_measures[i], y = olivier_lga_measures[i+3])) + 
  #   geom_point(aes(color = cohort_number)) + 
  #   facet_grid(~ cohort_number)
  #   labs(title = paste0(olivier_lga_measures[i], "_Raw_VS_Excel_U01_Kalivas", "\n")) + 
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  print(g)
  # print(g_cohort)
}
dev.off()

rewards_lga_tograph %>% dplyr::filter(rewards_raw != rewards_excel) %>% 
  select(labanimalid, exp, cohort, directory, filename, rewards_raw, rewards_excel)%>% 
  openxlsx::write.xlsx("lga_compare.xlsx")

rewards_lga_tograph %>% subset(rewards_raw != rewards_excel) %>% dim

