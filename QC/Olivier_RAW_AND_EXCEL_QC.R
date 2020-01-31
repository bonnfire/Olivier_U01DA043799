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

# uncomment for sha_rewards_old
# sha_rewards_new %>% 
#   mutate_at(vars(exp), as.factor) %>% 
#   mutate(sex = str_extract(labanimalid, "[M|F]")) %>% 
#   ggplot(aes(exp, rewards, group = labanimalid, color = sex)) + 
#   geom_line() +
#   facet_grid( ~ cohort) + 
#   labs(title = "SHA Rewards New (Raw Only) Directories, For Each Rat") + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 




# sha NEW and OLD combine 
sha_rewards <- rbindlist(list("sha_rewards_new" = sha_rewards_new, 
                              "sha_rewards_old" = sha_rewards_old), idcol = "directory", fill = T)
rewards_sha_tograph <- sha_rewards %>% merge(., allcohorts2 %>% select(labanimalid, rfid, matches("^sha")) %>% distinct() %>% 
                    gather(exp, rewards_excel, sha01:sha10) %>% mutate(exp = toupper(exp)),
                  by = c("labanimalid", "exp")) %>% 
  rename("rewards_raw"= "rewards") %>% 
  left_join(., WFU_OlivierCocaine_test_df[, c("rfid", "dob")], by = "rfid")

# add notes about missingness (file or dead)
sha_rewards_new <- sha_rewards_new %>%
  left_join(., allcohorts2[, c("labanimalid", "rfid")], by = "labanimalid") %>%
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

ratinfo_list_deaths_processed %>% dplyr::filter(grepl("Died",reasoning)) %>% dim ## the total number of dead animals should match here

olivier_sha_measures <- grep("rewards", names(rewards_sha_tograph), value = T) 
rewards_sha_tograph <- rewards_sha_tograph %>% 
  mutate_at(olivier_sha_measures, as.numeric)

# create plots 
pdf("olivier_sha.pdf", onefile = T)
for (i in 1:(length(olivier_sha_measures)/2)){
  g <-  ggplot(rewards_sha_tograph, aes_string(x = olivier_sha_measures[i], y = olivier_sha_measures[i+1])) + 
    geom_point(aes(color = directory)) + 
    labs(title = paste0(olivier_sha_measures[i], "_Raw_VS_Excel_U01_Olivier", "\n")) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
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

