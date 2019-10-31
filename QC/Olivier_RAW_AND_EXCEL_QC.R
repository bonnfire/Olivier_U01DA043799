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
         add_ind_cat = ifelse(addind>0.49|addind==0.49, "High AI", "Low AI")
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
##