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


selfadmin_df_gwascopy <- selfadmin_df %>% 
  mutate(avgsha7_10 = rowMeans(.[, sha07:sha10]),
         avglga11_14 = rowMeans(.[, lga11:lga14]),
         ) %>% # rowMeans(a[c('high', 'low')], na.rm=TRUE) XX PICK UP HERE -- NOT SURE IF NA SHOULD BE REMOVED FROM MEAN CALCULATIONS
  mutate(`acquis_sha` = ifelse(avgsha7_10 < 4, "No", "Yes"),
         `escalation` = ifelse(avglga11_14 > 4, "Vulnerable", "Resistant"))
# `add ind` =IF(BW2>=0.49,"High AI","Low AI")
# `fr` =IF(BT2>=0.49,"High FR","Low FR")
# `pr` =IF(BU2>=0.49,"High PR","Low PR")
# `shock`=IF(BV3>=0.49,"Res Shock","Vuln Shock")

rowMeans(a[c('high', 'low')], na.rm=TRUE)

## TO DO:  
# check if spleen and ceca are in the naive
# extract naive ids from wfu 
# extract ids from the spleen and ceca shipment files 
# check if matches 

setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/20190829_WFU_U01_ShippingMaster")
# extract the naive cases from the wfu data
# 