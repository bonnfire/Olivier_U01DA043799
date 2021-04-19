## create csv 

## wfu metadata 
setwd("~/Desktop/Database/csv files/u01_olivier_george_cocaine")
write.csv(rats_allcohorts_weights, "rats_allcohorts_weights_c01_11.csv", row.names = F)

write.csv(oliviercocaine_14_wfu_metadata, "~/Desktop/Database/csv files/u01_olivier_george_cocaine/oliviercocaine_14_wfu_metadata.csv", row.names = F)

## 

read.csv("~/Downloads/db_cocaine_selfadmin_phenotypes.csv", colClasses = "character") %>% 
  distinct(cohort) 


