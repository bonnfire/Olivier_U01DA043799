## CREATE DATABASE STRUCTURE VARIABLES 

library(data.table)

exp, varnames1, varnames2, 
list <- list("cohort1" = names(cohort1, 
             "cohort2" = names(cohort2,
             "cohort3" = cohort3, 
             "cohort4" = cohort4, 
             "cohort5" = cohort5)
allcohorts <- rbindlist(list, idcol = "cohort", cill = T)
