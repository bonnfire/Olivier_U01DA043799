## CREATE DATABASE STRUCTURE VARIABLES 

library(data.table)
library(dplyr)

# exp, varnames1, varnames2, 
list <- list("cohort1" = cohort1$tselfadmin, 
             "cohort2" = cohort2$tselfadmin,
             "cohort3" = cohort3$tselfadmin, 
             "cohort4" = cohort4$tselfadmin, 
             "cohort5" = cohort5$tselfadmin,
             "cohort7" = cohort7$tselfadmin,
             "cohort8" = cohort8$tselfadmin) # cohort 5 isn't working 


# listwnames <- lapply(list, function(x){
#   names(x) <- c(paste0("Var", 1:length(x)))
#   return(x)
# })

# allcohorts <- bind_rows(list, .id = "cohort") # creates blank inaccessible columns
allcohorts2 <- rbindlist(list, idcol = "cohort", fill= T) %>% 
  select(-date_labanimalid)

# make dataframe 
test <- names(allcohorts)
varnametally <- sapply(list, function(x){
  ifelse(test %in% names(x), 1, 0) 
  }) %>% t() %>% as.data.frame()
names(varnametally) <- test
varnametally$cohort <- rownames(varnametally)
rownames(varnametally) <- NULL

save(varnametally, file = "Oliviercohortvariables.RData")


test_measurevars <- grep("^(?!rf|date|comment|lab|cohort)", test, value = T, perl = T)
varnametally_measurevars <- sapply(list, function(x){
  measurevars <- grep("^(?!rf|date|comment|lab|cohort)", names(x), value = T, perl = T) 
  ifelse(test_measurevars %in% measurevars, 1, 0) 
}) %>% t() %>% as.data.frame()
names(varnametally_measurevars) <- test_measurevars
varnametally_measurevars <- cbind(cohort = rownames(varnametally_measurevars), varnametally_measurevars)
rownames(varnametally_measurevars) <- NULL

save(varnametally, file = "Oliviercohortvariables.RData")

