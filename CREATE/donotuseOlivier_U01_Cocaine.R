setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/Olivier_George_U01DA043799 (Cocaine)")

Olivier_U01_Cocaine_original <- read_excel()

# QC: sex vs shipment data sex  
## process WFU data 
WFU_Olivier_co_sexQC <- WFU_Olivier_co
WFU_Olivier_co_sexQC[[1]] <- data.frame(rfid = WFU_Olivier_co[[1]]$`Transponder ID` %>% as.character(),
                                        sex = WFU_Olivier_co[[1]]$`Sex`)
WFU_Olivier_co_sexQC[[2]] <- data.frame(rfid = WFU_Olivier_co[[2]]$`Transponder ID` %>% as.character(),
                                        sex = WFU_Olivier_co[[2]]$`Sex`)
WFU_Olivier_co_sexQC[[3]] <- data.frame(rfid = WFU_Olivier_co[[3]]$`Transponder ID` %>% as.character(),
                                        sex = WFU_Olivier_co[[3]]$`Sex`)
WFU_Olivier_co_sexQC[[4]] <- data.frame(rfid = WFU_Olivier_co[[4]]$`Transponder ID` %>% as.character(),
                                        sex = WFU_Olivier_co[[4]]$`Sex`)
WFU_Olivier_co_sexQC[[5]] <- data.frame(rfid = WFU_Olivier_co[[5]][6],
                                        sex = WFU_Olivier_co[[5]][5])
colnames(WFU_Olivier_co_sexQC[[5]]) <- c("rfid", "sex")
WFU_Olivier_co_sexQC[[5]] <- WFU_Olivier_co_sexQC[[5]][grep("^9", WFU_Olivier_co_sexQC[[5]]$rfid), ]
WFU_Olivier_co_sexQC[[6]] <- data.frame(rfid = WFU_Olivier_co[[6]][6],
                                        sex = WFU_Olivier_co[[6]][5])
colnames(WFU_Olivier_co_sexQC[[6]]) <- c("rfid", "sex")
WFU_Olivier_co_sexQC[[6]] <- WFU_Olivier_co_sexQC[[6]][grep("^9", WFU_Olivier_co_sexQC[[6]]$rfid), ]
WFU_Olivier_co_sexQC[[7]] <- data.frame(rfid = WFU_Olivier_co[[7]][6],
                                        sex = WFU_Olivier_co[[7]][5])
colnames(WFU_Olivier_co_sexQC[[7]]) <- c("rfid", "sex")
WFU_Olivier_co_sexQC[[7]] <- WFU_Olivier_co_sexQC[[7]][grep("^9", WFU_Olivier_co_sexQC[[7]]$rfid), ]
WFU_Olivier_co_sexQC[[8]] <- data.frame(rfid = WFU_Olivier_co[[8]][6],
                                        sex = WFU_Olivier_co[[8]][5])
colnames(WFU_Olivier_co_sexQC[[8]]) <- c("rfid", "sex")
WFU_Olivier_co_sexQC[[8]] <- WFU_Olivier_co_sexQC[[8]][grep("^9", WFU_Olivier_co_sexQC[[8]]$rfid), ]
WFU_Olivier_co_sexQC[[9]] <- data.frame(rfid = WFU_Olivier_co[[9]][6],
                                        sex = WFU_Olivier_co[[9]][5])
colnames(WFU_Olivier_co_sexQC[[9]]) <- c("rfid", "sex")
WFU_Olivier_co_sexQC[[9]] <- WFU_Olivier_co_sexQC[[9]][grep("^9", WFU_Olivier_co_sexQC[[9]]$rfid), ]

WFU_Olivier_co_sexQC <- bind_rows(WFU_Olivier_co_sexQC, .id = "shipmentdate")

## get this function to work
WFU_Olivier_co_sexQC[[7]] <- WFU_Olivier_co_sexQC[[7]][grep("^9", WFU_Olivier_co_sexQC[[7]]$rfid), ]
WFU_Olivier_co_sexQC[5:9] <- lapply(WFU_Olivier_co[5:9], function(x) x[-1])

colnames(WFU_Olivier_co_sexQC[5:9]) <- lapply(WFU_Olivier_co[5:9], function(x) x[1])

## process Olivier data 

library(dplyr)
library(data.table)
library(tidyverse)
library(readxl)
setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/Olivier_George_U01DA043799 (Cocaine)/Olivier_George_U01/DATA Updated")

u01.importxlsx2 <- function(xlname){
  df <- lapply(excel_sheets(path = xlname), read_excel, path = xlname)
  names(df) <- excel_sheets(xlname)
  return(df)
}
options(scipen = 100) # fixes scientific notation

filename <- "C01_cocaine.xlsx"
test <- u01.importxlsx2(filename)[[1]] %>%
  as.data.table %>% 
  na.omit(cols = seq_along('Rat')) 
# ind <- grep("^(?!Date)", names(test), perl = T) # actually keep dates in character form
# for (i in seq_along(ind)) {
#   set(test, NULL, ind[i], as.character(test[[ind[i]]]))
# }
test[, names(test) := lapply(.SD, as.character)] # turn all into character for preservation in transpose
# sapply(test, class)
# str(test) # reformat data types # NULL represents applying to all rows
test$Rat <- ifelse(grepl("..A", test$Rat), as.character(stringr::str_match(test$Rat,"..A")), test$Rat) # extract experiment name (ShA and LgA)
# test <- separate(test, "Rat", c("Rat", NA), sep = "-[[:space:]]") Code above preferred because of the inconsistent format
uniquify <- function(x) if (length(x) == 1) x else sprintf("%s%02d", x, seq_along(x)) # %s placeholder for post-specified vector and %02d integer values in 2 or more digits, the first being zero if # is <10
test$Rat <- ave(test$Rat, test$Rat, FUN = uniquify) # code from G. Grothendieck (Stack Overflow) # make experiment name unique

transposetest <- data.table::transpose(test)
colnames(transposetest) <- as.character(transposetest[1,])
transposetest <- transposetest[-1,] # data.table cannot delete rows by reference with transposetest[1 := NULL,]  

# add exp date columns 
nm <- names(transposetest)[-1] # use these columns to make date columns
nm1 <- paste("Date", nm, sep = "_") # make these date columns
transposetest[ , ( nm1 ) := lapply( .SD, function(x) c(grep("^\\d{4}\\-(0?[1-9]|1[012])\\-(0?[1-9]|[12][0-9]|3[01])$", x, value = T)) ) ,  .SDcols = nm ]
ind <- grep("^(Date)", names(transposetest), perl = T) # bring dates back to posixct
for (i in seq_along(ind)) {
  set(transposetest, NULL, ind[i], as.POSIXct(transposetest[[ind[i]]], tz = "UTC"))
}
transposetest <- transposetest[!grepl("^\\d{4}\\-", ShA01),] # remove the row with dates

# add general comments columns
rownumber <- which(is.na(transposetest$RFID)) %>% head(1)
colswithcomments <- colnames(transposetest)[which(!is.na(transposetest[rownumber,]))] # return columns that have comments
colswithcomments <- grep("^(?!Date)", colswithcomments, perl = T, value = T) # return non-date columns that have comments
nm2 <- paste("Comment", colswithcomments, sep = "_")
transposetest[ , ( nm2 ) := lapply( .SD, function(x) c(grep("^.", x[rownumber], value = T))) ,  .SDcols = colswithcomments ]
transposetest <- transposetest[-rownumber,] # remove the row with general comments

# extract [*DATA DICTIONARY*]
datadictionary <- test[, 1:3]
colnames(datadictionary) <- c("var_abbv", "var_type", "var_graphtext")
datadictionary$var_description <- NA
uniquify_graphtext <- function(x) if (length(x) == 1) x else paste0("day_", sprintf("%s%02d", x, seq_along(x))) # don't want this 
datadictionary$Rat <- ave(datadictionary$var_graphtext, datadictionary$var_graphtext, FUN = uniquify_graphtext)

transposetest <- transposetest[-c(1:2), ]
# set(test, 1:3, NULL) # remove data dictionary from data

# extract [*TABLE FOR SPECIFIC COMMENTS*]
rownumbers <- which(is.na(transposetest$RFID))
colswithspeccomments <- colnames(transposetest)[which(!is.na(transposetest[rownumbers[1],]))] # return columns that have comments
colswithspeccomments <- grep("^(?![Date|Comment])", colswithspeccomments, perl = T, value = T) 
specificcomments <- transposetest[rownumbers, ..colswithspeccomments] # extract special comments 
transposetest <- transposetest[-rownumbers, -..colswithspeccomments] # remove the special comments

# test if comments made it in 
comments <- grep("Comment", names(transposetest))
transposetest[, ..comments] ## checked, the comments all made it in

# add cohort column
# get from file name 
transposetest <- append(transposetest, list(data.table(Cohort = sub("_cocaine.*", "", filename)))) %>% as.data.table

# add rat id column
transposetest <- append(transposetest, list(data.table(LabAnimalID = grep("\\D\\d", names(test), value = T)))) %>% as.data.table

# ensure all na's are properly notated
for(j in seq_along(transposetest)){
  set(transposetest, i=which(transposetest[[j]]=="n/a"), j=j, value=NA)
}

# quant data should be numeric
ind <- grep("^(?![RFID|Date|Comment])", transposetest, perl = T, value = T)
for (i in seq_along(ind)) {
  set(transposetest, NULL, ind[i], as.numeric(transposetest[[ind[i]]]))
}

# variable name clean
setnames(transposetest, gsub(" ", "_", tolower(names(transposetest))))



# todo: 
# ask about series of 0 as na? 

# later use: create a function that returns various things
foo <- 12
bar <- c("a", "b", "e")
newList <- list("integer" = foo, "names" = bar)

########### pick up to create a list of datatables 

setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/Olivier_George_U01DA043799 (Cocaine)/Olivier_George_U01/DATA Updated")
path <- getwd()

mylist.names <- rep(paste0("cohort", 1:8), each = 2) 
mylist <- sapply(mylist.names, function(x) NULL)

# u01.importxlsx.olivier <- function(path){
#   xlname <- list.files(path = path, pattern = "*.xlsx", full.names = F)
#   lapply(seq_along(mylist, function(i){
#     mylist[[i]] <- lapply(xlname, function(x){
#       read_excel(path = x)
#       path_sheetnames <- excel_sheets(x)
#     })
#   }))
# }
# Olivier_co_all <- u01.importxlsx.olivier(path)

u01.importxlsx <- function(xlname){
  df <- lapply(excel_sheets(path = xlname), read_excel, path = xlname, col_names = FALSE)
  names(df) <- excel_sheets(xlname)
  return(df)
}  # diff from u01 function because of col_names and need to transpose for observations each row


xlname <- list.files(path = path, pattern = "*.xlsx", full.names = F)
for(i in 1:length(mylist)){
  for(j in 1:length(xlname)){
    data <- u01.importxlsx(xlname[j]) %>% 
      as.data.table
  }
  mylist[[i]] <- data[[1]]
  mylist[[i + 1]] <- data[[2]]
  names(mylist) <- gsub(" ", "", paste0(names(mylist), sapply(xlname, excel_sheets)))
}

# transpose data tables
dcast(melt(mydata, id.vars = "col0"), variable ~ col0)

sapply(xlname, excel_sheets)

# testing the setup for data tables
empty_list <- list()
dt <- u01.importxlsx("C01_cocaine.xlsx")[[1]] %>% 
  drop_na() %>%
  as.data.table()
empty_list <- append(empty_list, list(dt))
dt2 <- u01.importxlsx("C01_cocaine.xlsx")[[2]] %>% as.data.table()
empty_list <- append(empty_list, list(dt2))
empty_list[2]

