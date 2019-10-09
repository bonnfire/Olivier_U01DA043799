# Make function for Olivier's data (self admin)

library(dplyr)
library(data.table)
library(tidyverse)
library(readxl)

olivierfiles <- function(filename){
  options(scipen = 100) # prevent sci notation
  u01.importxlsx <- function(xlname){
    df <- lapply(excel_sheets(path = xlname), read_excel, path = xlname)
    names(df) <- excel_sheets(xlname)
    return(df)
  }
  
  selfadmin <- u01.importxlsx(filename)[[1]] %>%
    as.data.table %>% 
    na.omit(cols = seq_along('Rat'))
  
  # turn all into character for preservation in transpose
  selfadmin[, names(selfadmin) := lapply(.SD, as.character)] 
  
  # extract experiment name (ShA and LgA)
  selfadmin$Rat <- ifelse(grepl("^\\D{2}A", selfadmin$Rat), as.character(stringr::str_match(selfadmin$Rat,"..A")), selfadmin$Rat)
  
  # make experiment name unique (# code from G. Grothendieck (Stack Overflow) )
  uniquify <- function(x) if (length(x) == 1) x else sprintf("%s%02d", x, seq_along(x)) 
  selfadmin$Rat <- ave(selfadmin$Rat, selfadmin$Rat, FUN = uniquify) 
  
  ### transpose data 
  tselfadmin <- data.table::transpose(selfadmin)
  colnames(tselfadmin) <- as.character(tselfadmin[1,])
  tselfadmin <- tselfadmin[-1,] # data.table cannot delete rows by reference with tselfadmin[1 := NULL,]  
  
  # variable name clean
  setnames(tselfadmin, gsub(" ", "_", tolower(names(tselfadmin))))
  
  # add exp date columns
  nm <- names(tselfadmin)[-1] # use these columns to make date columns
  nm1 <- paste("date", nm, sep = "_") # make these date columns
  tselfadmin[ , ( nm1 ) := lapply( .SD, function(x) c(grep("^\\d{4}\\-(0?[1-9]|1[012])\\-(0?[1-9]|[12][0-9]|3[01])$", x, value = T)) ) ,  .SDcols = nm ]
  ind <- grep("^(date)", names(tselfadmin), perl = T) # bring dates back to posixct
  for (i in seq_along(ind)) {
    set(tselfadmin, NULL, ind[i], as.POSIXct(tselfadmin[[ind[i]]], tz = "UTC"))
  }
  tselfadmin <- tselfadmin[!grepl("^\\d{4}\\-", sha01),] # remove the row with dates
  
  # add general comments columns
  rownumber <- which(is.na(tselfadmin$rfid)) %>% head(1)
  colswithcomments <- colnames(tselfadmin)[which(!is.na(tselfadmin[rownumber,]))] # return columns that have comments
  colswithcomments <- grep("^(?!date)", colswithcomments, perl = T, value = T) # return non-date columns that have comments
  nm2 <- paste("comment", colswithcomments, sep = "_")
  tselfadmin[ , ( nm2 ) := lapply( .SD, function(x) c(grep("^.", x[rownumber], value = T))) ,  .SDcols = colswithcomments ]
  tselfadmin <- tselfadmin[-rownumber,] # remove the row with general comments
  
  # extract [*DATA DICTIONARY*]
  datadictionary <- selfadmin[, 1:3] # vertical formatting is preferred in selfadmin
  tselfadmin <- tselfadmin[-c(1:2), ] # remove data dictionary from data
  
  # extract [*TABLE FOR SPECIFIC COMMENTS*]
  rownumbers <- which(is.na(tselfadmin$rfid))
  colswithspeccomments <- colnames(tselfadmin)[which(!is.na(tselfadmin[rownumbers[1],]))] # return columns that have comments
  colswithspeccomments <- grep("^(?![date|comment])", colswithspeccomments, perl = T, value = T) # further prune columns that don't contain comment or date
  specificcomments <- tselfadmin[rownumbers, ..colswithspeccomments] # extract special comments
  tselfadmin <- tselfadmin[-rownumbers, -..colswithspeccomments] # remove the special comments
  
  # add cohort column
  tselfadmin <- append(tselfadmin, list(data.table(cohort = sub("_cocaine.*", "", filename)))) %>% as.data.table   # get from file name
  
  # add rat id column
  tselfadmin <- append(tselfadmin, list(data.table(labanimalid = grep("\\D\\d", names(selfadmin), value = T)))) %>% as.data.table
  
  # ensure all na's are properly notated
  # ind <- grep(pattern = "^[[:alpha:]]{2,3}[[:digit:]]{2}$", names(tselfadmin), perl = T)
  ind <- grep(pattern = "^(?!rf|date|comment|lab|cohort)", names(tselfadmin), perl = T)
  for(j in seq_along(ind)){
    set(tselfadmin, i=which(tselfadmin[[j]]=="n/a"), j=j, value=NA)
  }

  # quant data should be numeric
  for (i in seq_along(ind)) {
    set(tselfadmin, NULL, ind[i], as.numeric(tselfadmin[[ind[i]]]))
  }

  list <- list("tselfadmin" = tselfadmin,
               "datadictionary" = datadictionary,
               "specificcomments" = specificcomments)
  
  return(list)
} 


cohortfiles <- list.files()
for(i in 1:length(cohortfiles)){
  if(i == 4) next 
  assign(paste0("cohort", i), olivierfiles(cohortfiles[i]))
}

cohort4 <-  olivierfiles("C04_ cocaine.xlsx") # 4 breaks, 

#
cohort5ma <- u01.importxlsx("C05_cocaine.xlsx")[[1]]

selfadmin <- cohort4ma %>%
  as.data.table %>% 
  na.omit(cols = seq_along('Rat'))


selfadminRat <- selfadmin$Rat
selfadminRat <- c(selfadminRat, "Date_ShA01", "Date_ShA02", "Comment_ShA01", "LabanimalID", "Cohort")
ind <- grep(pattern = "^(?!RFID | Date | Comment | LabanimalID | Cohort)", selfadminRat, perl = T, value = T)

