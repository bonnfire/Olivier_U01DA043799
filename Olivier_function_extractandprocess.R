# Make function for Olivier's data (self admin)

library(dplyr)
library(data.table)
library(tidyverse)
library(readxl)

install.packages('splitstackshape')
install.packages('janitor')
library(splitstackshape)
library(janitor)

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
  if(any(grepl("^[[:digit:]]{5,}", selfadmin$Date))){
    selfadmin$Date <- as.character(as.POSIXct(as.numeric(selfadmin$Date) * (60*60*24), origin="1899-12-30", tz="UTC", format="%Y-%m-%d"))
  } # convert Excel character into dates
  
  # fix missing column names
  if(any(grepl("REWARDS|RAT" , names(selfadmin)))){
    setnames(selfadmin, old = c("REWARDS", "RAT"), new = c("Rat", "Rat"), skip_absent=TRUE) # prepares for next code
  }
  names(selfadmin)[length(names(selfadmin))] = "FLAG" # the last column should be FLAG
  
  # extract experiment name (ShA and LgA)
  selfadmin$Rat <- ifelse(grepl("^\\D{2}A", selfadmin$Rat), as.character(stringr::str_match(selfadmin$Rat,"^\\D{2}A")), selfadmin$Rat) # reg exp fixes issuse of extracting mA
  selfadmin$Rat <- gsub(" |[(]|[)]|-", "", selfadmin$Rat) # remove all unwanted characters
  
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
  ind <- grep("^(date)", names(tselfadmin), perl = T) # bring date characters back to posixct
  for (i in seq_along(ind)) {
    set(tselfadmin, NULL, ind[i], as.POSIXct(tselfadmin[[ind[i]]], tz = "UTC"))
  }
  tselfadmin <- tselfadmin[!grepl("^\\d{4}\\-", sha01),] # remove the row with dates
  
  # extract [*DATA DICTIONARY*]
  datadictionary <- selfadmin[, 1:3] # vertical formatting is preferred in selfadmin
  tselfadmin <- tselfadmin[-c(1:2), ] # remove data dictionary from data
  
  # add general comments columns
  # remove rows under rfid that don't have comments
  # rownumber <- which(is.na(tselfadmin$rfid)) %>% head(1) 
  
  # remove <- which(grepl("^[[:alpha:]]", tselfadmin$rfid)) %>% head(1) # define remove cases
  #   if(length(remove) != 0){
  #   tselfadmin <- tselfadmin[-remove, ]
  #   }
  
  rownumber <- which(is.na(tselfadmin$rfid)|grepl("^[^0-9]", tselfadmin$rfid)) %>% head(1)
  colswithcomments <- colnames(tselfadmin)[which(!is.na(tselfadmin[rownumber, ]))] # return columns that have comments
  colswithcomments <- grep("^(?!date)", colswithcomments, perl = T, value = T) # return non-date columns that have comments
  nm2 <- paste("comment", colswithcomments, sep = "_")
  tselfadmin[ , ( nm2 ) := lapply( .SD, function(x) c(grep("^.", x[rownumber], value = T))) ,  .SDcols = colswithcomments ]
  tselfadmin <- tselfadmin[-rownumber,] # remove the row with general comments

  
  # extract [*TABLE FOR SPECIFIC COMMENTS*]
  rownumbers <- which(is.na(tselfadmin$rfid))
  colswithspeccomments <- colnames(tselfadmin)[which(!is.na(tselfadmin[rownumbers[1],]))] # return columns that have comments
  colswithspeccomments <- grep("^(?![date|comment])", colswithspeccomments, perl = T, value = T) # further prune columns that don't contain comment or date
  specificcomments <- tselfadmin[rownumbers, ..colswithspeccomments] # extract special comments
  specificcommentsdates <- sapply(colswithspeccomments, grep, names(tselfadmin), value = T)
  specificcommentsdates <- sapply(specificcommentsdates, "[", 2) %>% as.character # extract the date column name from every item in list
  specificcommentsdates <- tselfadmin[1, ..specificcommentsdates] ## get dates
  tselfadmin <- tselfadmin[-rownumbers, ] # remove the special comments
  
  # reshape TABLE FOR SPECIFIC COMMENTS
  specificcomments <- rbindlist(list(specificcomments, as.list(names(specificcomments))), fill=FALSE) # preserve experiment name as rows before tranposing
  tspecificcomments <- data.table::transpose(specificcomments)
  tspecificcomments <- janitor::remove_empty(tspecificcomments, "cols") # remove empty columns
  specificcommentsdateslg <- melt(specificcommentsdates, variable.name = "experiment", value.name = "date")  # prepare dates # melt from wide to long 
  specificcommentsdateslg$experiment <- gsub("date_", "", specificcommentsdateslg$experiment)
  setnames(tspecificcomments, c("conflict", "comment", "flag", "experiment"))
  setkey(tspecificcomments, cols = experiment) # set key to prepare merge 
  tspecificcomments <- tspecificcomments[specificcommentsdateslg] # specific comments now has date info 
  tspecificcomments <- splitstackshape::cSplit(tspecificcomments, splitCols = "conflict", sep = ";", direction = "long") # breaks multiple cases in one entry into numerous rows
  tspecificcomments[, c("labanimalid", "conflict") := tstrsplit(conflict, ":", fixed=TRUE)] # data table splits text string
  tspecificcomments[, ("labanimalid") := sapply(.SD, function(x) c(stringr::str_extract(x, "[[:alnum:]]+[[:digit:]]+$"))) , .SDcols = "labanimalid"] # clean data
  tspecificcomments[, conflict := str_trim(conflict)] # should remove trailing whitespace
  
  # extract [*TABLE FOR COMPROMISED ANIMALS*]
  
  # extract [*TABLE FOR DEAD ANIMALS*]
  dead <- tspecificcomments[conflict %like% "(?i)died|dead"] # data.table %ilike% operator doesn't seem ready so using (?i) tag for case insensitive searches
  
  # extract [*TABLE FOR SWITCHED ANIMALS*]
  switches <- tspecificcomments[conflict %like% "(?i)switch(ed)?"]
  
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
               "specificcomments" = tspecificcomments,
               "dead" = dead,
               "switches" = switches)
  
  return(list)
} 


cohortfiles <- list.files(pattern = "*.xlsx")
# cohort 6 is missing
for(i in 1:length(cohortfiles)){
  # if(i == 4) next 
  if(i > 5){
    assign(paste0("cohort", i + 1), olivierfiles(cohortfiles[i]))
  } else
    assign(paste0("cohort", i), olivierfiles(cohortfiles[i]))
}



# function breaks at 
cohort2 <- olivierfiles("C02_cocaine.xlsx")
cohort3 <- olivierfiles("C03_cocaine.xlsx")
cohort4 <- olivierfiles("C04_cocaine.xlsx") # specific comment in non specific column 
cohort5 <- olivierfiles("C05_cocaine.xlsx") # note about renumbering -- can we ignore
cohort7 <- olivierfiles("C07_cocaine.xlsx")
cohort8 <- olivierfiles("C08_cocaine.xlsx") # comments aren't working because there is a row that needs to be removed
