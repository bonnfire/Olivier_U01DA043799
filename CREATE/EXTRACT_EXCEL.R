# Make function for Olivier's data (self admin)

library(dplyr)
library(data.table)
library(tidyverse)
library(readxl)

#install.packages('splitstackshape')
#install.packages('janitor')
library(splitstackshape)
library(janitor)
library(stringr)


# incorporate if statements into function 

# from Giordano's email 10/16/19
# LgA days varies bc "we want to keep running our rats before dissection and between treatments (PR or Shock)....[weekends or vacations  -> different LgA every cohort]." 
# Final analysis and addiction index = first 14 days of long access
# Starting Cohort 7 and 8, dropped two shock sessions (0.1 and 0.2) because the animals were responding normally to these; so you will only see 1 intensity
# Cohort 8 doesn't have a 1h session 

olivierfiles <- function(filename){
  # setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/Olivier_George_U01DA043799 (Cocaine)/Olivier_George_U01/DATA Updated")
  setwd("~/Dropbox (Palmer Lab)/Olivier_George_U01/DATA Updated")
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
  selfadmin$Rat <- ifelse(grepl("^PR", selfadmin$Rat), as.character(stringr::str_match(selfadmin$Rat,"PR")), selfadmin$Rat)
  selfadmin$Rat <- ifelse(grepl("PreShock[[:space:]][(LgA15)]", selfadmin$Rat), "LgA", selfadmin$Rat) # applies to cohort 8; technician error 
  selfadmin$Rat <- gsub(" |[(]|[)]|[-]", "", selfadmin$Rat) # remove all unwanted characters
  selfadmin$Rat <- gsub("^(1hr|PreShock).+", "preshock", selfadmin$Rat, ignore.case = T) # XX GET HIS INPUT ON THIS -- currently one HOUR admin; is preshock (in cohort7 also the same)
  selfadmin$Rat <- ifelse(grepl("^Shock", selfadmin$Rat, ignore.case = T), "Shock", selfadmin$Rat)
  # selfadmin$Rat <- ifelse(grepl("Shock", selfadmin$Rat, ignore.case = T), str_extract(selfadmin$Rat, "[0][.][1-3]"), selfadmin$Rat) # make the three shock variables uniform 
  # selfadmin$Rat <- ifelse(grepl("^[0][.][1-3]", selfadmin$Rat), sub("([0][.][1-3])", "Shock\\1", selfadmin$Rat), selfadmin$Rat) # must separate bc str_extract ignore case function is deprecated
  # XX SHOCK SHOULD STILL BE 0.1, 0.2, 0.3 (SINCE COHORT 7 AND 8 DON'T HAVE 0.1 AND 0.2)
  
  
  # make experiment name unique (# code from G. Grothendieck (Stack Overflow) )
  uniquify <- function(x) if (length(x) == 1) x else sprintf("%s%02d", x, seq_along(x)) 
  selfadmin$Rat <- ave(selfadmin$Rat, selfadmin$Rat, FUN = uniquify) 
  selfadmin$Rat <- ifelse(grepl("^Shock", selfadmin$Rat, ignore.case = T), gsub("(0)(\\d)", "\\1\\2", selfadmin$Rat), selfadmin$Rat) # remove the decimal pt in between because it causes problems in separation later ## note that later cohorts use shock03 as the only shock, so this may look like shock01
  # after it is uniquified shocks should turn into decimal numbers 
  
  if(grepl("7|8", filename)){
    selfadmin$Rat <- gsub("^Shock", "Shock03", selfadmin$Rat)
  } # if cohort 7 or 8, then turn the shock values to 03 #originally 0.3 but this creates problems in merging  
 
  
  ### transpose data 
  tselfadmin <- data.table::transpose(selfadmin)
  colnames(tselfadmin) <- as.character(tselfadmin[1,])
  tselfadmin$labanimalid = c(rep(NA, which(grepl("^\\d", tselfadmin$RFID)==T)[1]- 1) , grep("(F|M)\\d+", names(selfadmin), value = T), 
                             rep(NA, length(tselfadmin$RFID) - which(grepl("^\\d", tselfadmin$RFID)==T) %>% tail(1)))
  tselfadmin$labanimalid <- str_match(tselfadmin$labanimalid, "(M|F)\\d+")[,1] 
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
  tselfadmin <- tselfadmin[!grepl("^\\d{4}\\-", unlist(tselfadmin[,2])),] # remove the row with dates
  
  # extract [*DATA DICTIONARY*]
  datadictionary <- selfadmin[, 1:3] # vertical formatting is preferred in selfadmin
  colnames(datadictionary) <- c("var_abbv", "var_type", "var_graphtext")
  datadictionary$var_description <- NA
  uniquify_graphtext <- function(x) if (length(x) == 1) x else paste0(sprintf("%s on day %02d", x, seq_along(x)))
  datadictionary$var_graphtext <- ave(datadictionary$var_graphtext, datadictionary$var_graphtext, FUN = uniquify_graphtext)
  
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
  
  if(grepl("8", filename)){
    specificcommentsdates <- specificcommentsdates[2,] %>% as.character  # extract the date column name from every item in list
  } else{
    specificcommentsdates <- sapply(specificcommentsdates, "[", 2) %>% as.character # extract the date column name from every item in list
  }
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
  setcolorder(tspecificcomments, c("labanimalid", "date", "experiment", "conflict", "comment", "flag")) # reorder columns
  
  # extract [*TABLE FOR COMPROMISED ANIMALS*]
  
  # extract [*TABLE FOR DEAD ANIMALS*]
  dead <- tspecificcomments[conflict %like% "(?i)died|dead"] # data.table %ilike% operator doesn't seem ready so using (?i) tag for case insensitive searches
  
  # extract [*TABLE FOR SWITCHED ANIMALS*]
  switches <- tspecificcomments[conflict %like% "(?i)switch(ed)?"]
  
  # add cohort column
  # tselfadmin <- append(tselfadmin, list(data.table(cohort = sub("_cocaine.*", "", filename)))) %>% as.data.table   # get from file name # temporarily remove for the database names, cohort duplicate
  
  # add rat id column
  tselfadmin <- append(tselfadmin, list(data.table(labanimalid = grep("\\D\\d", names(selfadmin), value = T)))) %>% as.data.table
  tselfadmin[,ncol(tselfadmin)] <- NULL
  
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

  
  
  
  # PROCESS IRRITABILITY PORTION
  
  options(scipen = 100, digits = 2)
  irritability <- u01.importxlsx(filename)[[2]] %>%
    as.data.table 
  
  irritability[, names(irritability) := lapply(.SD, as.character)] 
  if(any(grepl("^[[:digit:]]{5,}", irritability$Date))){
    irritability$Date <- as.character(as.POSIXct(as.numeric(irritability$Date) * (60*60*24), origin="1899-12-30", tz="UTC", format="%Y-%m-%d"))
  } # convert Excel character into dates
  
  uniquify <- function(x) if (length(x) == 1) x else sprintf("%s%02d", x, seq_along(x)) 
  irritability$Rat <- ave(irritability$Rat, irritability$Rat, FUN = uniquify) # make unique measure values
  
  irritability <- rbindlist(list(irritability, as.list(names(irritability))), fill=FALSE) # preserve experiment name as rows before tranposing
  
  # print(irritability)
  
  tirritability <- data.table::transpose(irritability)
  colnames(tirritability) <- as.character(tirritability[1,]) # make colnames 
  tirritability <- tirritability[-1,] # remove first row that held colnames
  
  # make colnames clean
  setnames(tirritability, gsub(" ", "_", tolower(names(tirritability))))
  
  # make date and time values (removed a third value (scorer), only found in first cohort and all values were marked DC)
  if(any(grep("Scorer", tirritability$rat))){
    tirritability <- tirritability[!grepl("Scorer", rat),] # remove the row with dates
  }
  
  nm_irr <- grep("^(def|agg)", names(tirritability), value = T) # use these columns to make date columns
  nm1_irr <- paste("date", nm_irr, sep = "_") # make these date columns
  tirritability[ , ( nm1_irr ) := lapply( .SD, function(x) c(grep("^\\d{4}\\-(0?[1-9]|1[012])\\-(0?[1-9]|[12][0-9]|3[01])$", x, value = T)) ) ,  .SDcols = nm_irr ]
  ind_irr <- grep("^(date)", names(tirritability), perl = T) # bring date characters back to posixct
  for (i in seq_along(ind_irr)) {
    set(tirritability, NULL, ind_irr[i], as.POSIXct(tirritability[[ind_irr[i]]], tz = "UTC"))
  }
  tirritability <- tirritability[!grepl("^\\d{4}\\-", unlist(tirritability[,2])),] # remove the row with dates
  
  nm2_irr <- paste("timepoint", nm_irr, sep = "_") # make these time columns
  tirritability[ , ( nm2_irr ) := lapply( .SD, function(x) c(grep("^(before|after)", x, value = T)) ) ,  .SDcols = nm_irr ]
  tirritability <- tirritability[!grepl("^[[:alpha:]]", unlist(tirritability[,2])),] # remove the row with timepoints
  
  # clean up the var types (to numeric) and the character na (to NA) 
  ind <- grep(pattern = "^(def|agg)", names(tirritability))
  for(j in seq_along(ind)){
    set(tirritability, i=which(tirritability[[j]]=="n/a"), j=j, value=NA)
  }
  
  # quant data should be numeric
  for (i in seq_along(ind)) {
    set(tirritability, NULL, ind[i], as.numeric(tirritability[[ind[i]]]))
  }
  
  # extract [*DATA DICTIONARY*]
  irritability_datadictionary <- irritability[, 1:3] # vertical formatting is preferred in selfadmin
  tirritability <- tirritability[-c(1:2), ] # remove data dictionary from data
  
  ### PRINT ALL 
  
  list <- list("tselfadmin" = tselfadmin,
               "datadictionary" = datadictionary,
               "specificcomments" = tspecificcomments,
               "dead" = dead,
               "switches" = switches,
               "tirritability" = tirritability,
               "irritability_datadictionary" = irritability_datadictionary)
  
  return(list)
} 


## Before 08/05/2020
setwd("~/Dropbox (Palmer Lab)/Olivier_George_U01/DATA Updated")
cohortfiles <- list.files(pattern = "*.xlsx")
# cohort 6 is missing 
for(i in 1:length(cohortfiles)){
  # if(i == 4) next 
  if(i > 5){
    assign(paste0("cohort", i + 1), olivierfiles(cohortfiles[i]))
  } else
    assign(paste0("cohort", i), olivierfiles(cohortfiles[i]))
}

## After 08/05/2020
setwd("~/Dropbox (Palmer Lab)/Olivier_George_U01/GWAS Self Administration Data/Cocaine Data")
cohortfiles_sa_2 <- list.files(pattern = "*.xls[xm]")
# cohort 6 is missing 
for(i in 1:length(cohortfiles_sa_2)){
  # if(i == 4) next 
  if(i > 5){ # bc cohort 6 is skipped
    assign(paste0("cohort_sa", i + 1), olivierfiles(cohortfiles[i]))
  } else
    assign(paste0("cohort_sa", i), olivierfiles(cohortfiles[i]))
}

## XXXXXXXXXX go back to fix cohort8, the [[7]] bc of the formatting, it resuts in a find for  grep("^\\D", selfadmin_df$rfid, value = T)

# function breaks at 
cohort2 <- olivierfiles("C02_cocaine.xlsx")
cohort3 <- olivierfiles("C03_cocaine.xlsx")
cohort4 <- olivierfiles("C04_cocaine.xlsx") # specific comment in non specific column 
cohort5 <- olivierfiles("C05_cocaine.xlsx") # note about renumbering -- can we ignore
cohort7 <- olivierfiles("C07_cocaine.xlsx")
cohort8 <- olivierfiles("C08_cocaine.xlsx") # comments aren't working because there is a row that needs to be removed

### spot checking 
# cohort 5 typo in LGA19 1981 year should be 2018 AND
# cohort 5 conversion in date doesn't take the SHA02 excel string so this code leaves it as NA 
cohort5$tselfadmin[, date_lga19 := lubridate::ymd("2018-09-19")]
cohort5$tselfadmin[, date_sha02 := lubridate::ymd("2018-07-31")]

## don't use this line bc this converts the data.table to dataframe and throws off the rbindlist in create_dataabasestructure to make allcohorts2 # cohort5$tselfadmin %<>% mutate(date_lga19 = lubridate::ymd("2018-09-19"), date_sha02 = lubridate::ymd("2018-07-31"))


### EXTRACT THE MAPPING FILES FROM THEIR DROPBOX
setwd("~/Dropbox (Palmer Lab)/Olivier_George_U01/Rat Information/Cocaine")
rat_info_xl_filenames <- list.files(pattern = "*.xlsx")

rat_info_allcohort_xl_df <- lapply(rat_info_xl_filenames, function(x){
  path_sheetnames <- excel_sheets(x)
  df <- lapply(excel_sheets(path = x), read_excel, path = x) # including this, just in case the info file ever moves out of order
  names(df) <- path_sheetnames
  info_name <- grep("info", path_sheetnames, ignore.case = T, value = T) # allows for small changes, like info sheet vs info sheets
  df_info <- df[[info_name]]
  
  df_info <- df_info %>% mutate_all(as.character) # prevent any mismatched class attributes
  df_info <- df_info %>%
    mutate(naive = NA, 
           naive = replace(naive, is.na(RAT), "Empty"),
           naive = replace(naive, grepl("Naive", RAT, ignore.case=F), "Naive")) %>% 
    tidyr::fill(naive) %>% 
    mutate(naive = replace(naive, !grepl("Naive", naive), NA)) %>% 
    subset(grepl("^\\d", RFID)) # tackle the naive cases

  return(df_info)
})
names(rat_info_allcohort_xl_df) <- rat_info_xl_filenames
rat_info_allcohort_xl_df %<>% rbindlist(fill = T, idcol = "cohort") %<>% 
  mutate(cohort = str_extract(cohort, "C\\d{2}")) %<>%
  clean_names() %<>%
  mutate(labanimalid = str_extract(rat, "[MF]\\d+.*"))


# which animals are in the excel sheets but not in the wfu master tables 
anti_join(rat_info_allcohort_xl_df[, c("cohort","rfid", "rat", "naive")], 
          WFU_OlivierCocaine_test_df[,c("cohort", "rfid", "sires")],  by = "rfid") %>% View()
# which animals are in the wfu master tables but not in excel sheets 
anti_join(rat_info_allcohort_xl_df[, c("cohort","rfid", "rat", "naive")], 
          WFU_OlivierCocaine_test_df[,c("cohort", "rfid", "sires")],  by = "rfid") %>% View()









### EXTRACT THE COMPUTER NOTES FROM THEIR DROPBOX
setwd("~/Dropbox (Palmer Lab)/GWAS (1)")
computernotes_coc <- u01.importxlsx("computer notes.xlsx")[[1]] %>% 
  gather(exp, computernote, SHA01:cohort_notes) %>% 
  clean_names() %>% 
  dplyr::filter(grepl("^C", cohort)) %>% 
  mutate(computernote = replace(computernote, computernote == "NA", "Did not run as part of protocol"))
  # naniar::replace_with_na(replace = list(computernote = "NA"))





## EXTRACT THE RAT WEIGHTS 
# remove the weights associated before the shipping dates 
rats_allcohorts_weights <- rat_info_allcohort_xl_df %>% select(cohort, rat, rfid, d_o_b, matches("arrival_date$|surgery_date$|weight")) %>% 
                                                                 split(., .$cohort) %>% 
                                                                 lapply(function(x){
                                                                   df <- x %>% select_if(~sum(!is.na(.)) > 0) %>% 
                                                                     mutate(surgery_date = openxlsx::convertToDate(surgery_date)) %>% 
                                                                            # ,
                                                                            # arrival_date = as.Date(arrival_date)) %>% 
                                                                     select(-matches("pre_shipment")) ## remove if varname explicitly includes the preshipment label
                                                                   names(df) <- gsub("sugery", "surgery", names(df)) # fix the surgery name typo
                                                                   names(df) <- gsub("d_o_b", "dob", names(df)) # fix the surgery name typo
                                                                   
                                                                   # separate the dates from the variable names 
                                                                   weights_dates <- gsub("weight_\\d+_", "", names(df)) %>% grep("\\d+_\\d+_\\d+", ., value = T) %>% t() %>% data.frame()
                                                                   names(weights_dates) <- paste0("dateweight_", 1:length(weights_dates))
                                                                   
                                                                   
                                                                   # make the variable names contain only numbers
                                                                   names(df) <- gsub("(weight(_\\d+_surgery)?)_\\d+_\\d+_.*", "\\1", names(df)) # leave only weight and weight_surgery as varnames 
                                                                   names(df) <- gsub("weight_\\d+_surgery", "surgery_weight", names(df))
                                                                   # names(df) <- gsub("weight_\\d+_surgery", "weight_surgery", names(df))
                                                                   df <- df %>% clean_names # append numbers after weight
                                                                   names(df) <- gsub("^weight$", "weight_1", names(df))
                                                                   
                                                                   bound_df <- cbind(df, weights_dates)
                                                                   bound_df <- bound_df %>% 
                                                                     mutate_at(vars(matches("date")), ~gsub("_", "-",.)) %>%
                                                                     mutate_at(vars(starts_with("date")), ~as.Date(., "%m-%d-%Y") %>% as.character) %>% 
                                                                     mutate_all(as.character) %>% 
                                                                     naniar::replace_with_na_all(~.x %in% c("n/a", "na", "NA"))  
                                                                   # %>% 
                                                                     # mutate_at(vars(matches("date|dob")), as.Date) %>% 
                                                                     # mutate(age_surgery = difftime(surgery_date, dob, units = c("days")))
                                                                   # %>% 
                                                                     # mutate_at(vars(matches("date")), list(age = ~difftime(., dob, units = c("days"))))
                                                                                                           # %>% as.numeric %>% as.character))

                                                                   return(bound_df)
                                                                   
                                                                 }) %>% rbindlist(fill = T, use.names = T)%>% 
  mutate_at(vars(matches("date|dob")), as.Date) %>% 
  
# %>% 
  # mutate_at(vars(matches("date")), list(age = ~difftime(., dob, units = c("days")) %>% as.numeric)) %>% 
  mutate_all(as.character)
# %>% 
#   mutate(weight_3 = replace(weight_3, rfid == "933000320045777", "180"),
#          ) %>% ## confirmed with Brent 08/05/2020

names(rats_allcohorts_weights) <- gsub("date_(.*_age$)", "\\1", names(rats_allcohorts_weights)) 
names(rats_allcohorts_weights) <- gsub("surgery_date_age", "surgery_age", names(rats_allcohorts_weights)) 

# quick qc before uploading
# checking min and max 
rats_allcohorts_weights %>% mutate_at(vars(starts_with("weight")), as.numeric) %>% 
  summary
rats_allcohorts_weights %>% mutate_at(vars(ends_with("age")), as.numeric) %>% 
  summary

rats_allcohorts_weights %>% mutate_at(vars(starts_with("weight")), as.numeric) %>% 
  subset(weight_3 < 30)
rats_allcohorts_weights %>% mutate_at(vars(starts_with("weight")), as.numeric) %>% 
  subset(weight_1_age < 20)



rats_allcohorts_weights %>% subset(is.na(surgery_date)&!is.na(weight_surgery)) %>% str




# trying to get rid of the preshipment weights 
rats_allcohorts_weights_long <- rats_allcohorts_weights %>% 
  gather(key = "datename", value = "date", starts_with("date")) %>% 
  subset(!is.na(date)) %>% 
  gather("weightname", "weight", starts_with("weight")) %>% 
  mutate(datename = gsub('[^0-9]', "", datename), weightname = gsub('[^0-9]', "", weightname)) %>% 
  dplyr::filter(datename == weightname) %>%  select(-ends_with("name")) %>% arrange(rat)

## fix dateweight_1 and no arrival date info


setwd("~/Desktop/Database/csv files/u01_olivier_george_cocaine")
# run this code once confirmed preshipment weights
# write.csv(rats_allcohorts_weights, file = "rats_cohorts01_11_weights.csv")






## extract the intermediate values
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

cocaine_intermediates_xl_lga_df <- cocaine_intermediates_xl_lga %>% 
  lapply(function(x){
    names(x) <- c("labanimalid", paste0("lga_", str_pad(1:14, "2", "left", "0"))) 
    return(x)
  }) %>% 
  rbindlist(idcol = "cohort", fill = T) %>%  
  mutate(labanimalid = str_extract(labanimalid, "[MF]\\d+")) %>% 
  mutate_at(vars(-matches("cohort|labanimal")), as.numeric) %>% 
  subset(!is.na(labanimalid)) %>% 
  mutate_at(vars(matches("lga_1[1-4]$")), list(esc = ~.-lga_01)) %>% 
  mutate(esc11_14_mean = rowMeans(select(., ends_with("_esc")), na.rm = TRUE)) %>% 
  select(cohort, labanimalid, matches("_esc$"), matches("_mean$"))


## pr 

cocaine_intermediates_xl_pr <- u01.importxlsx("Addiction indices for C01-C07 Cocaine.xlsx") %>% 
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
cocaine_intermediates_xl_pr$C01<- cocaine_intermediates_xl_pr$C01 %>% select(matches("rats$"), matches("^day_\\d_3$"))
cocaine_intermediates_xl_pr$C02<- cocaine_intermediates_xl_pr$C02 %>% select(matches("rats$"), matches("^day_\\d_3$"))
cocaine_intermediates_xl_pr$C03<- cocaine_intermediates_xl_pr$C03 %>% select(matches("rats_1$"), matches("^pr")) 
cocaine_intermediates_xl_pr$C04<- cocaine_intermediates_xl_pr$C04 %>% select(matches("rat$"), matches("^pr")) 
cocaine_intermediates_xl_pr$C05<- cocaine_intermediates_xl_pr$C05 %>% select(matches("rats_1$"), matches("^pr"))
cocaine_intermediates_xl_pr$C07<- cocaine_intermediates_xl_pr$C07 # no change necessary

cocaine_intermediates_xl_pr_df <- cocaine_intermediates_xl_pr %>% 
  lapply(function(x){
    names(x) <- c("labanimalid", paste0("pr_", str_pad(1:2, "2", "left", "0"))) 
    return(x)
  }) %>% 
  rbindlist(idcol = "cohort", fill = T) %>%  
  mutate(labanimalid = str_extract(labanimalid, "[MF]\\d+")) %>% 
  mutate_at(vars(-matches("cohort|labanimal")), as.numeric) %>% 
  subset(!is.na(labanimalid)) %>% 
  mutate(pr_max = pmax(pr_01, pr_02)) %>% 
  mutate(pr_mean = rowMeans(select(., starts_with("pr")), na.rm = TRUE)) %>% 
  select(labanimalid, pr_max, pr_mean)

## shock 
cocaine_intermediates_xl_shock <- u01.importxlsx("Addiction indices for C01-C07 Cocaine.xlsx") %>% 
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
      select(starts_with("rat"), matches("x0_"), starts_with("shock"), matches("BSL"))
    return(x)
  }) 
cocaine_intermediates_xl_shock$C01<- cocaine_intermediates_xl_shock$C01 %>% select(matches("rats$"), matches("x0_"), -matches("0_1"))
cocaine_intermediates_xl_shock$C02<- cocaine_intermediates_xl_shock$C02 %>% select(matches("rats$"), matches("x0_"), -matches("0_1"))
cocaine_intermediates_xl_shock$C03<- cocaine_intermediates_xl_shock$C03 %>% select(matches("rats_1$"), starts_with("shock"), -matches("0_1"))
cocaine_intermediates_xl_shock$C04<- cocaine_intermediates_xl_shock$C04 %>% select(matches("rat$"), starts_with("shock"), -matches("0_1"))
cocaine_intermediates_xl_shock$C05<- cocaine_intermediates_xl_shock$C05 %>% select(matches("rats_1$"), starts_with("shock"), -matches("0_1"))
cocaine_intermediates_xl_shock$C07<- cocaine_intermediates_xl_shock$C07 %>% select(matches("rats$"), matches("bsl"), matches("^shock$"), -matches("0_1"))

cocaine_intermediates_xl_shock_df <- cocaine_intermediates_xl_shock %>% 
  lapply(function(x){
    names(x) <- c("labanimalid", paste0("shock_", str_pad(2:3, "2", "left", "0"))) 
    return(x)
  }) %>% 
  rbindlist(idcol = "cohort", fill = T) %>%  
  mutate(labanimalid = str_extract(labanimalid, "[MF]\\d+")) %>% 
  mutate_at(vars(-matches("cohort|labanimal")), as.numeric) %>% 
  subset(!is.na(labanimalid)) %>% 
  select(labanimalid, shock_03)


cocaine_intermediates_xl_bind <- left_join(cocaine_intermediates_xl_lga_df, cocaine_intermediates_xl_pr_df, by = "labanimalid") %>% 
  left_join(cocaine_intermediates_xl_shock_df, by = "labanimalid")

rm(list=ls(pattern="irr")) # conditionally clean the environment
rm(list = ls())
