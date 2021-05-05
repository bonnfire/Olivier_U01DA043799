## replicating cocaine excel extraction

# run this line after running all sections after cohort1
# create object for all values
oliviercocaine_excel_all <- list(
  "C01"=selfadmin_xl_cohort1,
  "C03"=selfadmin_xl_cohort3,
  "C04"=selfadmin_xl_cohort4,
  "C05"=selfadmin_xl_cohort5,
  "C06"=selfadmin_xl_cohort6,
  "C07"=selfadmin_xl_cohort7,
  "C09"=selfadmin_xl_cohort9,
  "C10"=selfadmin_xl_cohort10
) %>% rbindlist(idcol = "cohort", fill = T) %>%
  select_if(~sum(!is.na(.)) > 0) %>% 
  clean_names() %>% 
  rename("labanimalid" = "rat") %>% 
  mutate(measurement = gsub("^ACTIVE.*", "active", measurement),
         measurement = gsub("^INACTIVE.*", "inactive", measurement),
         measurement = gsub("^PR.*", "pr_breakpoint", measurement),
         measurement = gsub("^REWARDS.*", "rewards", measurement)) 

oliviercocaine_excel_all <- oliviercocaine_excel_all %>% 
  subset(!grepl("\\D+:", labanimalid)) %>% 
  naniar::replace_with_na_all(condition = ~.x %in% c("N/A")) %>% 
  subset(!(is.na(labanimalid)&cohort %in% c("C01", "C03", "C04"))) # add to investigate other cohorts before removing, ensure with naniar::vis_miss()


rm(list = ls(pattern = "selfadmin_(xl|rewards)|cohortinfo|comments_df|selfadmin(_df)?|nm(1)?|dates"))


## to extract the mapping excels 
cohortinfofiles <- list.files(path = "~/Dropbox (Palmer Lab)/Olivier_George_U01/Rat Information/Cocaine", pattern = "*.xlsx", full.names = T)


u01.importxlsx <- function(xlname){
  df <- lapply(excel_sheets(path = xlname), read_excel, path = xlname)
  names(df) <- excel_sheets(xlname)
  return(df)
}

cocaine_metadata <- lapply(cohortinfofiles, openxlsx::read.xlsx)
names(cocaine_metadata) <- cohortinfofiles
cocaine_metadata_df <- cocaine_metadata %>% 
  rbindlist(fill = T, idcol = "cohort") %>% 
  subset(grepl("^\\d{15}$", RFID)) %>% 
  mutate(cohort = gsub(".*(C\\d+).*", "\\1", cohort),
         sex = gsub(".*([MF]).*", "\\1", RAT)) %>% 
  clean_names 
# %>% 
# select(cohort, rat, rfid, sex)

cohortfiles_xl_c01_11 <- list.files(path = "~/Dropbox (Palmer Lab)/Olivier_George_U01/GWAS Self Administration Data/Cocaine Data", pattern = "*.xls(x|m)", full.names = T) %>% grep("C(0[1-9]|1[0|1])", ., value = T)


########################
# COHORT 1 
########################
setwd("~/Dropbox (Palmer Lab)/Olivier_George_U01/GWAS Self Administration Data/Cocaine Data")
filename <- cohortfiles_xl_c01_11[1]
selfadmin <- u01.importxlsx(filename)[[1]] %>%
  as.data.table
selfadmin[ selfadmin == "n/a" ] <- NA # change all character n/a to actual NA
# set correct column names 
# create date columns
dates <- grep("^\\d+", names(selfadmin), value = T) # use these columns to make date columns # ignore the ...\\d columns
dates <- as.character(as.POSIXct(as.numeric(dates) * (60*60*24), origin="1899-12-30", tz="UTC", format="%Y-%m-%d")) # convert Excel character into dates

setnames(selfadmin, toupper(as.character(selfadmin[4,]) )) # now that dates are moved into separate vector, remove from the column names 
# setnames(selfadmin,  sub("(PR\\d+)(TRTMEN)?(T)[0]?([1-9])", paste0("\\1_\\3", str_pad("\\4", 3, "left", "0")), names(selfadmin)) )
selfadmin <- selfadmin[-4,]
selfadmin <- remove_empty(selfadmin, "cols") # janitor::remove_empty_cols() deprecated
# selfadmin <- selfadmin[1:grep("average", selfadmin$RAT)[1],] # subset only the rewards table


nm <- names(selfadmin)[-c(1:2)]  
nm <- ifelse(grepl("^\\D{2}A", nm), as.character(stringr::str_match(nm,"^\\D{2}A")), nm) # reg exp fixes issuse of extracting mA
nm <- ifelse(grepl("^PR", nm), as.character(stringr::str_match(nm,"PR")), nm)
nm <- ifelse(grepl("PreShock[[:space:]][(LgA15)]", nm), "LgA", nm) # applies to cohort 8; technician error 
nm <- gsub(" |[(]|[)]|[-]", "", nm) # remove all unwanted characters
nm <- gsub("^(1hr|PreShock).+", "preshock", nm, ignore.case = T) # XX GET HIS INPUT ON THIS -- currently one HOUR admin; is preshock (in cohort7 also the same)
nm <- ifelse(grepl("^Shock", nm, ignore.case = T), "Shock", nm) # make date columns for this vector of exp names 
nm <- nm[-35]


nm1 <- paste("date", nm, sep = "_") # make these date columns
selfadmin[ , ( nm1 ) := lapply( dates, function(x) rep(x, each = .N) ) ] # make the date columns 

comments_df <- selfadmin[which(selfadmin$RAT %in% c("COMMENT", "CONFLICT", "RESOLUTION"))] #extact comments
comments_df <- comments_df %>% select(-matches("RFID|date")) %>% t() %>% 
  as.data.frame() %>% 
  rownames_to_column()
setnames(comments_df, append("EXP", comments_df[1, 2:4] %>% t() %>% unlist() %>% as.character) %>% tolower)
comments_df <- comments_df[-1,]
# selfadmin <- selfadmin[!which(selfadmin$RAT %in% c("COMMENT", "CONFLICT", "RESOLUTION"))]

selfadmin_exps <- grep("REWARDS|ACTIVE|INACTIVE|PR", selfadmin$RAT)
selfadmin_split <- split(selfadmin, cumsum(1:nrow(selfadmin) %in% selfadmin_exps))
names(selfadmin_split) <- lapply(selfadmin_split, function(x){ x$RAT %>% head(1)}) %>% unlist() %>% as.character()
selfadmin_split <- lapply(selfadmin_split, function(x){ x %>% dplyr::filter(grepl("^\\d", RFID))})
selfadmin_df <- selfadmin_split %>% rbindlist(idcol = "measurement") %>% dplyr::filter(measurement != "COMMENT")
selfadmin_rewards_cohort1 <- selfadmin_df %>% dplyr::filter(measurement == "ACTIVE")
selfadmin_xl_cohort1 <- selfadmin_df %>% ## Lani notes on 03/12 "Rewards data for this cohort is under Active; Inactive data is inactive"
  mutate(measurement = replace(measurement, measurement == "ACTIVE", "REWARDS"))

# fixing shifted rfid's until Lani updates the Excel files
# selfadmin_xl_cohort1 <- selfadmin_xl_cohort1 %>% 
#   left_join(cocaine_metadata_df %>% select(-cohort, -sex), by = c("RAT" = "rat")) %>% 
#   rowwise() %>%
#   mutate(RFID = replace(RFID, RFID != rfid, rfid)) %>% 
#   select(-rfid) %>% 
#   ungroup()

selfadmin_xl_cohort1 %>% distinct(RAT, RFID) %>% get_dupes(RAT)


########################
# COHORT 3
########################  

# this file doesn't have rfid's so will need to extract from the replacement table
filename <- cohortfiles_xl_c01_11[2]
selfadmin <- u01.importxlsx(filename)[[1]] %>%
  as.data.table
selfadmin[ selfadmin == "n/a" ] <- NA # change all character n/a to actual NA

# EXTRACT COMMENTS bc comments are stuck as header
comments_df <- selfadmin[1:2,] %>% t() %>% as.data.frame() 
setDT(comments_df, keep.rownames = TRUE)[]
setnames(comments_df, comments_df[1,] %>% t() %>% unlist() %>% as.character %>% tolower)
comments_df <- comments_df[-1, ]
comments_df$comments <- comments_df$comments %>% gsub("[.]{3}\\d", NA,. ) 

# set correct column names 
# create date columns
dates <- grep("^\\d+", selfadmin[3,] %>% unlist() %>% as.character(), value = T) # use these columns to make date columns # ignore anything else
dates <- as.character(as.POSIXct(as.numeric(dates) * (60*60*24), origin="1899-12-30", tz="UTC", format="%Y-%m-%d")) # convert Excel character into dates
setnames(selfadmin, toupper(as.character(selfadmin[4,]) )) # now that dates are moved into separate vector, remove from the column names 
setnames(selfadmin,  mgsub::mgsub(names(selfadmin), c("LGA15.*", "^PRT.+1$", "^PRT.+2$", "^PRT.+3$", "^PRT.+4$"), c("LGA15", "PR03_T01", "PR04_T02", "PR05_T03", "PR06_T04")))
names(selfadmin)[1] <- "RAT"
selfadmin <- remove_empty(selfadmin, "cols") # janitor::remove_empty_cols() deprecated
# selfadmin <- selfadmin[1:grep("average", selfadmin$RAT)[1],] # subset only the rewards table
# selfadmin <- selfadmin[-1,]

nm <- names(selfadmin)[-1] # make date columns for this vector of exp names  ## MISSING THE RFID COLUMN SO THAT IS WHY [-1] INSTEAD OF [-c(1:2)]
nm1 <- paste("date", nm, sep = "_") # make these date columns
selfadmin[ , ( nm1 ) := lapply( dates, function(x) rep(x, each = .N) ) ] # make the date columns 

selfadmin_exps <- grep("REWARDS|ACTIVE|INACTIVE|PR$", selfadmin$RAT)
selfadmin_split <- split(selfadmin, cumsum(1:nrow(selfadmin) %in% selfadmin_exps))
names(selfadmin_split) <- lapply(selfadmin_split, function(x){ x$RAT %>% head(1)}) %>% unlist() %>% as.character()
selfadmin_split <- lapply(selfadmin_split, function(x){ x %>% dplyr::filter(grepl("^[MF]\\d+", RAT))})
selfadmin_df <- selfadmin_split %>% rbindlist(idcol = "measurement")
selfadmin_rewards_cohort3 <- selfadmin_df %>% dplyr::filter(measurement == "REWARDS")
selfadmin_xl_cohort3 <- selfadmin_df

# EXTRACT RFID
filename <- cohortinfofiles[2]
cohortinfo <- u01.importxlsx(filename)[[1]] %>% as.data.frame()
names(cohortinfo) <-  mgsub::mgsub(names(cohortinfo),
                                   c(" |\\.", "#", "Transponder ID", "Date of Wean|Wean Date","Animal", "Shipping|Ship", "Dams"),
                                   c("", "Number", "RFID", "DOW","LabAnimal", "Shipment", "Dames")) 
# join rfid onto the rewards
selfadmin_rewards_cohort3 <- selfadmin_rewards_cohort3 %>% 
  left_join(cohortinfo[, c("RAT", "RFID")], ., by = c("RAT"))

selfadmin_xl_cohort3 <- selfadmin_xl_cohort3 %>% 
  left_join(cohortinfo[, c("RAT", "RFID")], ., by = c("RAT")) 


########################
# COHORT 4
########################
# this file doesn't have rfid's so will need to extract from the replacement table

filename <- cohortfiles_xl_c01_11[3]
selfadmin <- u01.importxlsx(filename)[[1]] %>%
  as.data.table
selfadmin[ selfadmin == "n/a" ] <- NA # change all character n/a to actual NA

# set correct column names 
# create date columns
dates <- grep("^\\d+", names(selfadmin), value = T) # use these columns to make date columns # ignore the ...\\d columns
dates <- as.character(as.POSIXct(as.numeric(dates) * (60*60*24), origin="1899-12-30", tz="UTC", format="%Y-%m-%d")) # convert Excel character into dates

setnames(selfadmin, toupper(as.character(selfadmin[1,]) )) # now that dates are moved into separate vector, remove from the column names 
setnames(selfadmin,  mgsub::mgsub(names(selfadmin), c("PR2", "^T.+1$", "^T.+2$", "^T.+3$", "^T.+4$", "LGA([1-9]{1})$"), c("PR02", "PR03_T01", "PR04_T02", "PR05_T03", "PR06_T04", "LGA0\\1")))
names(selfadmin)[1] <- "RAT"
selfadmin <- selfadmin[-1,]
selfadmin <- remove_empty(selfadmin, "cols") # janitor::remove_empty_cols() deprecated
# selfadmin <- selfadmin[1:grep("average", selfadmin$RAT)[1],] # subset only the rewards table

nm <- names(selfadmin)[-1] # make date columns for this vector of exp names  ## MISSING THE RFID COLUMN SO THAT IS WHY [-1] INSTEAD OF [-c(1:2)] 
nm1 <- paste("date", nm, sep = "_") # make these date columns
selfadmin[ , ( nm1 ) := lapply( dates, function(x) rep(x, each = .N) ) ] # make the date columns 

#extract comments
comments_df <- selfadmin[which(selfadmin$RAT %in% c("COMMENT", "CONFLICT", "RESOLUTION"))]
comments_df <- comments_df %>% select(-matches("RFID|date")) %>% t() %>% 
  as.data.frame() %>% 
  rownames_to_column()
setnames(comments_df, append("EXP", comments_df[1, 2:4] %>% t() %>% unlist() %>% as.character) %>% tolower)
comments_df <- comments_df[-1,]
# selfadmin <- selfadmin[!which(selfadmin$RAT %in% c("COMMENT", "CONFLICT", "RESOLUTION"))]

selfadmin_exps <- grep("REWARDS|ACTIVE|INACTIVE|PR$", selfadmin$RAT)
selfadmin_split <- split(selfadmin, cumsum(1:nrow(selfadmin) %in% selfadmin_exps))
names(selfadmin_split) <- lapply(selfadmin_split, function(x){ x$RAT %>% head(1)}) %>% unlist() %>% as.character()
selfadmin_split <- lapply(selfadmin_split, function(x){ x %>% dplyr::filter(grepl("^[MF]\\d+", RAT))})
selfadmin_df <- selfadmin_split %>% rbindlist(idcol = "measurement") %>% dplyr::filter(measurement != "COMMENT")
selfadmin_rewards_cohort4 <- selfadmin_df %>% dplyr::filter(measurement == "REWARDS")
selfadmin_xl_cohort4 <- selfadmin_df



# EXTRACT RFID
filename <- cohortinfofiles[3]
cohortinfo <- u01.importxlsx(filename)[[1]] %>% as.data.frame()
names(cohortinfo) <-  mgsub::mgsub(names(cohortinfo),
                                   c(" |\\.", "#", "Transponder ID", "Date of Wean|Wean Date","Animal", "Shipping|Ship", "Dams"),
                                   c("", "Number", "RFID", "DOW","LabAnimal", "Shipment", "Dames")) 
# join rfid onto the rewards
selfadmin_rewards_cohort4 <- selfadmin_rewards_cohort4 %>% 
  left_join(cohortinfo[, c("RAT", "RFID")], ., by = c("RAT"))

selfadmin_xl_cohort4 <- selfadmin_xl_cohort4 %>% 
  left_join(cohortinfo[, c("RAT", "RFID")], ., by = c("RAT"))

########################
# COHORT 5
########################
# filename <- cohortfiles_xl_c01_11[4]

filename <- "~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/U01/Olivier_George_U01DA044451 (cocaine)/excel_and_csv_filesCocaine GWAS C05 Data.xlsx" # XX temp for Lani 
selfadmin <- u01.importxlsx(filename)[[1]] %>%
  # u01.importxlsx(filename) %>%
  as.data.table
selfadmin[ selfadmin == "n/a" ] <- NA # change all character n/a to actual NA

# set correct column names 
# create date columns
dates <- grep("^\\d+", names(selfadmin), value = T) # use these columns to make date columns # ignore the ...\\d columns
dates <- as.character(as.POSIXct(as.numeric(dates) * (60*60*24), origin="1899-12-30", tz="UTC", format="%Y-%m-%d")) # convert Excel character into dates

setnames(selfadmin, toupper(as.character(janitor::make_clean_names(selfadmin[1,]))) %>% gsub("(\\D{2})_(\\D)_", "\\1\\2_", .) ) # now that dates are moved into separate vector, remove from the column names 
setnames(selfadmin,  mgsub::mgsub(names(selfadmin), c("^PR$", "PR_2", "^TREATMENT_PR$", "^T.+2$", "^T.+3$", "^T.+4$", "SHA.*?([0-9]{2})$", "LGA.*?([0-9]{2})$"), c("PR01", "PR02", "PR03_T01", "PR04_T02", "PR05_T03", "PR06_T04", "SHA\\1", "LGA\\1")))
names(selfadmin)[1:2] <- c("RAT", "RFID")
selfadmin <- selfadmin[-1,]
selfadmin <- remove_empty(selfadmin, "cols") # janitor::remove_empty_cols() deprecated
# selfadmin <- selfadmin[1:grep("average", selfadmin$RAT)[1],] # subset only the rewards table

nm <- names(selfadmin)[-c(1:2)] # make date columns for this vector of exp names  ## MISSING THE RFID COLUMN SO THAT IS WHY [-1] INSTEAD OF [-c(1:2)] 
nm1 <- paste("date", nm, sep = "_") # make these date columns
selfadmin[ , ( nm1 ) := lapply( dates, function(x) rep(x, each = .N) ) ] # make the date columns 

#extract comments
comments_df <- selfadmin[which(selfadmin$RAT %in% c("COMMENT", "CONFLICT", "RESOLUTION"))]
comments_df <- comments_df %>% select(-matches("RFID|date")) %>% t() %>% 
  as.data.frame() %>% 
  rownames_to_column()
setnames(comments_df, append("EXP", comments_df[1, 2:3] %>% t() %>% unlist() %>% as.character) %>% tolower)
comments_df <- comments_df[-1,]
# selfadmin <- selfadmin[!which(selfadmin$RAT %in% c("COMMENT", "CONFLICT", "RESOLUTION"))]

selfadmin_exps <- grep("REWARDS|ACTIVE|INACTIVE|PR$", selfadmin$RAT)
selfadmin_split <- split(selfadmin, cumsum(1:nrow(selfadmin) %in% selfadmin_exps))
names(selfadmin_split) <- lapply(selfadmin_split, function(x){ x$RAT %>% head(1)}) %>% unlist() %>% as.character()
selfadmin_split <- lapply(selfadmin_split, function(x){ x %>% dplyr::filter(grepl("^[MF]\\d+", RAT))})
selfadmin_df <- selfadmin_split %>% rbindlist(idcol = "measurement") %>% dplyr::filter(measurement != "COMMENT")
selfadmin_rewards_cohort5 <- selfadmin_df %>% dplyr::filter(measurement == "REWARDS")
selfadmin_xl_cohort5 <- selfadmin_df


########################
# COHORT 6
########################

filename <- cohortfiles_xl_c01_11[5]
selfadmin <- u01.importxlsx(filename)[[1]] %>%
  as.data.table
selfadmin[ selfadmin == "n/a" ] <- NA # change all character n/a to actual NA

# set correct column names 
# create date columns
dates <- grep("^\\d+", names(selfadmin), value = T) # use these columns to make date columns # ignore the ...\\d columns
dates <- as.character(as.POSIXct(as.numeric(dates) * (60*60*24), origin="1899-12-30", tz="UTC", format="%Y-%m-%d")) # convert Excel character into dates

setnames(selfadmin, toupper(as.character(selfadmin[1,]) )) # now that dates are moved into separate vector, remove from the column names 
setnames(selfadmin,  mgsub::mgsub(names(selfadmin), c("PR1", "PR2", "^T.+1$", "^T.+2$", "^T.+3$", "^T.+4$", "SHA.*?([0-9]{2})$", "LGA.*?([0-9]{2})$"), c("PR01", "PR02", "PR03_T01", "PR04_T02", "PR05_T03", "PR06_T04", "SHA\\1", "LGA\\1")))
names(selfadmin)[1:2] <- c("RAT", "RFID")
selfadmin <- selfadmin[-1,]
selfadmin <- remove_empty(selfadmin, "cols") # janitor::remove_empty_cols() deprecated
# selfadmin <- selfadmin[1:grep("average", selfadmin$RAT)[1],] # subset only the rewards table

nm <- names(selfadmin)[-c(1:2)] # make date columns for this vector of exp names  ## MISSING THE RFID COLUMN SO THAT IS WHY [-1] INSTEAD OF [-c(1:2)] 
nm1 <- paste("date", nm, sep = "_") # make these date columns
selfadmin[ , ( nm1 ) := lapply( dates, function(x) rep(x, each = .N) ) ] # make the date columns 

#extract comments
comments_df <- selfadmin[which(selfadmin$RAT %in% c("COMMENT", "CONFLICT", "RESOLUTION"))]
comments_df <- comments_df %>% select(-matches("RFID|date")) %>% t() %>% 
  as.data.frame() %>% 
  rownames_to_column()
setnames(comments_df, append("EXP", comments_df[1, 2:3] %>% t() %>% unlist() %>% as.character) %>% tolower)
comments_df <- comments_df[-1,]
# selfadmin <- selfadmin[!which(selfadmin$RAT %in% c("COMMENT", "CONFLICT", "RESOLUTION"))]

selfadmin_exps <- grep("REWARDS|ACTIVE|INACTIVE|PR$", selfadmin$RAT)
selfadmin_split <- split(selfadmin, cumsum(1:nrow(selfadmin) %in% selfadmin_exps))
names(selfadmin_split) <- lapply(selfadmin_split, function(x){ x$RAT %>% head(1)}) %>% unlist() %>% as.character()
selfadmin_split <- lapply(selfadmin_split, function(x){ x %>% dplyr::filter(grepl("^[MF]\\d+", RAT))})
selfadmin_df <- selfadmin_split %>% rbindlist(idcol = "measurement") %>% dplyr::filter(measurement != "COMMENT")
selfadmin_rewards_cohort6 <- selfadmin_df %>% dplyr::filter(measurement == "REWARDS")

# move columns over because of formatting in excel
selfadmin_df <- selfadmin_df %>% 
  rowwise() %>% 
  mutate(PR01 = replace(PR01, measurement == "PR", SHA01), 
         PR02 = replace(PR02, measurement == "PR", SHA02), 
         PR03_T01 = replace(PR03_T01, measurement == "PR", SHA03),
         PR04_T02 = replace(PR04_T02, measurement == "PR", SHA04),
         PR05_T03 = replace(PR05_T03, measurement == "PR", PR01),
         PR06_T04 = replace(PR06_T04, measurement == "PR", LGA01)) %>% 
  mutate(SHA01 = replace(SHA01, measurement == "PR", NA), 
         SHA02 = replace(SHA02, measurement == "PR", NA), 
         SHA03 = replace(SHA03, measurement == "PR", NA),
         SHA04 = replace(SHA04, measurement == "PR", NA),
         LGA01 = replace(LGA01, measurement == "PR", NA)) %>% 
  ungroup()
selfadmin_xl_cohort6 <- selfadmin_df


########################
# COHORT 7
########################

filename <- cohortfiles_xl_c01_11[6]
selfadmin <- u01.importxlsx(filename)[[1]] %>%
  as.data.table
selfadmin[ selfadmin == "n/a" ] <- NA # change all character n/a to actual NA

# set correct column names 
# create date columns
dates <- grep("^\\d+", names(selfadmin), value = T) # use these columns to make date columns # ignore the ...\\d columns
dates <- as.character(as.POSIXct(as.numeric(dates) * (60*60*24), origin="1899-12-30", tz="UTC", format="%Y-%m-%d")) # convert Excel character into dates

setnames(selfadmin, toupper(as.character(selfadmin[1,]) )) # now that dates are moved into separate vector, remove from the column names 
setnames(selfadmin,  mgsub::mgsub(names(selfadmin), c("PR1", "PR2", "^T.+1$", "^T.+2$", "^T.+3$", "^T.+4$", "SHA.*?([0-9]{2})$", "LGA.*?([0-9]{2})$"), c("PR01", "PR02", "PR03_T01", "PR04_T02", "PR05_T03", "PR06_T04", "SHA\\1", "LGA\\1")))
names(selfadmin)[1:2] <- c("RAT", "RFID")
selfadmin <- selfadmin[-1,]
selfadmin <- remove_empty(selfadmin, "cols") # janitor::remove_empty_cols() deprecated
# selfadmin <- selfadmin[1:grep("average", selfadmin$RAT)[1],] # subset only the rewards table

nm <- names(selfadmin)[-c(1:2)] # make date columns for this vector of exp names  ## MISSING THE RFID COLUMN SO THAT IS WHY [-1] INSTEAD OF [-c(1:2)] 
nm1 <- paste("date", nm, sep = "_") # make these date columns
selfadmin[ , ( nm1 ) := lapply( dates, function(x) rep(x, each = .N) ) ] # make the date columns 

#extract comments
comments_df <- selfadmin[which(selfadmin$RAT %in% c("COMMENT", "CONFLICT", "RESOLUTION"))]
comments_df <- comments_df %>% select(-matches("RFID|date")) %>% t() %>% 
  as.data.frame() %>% 
  rownames_to_column()
setnames(comments_df, append("EXP", comments_df[1, 2:3] %>% t() %>% unlist() %>% as.character) %>% tolower)
comments_df <- comments_df[-1,]
# selfadmin <- selfadmin[!which(selfadmin$RAT %in% c("COMMENT", "CONFLICT", "RESOLUTION"))]

selfadmin_exps <- grep("REWARDS|ACTIVE|INACTIVE|PR$", selfadmin$RAT)
selfadmin_split <- split(selfadmin, cumsum(1:nrow(selfadmin) %in% selfadmin_exps))
names(selfadmin_split) <- lapply(selfadmin_split, function(x){ x$RAT %>% head(1)}) %>% unlist() %>% as.character()
selfadmin_split <- lapply(selfadmin_split, function(x){ x %>% dplyr::filter(grepl("^[MF]\\d+", RAT))})
selfadmin_df <- selfadmin_split %>% rbindlist(idcol = "measurement") %>% dplyr::filter(measurement != "COMMENT")
selfadmin_rewards_cohort7 <- selfadmin_df %>% dplyr::filter(measurement == "REWARDS")
selfadmin_xl_cohort7 <- selfadmin_df



########################
# COHORT 8
########################
## Cancelled bc of COVID 


########################
# COHORT 9
########################

filename <- cohortfiles_xl_c01_11[7]
selfadmin <- u01.importxlsx(filename)[[1]] %>%
  as.data.table
selfadmin[ selfadmin == "n/a" ] <- NA # change all character n/a to actual NA

# set correct column names 
# create date columns
dates <- grep("^\\d+", names(selfadmin), value = T) # use these columns to make date columns # ignore the ...\\d columns
dates <- as.character(as.POSIXct(as.numeric(dates) * (60*60*24), origin="1899-12-30", tz="UTC", format="%Y-%m-%d")) # convert Excel character into dates

setnames(selfadmin, toupper(as.character(selfadmin[1,]) )) # now that dates are moved into separate vector, remove from the column names 
setnames(selfadmin,  mgsub::mgsub(names(selfadmin), c("PR1", "PR2", "^T.+1$", "^T.+2$", "^T.+3$", "^T.+4$", "SHA.*?([0-9]{2})$", "LGA.*?([0-9]{2})$"), c("PR01", "PR02", "PR03_T01", "PR04_T02", "PR05_T03", "PR06_T04", "SHA\\1", "LGA\\1")))
names(selfadmin)[1:2] <- c("RAT", "RFID")
selfadmin <- selfadmin[-1,]
selfadmin <- remove_empty(selfadmin, "cols") # janitor::remove_empty_cols() deprecated
# selfadmin <- selfadmin[1:grep("average", selfadmin$RAT)[1],] # subset only the rewards table

nm <- names(selfadmin)[-c(1:2)] # make date columns for this vector of exp names  ## MISSING THE RFID COLUMN SO THAT IS WHY [-1] INSTEAD OF [-c(1:2)] 
nm1 <- paste("date", nm, sep = "_") # make these date columns
selfadmin[ , ( nm1 ) := lapply( dates, function(x) rep(x, each = .N) ) ] # make the date columns 

#extract comments
comments_df <- selfadmin[which(selfadmin$RAT %in% c("COMMENT", "CONFLICT", "RESOLUTION"))]
comments_df <- comments_df %>% select(-matches("RFID|date")) %>% t() %>% 
  as.data.frame() %>% 
  rownames_to_column()
setnames(comments_df, append("EXP", comments_df[1, 2:3] %>% t() %>% unlist() %>% as.character) %>% tolower)
comments_df <- comments_df[-1,]
# selfadmin <- selfadmin[!which(selfadmin$RAT %in% c("COMMENT", "CONFLICT", "RESOLUTION"))]

selfadmin_exps <- grep("REWARDS|ACTIVE|INACTIVE|PR$", selfadmin$RAT)
selfadmin_split <- split(selfadmin, cumsum(1:nrow(selfadmin) %in% selfadmin_exps))
names(selfadmin_split) <- lapply(selfadmin_split, function(x){ x$RAT %>% head(1)}) %>% unlist() %>% as.character()
selfadmin_split <- lapply(selfadmin_split, function(x){ x %>% dplyr::filter(grepl("^[MF]\\d+", RAT))})
selfadmin_df <- selfadmin_split %>% rbindlist(idcol = "measurement") %>% dplyr::filter(measurement != "COMMENT")
selfadmin_rewards_cohort9 <- selfadmin_df %>% dplyr::filter(measurement == "REWARDS")
selfadmin_xl_cohort9 <- selfadmin_df



########################
# COHORT 10
########################

filename <- cohortfiles_xl_c01_11[8]
selfadmin <- u01.importxlsx(filename)[[1]] %>%
  as.data.table
selfadmin[ selfadmin == "n/a" ] <- NA # change all character n/a to actual NA

# set correct column names 
# create date columns
dates <- grep("^\\d+", names(selfadmin), value = T) # use these columns to make date columns # ignore the ...\\d columns
dates <- as.character(as.POSIXct(as.numeric(dates) * (60*60*24), origin="1899-12-30", tz="UTC", format="%Y-%m-%d")) # convert Excel character into dates

setnames(selfadmin, toupper(as.character(selfadmin[1,]) )) # now that dates are moved into separate vector, remove from the column names 
setnames(selfadmin,  mgsub::mgsub(names(selfadmin), c("PR1", "PR2", "^T.+1$", "^T.+2$", "^T.+3$", "^T.+4$", "SHA.*?([0-9]{2})$", "LGA.*?([0-9]{2})$"), c("PR01", "PR02", "PR03_T01", "PR04_T02", "PR05_T03", "PR06_T04", "SHA\\1", "LGA\\1")))
names(selfadmin)[1:2] <- c("RAT", "RFID")
selfadmin <- selfadmin[-1,]
selfadmin <- remove_empty(selfadmin, "cols") # janitor::remove_empty_cols() deprecated
# selfadmin <- selfadmin[1:grep("average", selfadmin$RAT)[1],] # subset only the rewards table

nm <- names(selfadmin)[-c(1:2)] # make date columns for this vector of exp names  ## MISSING THE RFID COLUMN SO THAT IS WHY [-1] INSTEAD OF [-c(1:2)] 
nm1 <- paste("date", nm, sep = "_") # make these date columns
selfadmin[ , ( nm1 ) := lapply( dates, function(x) rep(x, each = .N) ) ] # make the date columns 

#extract comments
comments_df <- selfadmin[which(selfadmin$RAT %in% c("COMMENT", "CONFLICT", "RESOLUTION"))]
comments_df <- comments_df %>% select(-matches("RFID|date")) %>% t() %>% 
  as.data.frame() %>% 
  rownames_to_column()
setnames(comments_df, append("EXP", comments_df[1, 2:3] %>% t() %>% unlist() %>% as.character) %>% tolower)
comments_df <- comments_df[-1,]
# selfadmin <- selfadmin[!which(selfadmin$RAT %in% c("COMMENT", "CONFLICT", "RESOLUTION"))]

selfadmin_exps <- grep("REWARDS|ACTIVE|INACTIVE|PR$", selfadmin$RAT)
selfadmin_split <- split(selfadmin, cumsum(1:nrow(selfadmin) %in% selfadmin_exps))
names(selfadmin_split) <- lapply(selfadmin_split, function(x){ x$RAT %>% head(1)}) %>% unlist() %>% as.character()
selfadmin_split <- lapply(selfadmin_split, function(x){ x %>% dplyr::filter(grepl("^[MF]\\d+", RAT))})
selfadmin_df <- selfadmin_split %>% rbindlist(idcol = "measurement") %>% dplyr::filter(measurement != "COMMENT")
selfadmin_rewards_cohort10 <- selfadmin_df %>% dplyr::filter(measurement == "REWARDS")
selfadmin_xl_cohort10 <- selfadmin_df



########################
# COHORT 11
########################

filename <- cohortfiles_xl_c01_11[8]
selfadmin <- u01.importxlsx(filename)[[1]] %>%
  as.data.table
selfadmin[ selfadmin == "n/a" ] <- NA # change all character n/a to actual NA

# set correct column names 
# create date columns
dates <- grep("^\\d+", names(selfadmin), value = T) # use these columns to make date columns # ignore the ...\\d columns
dates <- as.character(as.POSIXct(as.numeric(dates) * (60*60*24), origin="1899-12-30", tz="UTC", format="%Y-%m-%d")) # convert Excel character into dates

setnames(selfadmin, toupper(as.character(selfadmin[1,]) )) # now that dates are moved into separate vector, remove from the column names 
setnames(selfadmin,  mgsub::mgsub(names(selfadmin), c("PR1", "PR2", "^T.+1$", "^T.+2$", "^T.+3$", "^T.+4$", "SHA.*?([0-9]{2})$", "LGA.*?([0-9]{2})$"), c("PR01", "PR02", "PR03_T01", "PR04_T02", "PR05_T03", "PR06_T04", "SHA\\1", "LGA\\1")))
names(selfadmin)[1:2] <- c("RAT", "RFID")
selfadmin <- selfadmin[-1,]
selfadmin <- remove_empty(selfadmin, "cols") # janitor::remove_empty_cols() deprecated
# selfadmin <- selfadmin[1:grep("average", selfadmin$RAT)[1],] # subset only the rewards table

nm <- names(selfadmin)[-c(1:2)] # make date columns for this vector of exp names  ## MISSING THE RFID COLUMN SO THAT IS WHY [-1] INSTEAD OF [-c(1:2)] 
nm1 <- paste("date", nm, sep = "_") # make these date columns
selfadmin[ , ( nm1 ) := lapply( dates, function(x) rep(x, each = .N) ) ] # make the date columns 

#extract comments
comments_df <- selfadmin[which(selfadmin$RAT %in% c("COMMENT", "CONFLICT", "RESOLUTION"))]
comments_df <- comments_df %>% select(-matches("RFID|date")) %>% t() %>% 
  as.data.frame() %>% 
  rownames_to_column()
setnames(comments_df, append("EXP", comments_df[1, 2:3] %>% t() %>% unlist() %>% as.character) %>% tolower)
comments_df <- comments_df[-1,]
# selfadmin <- selfadmin[!which(selfadmin$RAT %in% c("COMMENT", "CONFLICT", "RESOLUTION"))]

selfadmin_exps <- grep("REWARDS|ACTIVE|INACTIVE|PR$", selfadmin$RAT)
selfadmin_split <- split(selfadmin, cumsum(1:nrow(selfadmin) %in% selfadmin_exps))
names(selfadmin_split) <- lapply(selfadmin_split, function(x){ x$RAT %>% head(1)}) %>% unlist() %>% as.character()
selfadmin_split <- lapply(selfadmin_split, function(x){ x %>% dplyr::filter(grepl("^[MF]\\d+", RAT))})
selfadmin_df <- selfadmin_split %>% rbindlist(idcol = "measurement") %>% dplyr::filter(measurement != "COMMENT")
selfadmin_rewards_cohort10 <- selfadmin_df %>% dplyr::filter(measurement == "REWARDS")
selfadmin_xl_cohort10 <- selfadmin_df


########################
# COHORT 12
########################

## Cancelled because of COVID 



