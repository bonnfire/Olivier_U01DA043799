## extract raw 
setwd("~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS")

# after cohort 5, there are only new files

## USEFUL FUNCTIONS

# FOR ~NEW~ DIRECTORIES
# #extract names to be assigned for various tables later
process_subjects_new <- function(x){
  
  read_subjects_new <- function(x){
    subjects <- fread(paste0("awk '/Subject/{print $2}' ", "'", x, "'"),fill = T,header=F)
    subjects$filename <- x
    return(subjects)
  }

  read_meta_subjects_new <- function(x){
    date_time <- fread(paste0("awk '/Start/{print $3}' ", "'", x, "'"),fill = T,header=F)
    return(date_time)
  }
  date_time <- lapply(x, read_meta_subjects_new) %>% rbindlist() 
  
  date <- date_time %>% dplyr::filter(grepl("/", V1)) %>% unlist() %>% as.character()
  time <- date_time %>% dplyr::filter(grepl(":", V1)) %>% unlist() %>% as.character()
  date_time_bind <- data.frame(date = date, time = time) %>% 
    mutate(date_time = paste0(date, "_", time))
  
  names_sha_append <- lapply(x, read_subjects_new) %>% rbindlist() %>% rename("labanimalid"="V1") %>%
    cbind(., date_time_bind) %>% 
    mutate(labanimalid = paste0(str_extract(toupper(labanimalid), "[MF]\\d{1,3}"), "_",
                                str_extract(filename, "C\\d+"), "_",
                                sub('.*HS', '', toupper(filename)), "_",
                                sub(".*/.*/.*/", '', filename), "_",
                                date_time)) %>%  # subject id, cohort, experiment, file/location perhaps
  select(-c(date_time, date, time))
  
  return(names_sha_append)
  
}

# bc some animals don't have id's under subject, mistake in manual entry; extract subject information from experiment or group? 

# #notice system call, doing it recurisvely for the entire directory; for these cases, the group number is not 0 and the subject is found elsewhere.  # 12/18 only lga
groups_files <- system("grep -ir \"Group:\" | grep -v \"Group: 0\"", intern = TRUE) %>% 
  gsub("\r", "", .) %>% 
  as.data.frame() %>% 
  rename("filename" = ".") %>% 
  separate(filename, c("filename", "labanimalid"), sep = ":", extra = "merge") %>%  # only split specified number of times
  # mutate(labanimalid = paste0(str_match(labanimalid, "[FM]\\d{1,3}"), "_", str_extract(filename, "C\\d+"), "_", sub('.*HS', '', toupper(filename)), "_", sub(".*/.*/.*/", '', filename) ),
  #        comment = "labanimalid extracted from group in raw files") %>%
  mutate(labanimalid = paste0(str_match(labanimalid, "[FM]\\d{1,3}"), "_", str_extract(filename, "C\\d+"), "_", sub('.*HS', '', toupper(filename)), "_", sub(".*/.*/.*/", '', filename) ),
         filename = paste0("./", filename))

# FOR ~NEW~ DIRECTORIES
# extracting left/right responses/timestamps 
# create the dataframe with vector in it
# check that the column we are removing ends with 5 or 0 and then remove
# write a df containing the fread statements and function for extracting the different dfs
read_fread_new <- function(x, varname){
  
  fread_statements <- data.frame(varname = c("leftresponses", "rightresponses", "rewards", "lefttimestamps", "righttimestamps", "rewardstimestamps"),
                                 statement = c("awk '/L:/{flag=1;next}/R:/{flag=0}flag' ",
                                               "awk '/R:/{flag=1;next}/U:/{flag=0}flag' ",
                                               "awk '/W:/{flag=1;next}/Y:/{flag=0}flag' ", 
                                               "awk '/U:/{flag=1;next}/V:/{flag=0}flag' ",
                                               "awk '/Y:/{flag=1;next}/^$/{flag=0}flag' ",
                                               "awk '/V:/{flag=1;next}/W:/{flag=0}flag' "))
  statement <- fread_statements[which(fread_statements$varname == varname),]$statement
  rawdata <- fread(paste0(statement, "'", x, "'"), fill = T)
  data_indices <- grep("^0:$", rawdata$V1)
  split_data <- split(rawdata, cumsum(1:nrow(rawdata) %in% data_indices))

  keepzeroes <- c("leftresponses", "rightresponses", "rewards") # preserve bin sequences
  
  if(varname %in% keepzeroes){
    processeddata <- lapply(split_data, function(x){
      indexremoved <- x[,-1]
      processeddata_df <- data.frame(counts = as.vector(t(data.matrix(indexremoved)))) %>% # transpose to get by row
        mutate(bin = ifelse(row_number() == 1, "total", as.character(row_number() - 1)))
      return(processeddata_df)
    })
  }
  else{
    processeddata <- lapply(split_data, function(x){
      indexremoved <- x[,-1]
      nonzerorows <- indexremoved[rowSums(indexremoved) > 0, ] # remove excessively trailing 0's 
      processeddata_df <- data.frame(timestamps = as.vector(t(data.matrix(nonzerorows)))) # transpose to get by row
      if(any(processeddata_df$timestamps > 7500)){
        processeddata_df %<>% 
          mutate(bin = cut(timestamps, breaks=seq(from = 0, length.out = 73, by = 300), right = T, labels = seq(from = 1, to = 72, by =1))) %<>% 
          dplyr::filter(timestamps != 0)
      }
      else{
        processeddata_df %<>% 
          mutate(bin = cut(timestamps, breaks=seq(from = 0, length.out = 25, by = 300), right = T, labels = seq(from = 1, to = 24, by =1))) %<>% 
          dplyr::filter(timestamps != 0)
      }
      return(processeddata_df)
    }) 
  }
  

  return(processeddata)
}

# FOR ~OLD~ DIRECTORIES

## know how many subjects to expect in each filename
## find . -name "SHA" -exec grep -ira1 "NumberOfSubjects" {} +

process_subjects_old <- function(x){
  
  read_subjects_old <- function(x){
    subject_old <- fread(paste0("grep -iEA1 --no-group-separator \"ratnumber|boxnumber\" ", "'", x, "'", "| grep -vE \"Rat|Box\""), header = F)
    subject_old$filename <- x 
    return(subject_old)
  }
  
  subjects_old <- lapply(x, read_subjects_old) %>% rbindlist() %>% 
    rename("labanimalid" = "V1") %>% 
    mutate(box_id = ifelse((row_number() %% 2) == 0, "labanimalid", "box"))
  
  box <- subjects_old %>% dplyr::filter(box_id == "box") %>% select(labanimalid) %>% unlist() %>% as.character()
  labanimalid <- subjects_old %>% dplyr::filter(box_id == "labanimalid") %>% select(labanimalid) %>% unlist() %>% as.character()
  filename <- subjects_old %>% dplyr::filter(box_id == "box") %>% select(filename) %>% unlist() %>% as.character()
  
  box_id_bind <- data.frame(box = box, 
                            labanimalid = labanimalid, 
                            filename = filename) %>% 
    mutate_all(as.character) %>% 
    mutate(labanimalid = replace(labanimalid, labanimalid=="999", "F000")) %>% # create placeholder for the problematic cases
    mutate(date = str_extract(filename, "\\d{8}(-\\d+)?") %>% lubridate::ymd(),
           cohort = str_extract(filename, "C\\d+"), 
           # cohort = gsub(".*C[0]?(\\d)+/.*", "cohort\\1", filename),
           exp = sub("-.*", "", sub(".*HS([^.]+)[-].*", "\\1", toupper(filename)))) %>% 
    merge(., cohorts_exp_date, by = c("cohort", "exp")) %>% 
    mutate(valid = ifelse(date == excel_date, "valid", "invalid")) %>%
    mutate(labanimalid = paste0(str_match(toupper(labanimalid), "[FM]\\d{1,3}"), "_",
                                box, "_",
                                str_extract(filename, "C\\d+"), "_",
                                sub("-.*", "", sub(".*HS([^.]+)[-].*", "\\1", toupper(filename))), "_",
                                sub("C.*", "", sub(".*/.*/.*/.*/", "", filename)), "_",
                                str_extract(filename, "\\d{8}"), "_",
                                valid)) %>%  # subject id, box, cohort, experiment, computer, date
    select(one_of("labanimalid", "filename"))
  # 
  # %>% dplyr::filter(!grepl("^NA", labanimalid))
  
  return(box_id_bind)
}

# FOR ~OLD~ DIRECTORIES
read_fread_old <- function(x, varname){
  fread_old_statements <- data.frame(varname = c("leftresponses", "rightresponses", "rewards"),
                                 statement = c("awk '/^BinsInActiveResponses/{flag=1;next}/endl/{flag=0}flag' ",
                                               "awk '/^ResponsesActBins/{flag=1;next}/endl/{flag=0}flag' ",
                                               "awk '/totalRewards/{flag=1;next}/TotalResponses/{flag=0}flag' ")) 
                                               # "awk '/BinRewards/{flag=1;next}/endl/{flag=0}flag' "))  #### 	 In=L Act=R  Rew=W InTS=U	ActTS=Y  RewTS=V  RewIRI=Z 	
  statement <- fread_old_statements[which(fread_old_statements$varname == varname),]$statement
  rawdata <- fread(paste0(statement, "'", x, "'"), fill = T, header = F)
  rawdata$filename <- x
  return(rawdata)
}

## all ts + iri labels (4) are appended w 2000, count labels (3) have just 24

read_iri_old <- function(x){
  iri <- fread(paste0("awk '/^IRI\\r/{flag=1;next}/endl/{flag=0}flag' ", "'", x, "'"), fill = T, header = F)
  iri$filename <- x
  names(iri)[1] <- "iritime"
  
  iricode <- fread(paste0("awk '/^IRICode\\r/{flag=1;next}/endl/{flag=0}flag' ", "'", x, "'"), fill = T, header = F)
  iricode$filename <- x
  names(iricode)[1] <- "iricode"
  

  rawdata <- cbind(iri, iricode)
  rawdata <- rawdata[,-2]
  
  return(rawdata)
}

convert_iri_matrix_to_df <- function(x){
  if((x[2,1] == nrow(x) -2) == T){
    x <- x[-c(1:2),]
    x <- x %>% 
      mutate(iritime = as.numeric(iritime),
             iritime = iritime/100)
  }
  
  timer = x$iritime[1]
  
  for(i in 1:nrow(x)){
    timer = timer + x$iritime[i]
    if(x$iricode[i] %in% c(1,4)){
      x$activeTS[i] = timer
    } else x$activeTS[i] = NA 
    
    if(x$iricode[i] %in% c(3,5)){
      x$inactiveTS[i] = timer
    } else x$inactiveTS[i] = NA 
    
    
    if(x$iricode[i] == 1){
      x$rewardTS[i] = timer
    } else x$rewardTS[i] = NA 
    
  }
  
  return(x)}
  
# FOR ~NEW~ DIRECTORIES
# NON EXPERIMENT SPECIFIC 
## TO VALIDATE ENTRIES  
## extract date and start time/end time to determine valid sessions
read_date_time_subject <- system("grep -a7r --no-group-separator \"Start Date: \" . | grep -E \"(Start Date|End|Subject|Box|Start Time|End Time):\"", intern = T)
read_date_time_subject <- gsub("\\r", "", read_date_time_subject)
read_date_time_subject <- read_date_time_subject[!grepl("/LGA/:", read_date_time_subject)] # remove the duplicate file # removed on 1/27

date_time_subject <- data.frame(labanimalid = gsub(".*Subject: ", "", grep("Subject", read_date_time_subject, value = T)) %>% toupper,
                                   cohort = str_match(grep("Subject", read_date_time_subject, value = T), "C\\d{2}") %>% unlist() %>% as.character(),  
                                   exp = toupper(sub('.*HS', '', grep("Subject", read_date_time_subject, value = T) %>% gsub("-Subject.*", "", .))),
                                   start_date = gsub(".*Start Date: ", "", grep("Start Date:", read_date_time_subject, value = T)),
                                   box = gsub(".*Box: ", "", grep("Box", read_date_time_subject, value = T)),
                                   start_time = gsub(".*Start Time: ", "", grep("Start Time", read_date_time_subject, value = T)),
                                   end_time = gsub(".*End Time: ", "", grep("End Time", read_date_time_subject, value = T)),
                                   filename = sub(".*/.*/.*/", '', grep("Subject", read_date_time_subject, value = T)) %>% gsub("-Subject.*", "", .),
                                   directory = str_match(grep("Subject", read_date_time_subject, value = T) %>% gsub("-Subject.*", "", .), "New_medassociates|Old") %>% unlist() %>% as.character()
)

date_time_subject_mut <- date_time_subject %>% 
  mutate(start_date = lubridate::mdy(format(as.Date(start_date, "%m/%d/%y"), "%m/%d/20%y")),
         start_time = chron::chron(times = start_time),
         end_time = chron::chron(times = end_time), 
         experiment_duration = end_time - start_time,
         experiment_duration = 60 * 24 * as.numeric(chron::times(experiment_duration))
         ) %>% 
  mutate_if(is.factor, as.character) %>% 
  mutate(labanimalid = replace(labanimalid, labanimalid == "M7678", "M768"),
         labanimalid = replace(labanimalid, labanimalid == "X", "F507"),
         labanimalid = replace(labanimalid, labanimalid == "717", "F717")) %>% # verified by box
  mutate(labanimalid = if_else(grepl("^C0", labanimalid), str_match(labanimalid,"[FM]\\d{1,3}" %>% unlist() %>% as.character()), labanimalid %>% as.character()))

## problems in being too lax in accepting all forms of subjects 
# gsub(".*Subject: ", "", grep("Subject", read_date_time_subject, value = T)) %>% toupper %>% table() # before processing
# date_time_subject_mut[str_detect(date_time_subject_mut$labanimalid, "^(M|F)\\d{4}", negate = F),]
# date_time_subject_df %>% subset(labanimalid == "0") %>% group_by(filename) %>% dplyr::filter(n() > 5) # more than 5 is most likely a broken file but less than five is most likely a dead rat? 
date_time_subject_mut$labanimalid %>% table()

#trying to fix the subject 0
subject0 <- date_time_subject_mut %>% split(., .$cohort) %>% lapply(., function(x){
  x <- x %>% 
    dplyr::filter(!grepl("SHOCK", exp)) %>% 
    mutate(room = ifelse(grepl("[[:alnum:]]+C\\d{2}HS", filename), gsub("C\\d{2}HS.*", "", filename), NA)) %>% 
    arrange(room, as.numeric(box)) %>% 
    dplyr::filter(!grepl("[MF]", labanimalid)|lead(!grepl("[MF]", labanimalid))|lag(!grepl("[MF]", labanimalid))) %>% 
    mutate(dbcomment = ifelse(!grepl("[MF]", labanimalid), "box info used to fill labanimalid", NA)) %>% 
    group_by(box) %>% mutate(labanimalid = labanimalid[grepl("[MF]", labanimalid)][1],
                             labanimalid = replace(labanimalid, cohort == "C04"&box == "2", "F402"),
                             labanimalid = replace(labanimalid, cohort == "C04"&box == "4", "F404")) %>%  # spot checking for deaths
    arrange(labanimalid, start_date) 
    return(x)
  }) %>%  rbindlist(., idcol = "cohort")

# replace rbindlist... with openxlsx::write.xlsx(., "labanimalid_assign_bybox.xlsx") to create the excel sheets that I sent to their lab 
subject0[[3]] <- NULL
subject0 %>% openxlsx::write.xlsx(., "labanimalid_assign_bybox.xlsx")

# trying to fix subject 0 for shock files 
date_time_subject_mut %>% split(., .$cohort) %>% lapply(., function(x){
  x <- x %>% 
    dplyr::filter(grepl("SHOCK", exp)) %>% 
    mutate(room = ifelse(grepl("[[:alnum:]]+C\\d{2}HS", filename), gsub("C\\d{2}HS.*", "", filename), NA)) %>% 
    dplyr::filter(!grepl("[MF]", labanimalid)|lead(!grepl("[MF]", labanimalid))|lag(!grepl("[MF]", labanimalid))) %>% 
    arrange(room, as.numeric(box))
  return(x)
}) 

# remove labanimalid0 or blank subset from original df and then insert the corrected ones (keep the dbcomment variable)
date_time_subject_df <- date_time_subject_mut %>% 
  dplyr::filter(grepl("[MF]", labanimalid)) %>% 
  plyr::rbind.fill(., subject0) %>% # rbind with the added function of creating an NA column for nonmatching columns bw dfs A and B
  arrange(cohort, start_date, as.numeric(box)) %>% 
  distinct() %>%  # needed bc otherwise subject0 will "double count" the "reference" rows
  mutate(exp = gsub("-.*", "", exp),
         exp = replace(exp, as.numeric(str_extract(cohort, "\\d+")) > 6&grepl("SHOCK", exp), "SHOCK03")) # for all cohorts later than cohort 6, they only use one shock value, but we can change it to shock03 so that we can have uniformity 
## XX Note: remove all SHOCK03 and Pre Shock - Olivier (from meeting 1/24)


# waiting on their response for these cases (trying to assign labanimalid [MF]\\d{4,})
## date_time_subject_df %>% arrange(cohort, as.numeric(box)) %>% dplyr::filter(!grepl("[MF]\\d{1,3}(?!\\d+?)", labanimalid, perl = T )|lead(!grepl("[MF]\\d{1,3}(?!\\d+?)", labanimalid, perl = T ))|lag(!grepl("[MF]\\d{1,3}(?!\\d+?)", labanimalid, perl = T )))
## OR  date_time_subject_df %>% dplyr::filter(grepl("[MF]\\d{4,}", labanimalid, perl = T ))

## date_time_subject_df_comp %>% dplyr::filter(labanimalid == "F16", exp == "SHA06", box %in% c("8", "16"))

## PICK BACK UP 
# before merging with excel dates
# include more dbcomments 
# fix strange filenames -2, -3
date_time_subject_df$exp %>% table()
date_time_subject_df %>% 
  dplyr::filter(filename %in% c(grep("-", date_time_subject_df$filename, value = T), 
                                gsub("-.*", "", grep("-", date_time_subject_df$filename, value = T)))) %>% 
  arrange(labanimalid) %>% group_by(labanimalid, exp) %>% dplyr::filter(n()>1)
## interesting... keep this here just in case it comes up (fixed, needed more conditions in the replace logic)
date_time_subject_df %>% subset(as.numeric(str_extract(cohort, "\\d+")) > 6) %>% select(exp, cohort) %>% table()
## PICK BACK UP 
# date_time_subject_df <- date_time_subject_df %>% 
#   mutate(dbcomment = replace())

# include correct dates as another check (dates extracted from CREATE_DATABASESTRUCTURE allcohorts2 object)
# allcohorts2 %>% select(matches("date|cohort")) %>% distinct()
# reformat exported excel object to prepare for merge and check with in file dates ((wrong dates should be noted and possbily removed))

cohorts_exp_date <- allcohorts2 %>% 
  mutate(date_lga19 = replace(date_lga19, cohort == "cohort5", lubridate::ymd("2018-09-19")),
         date_sha02 = replace(date_sha02, cohort == "cohort5", lubridate::ymd("2018-07-31"))) %>% select(matches("date|cohort")) %>% distinct() %>% # for record keeping, make sure to make this change on the actual excel! 
gather(v, value, date_sha01:date_lga23) %>% 
  separate(v, c("date", "exp")) %>% 
  arrange(cohort) %>% 
  select(-date) %>% 
  mutate(cohort = paste0("C", str_pad(gsub("COHORT", "", toupper(cohort)),  2, "left","0")),
         exp = toupper(exp),
         value = lubridate::ymd(value)) %>% 
  rename("excel_date" = "value")

date_time_subject_df_comp <- left_join(date_time_subject_df, cohorts_exp_date, by = c("cohort", "exp")) %>%
  mutate(valid = case_when(
    grepl("SHOCK", exp) & experiment_duration > 58 & excel_date == start_date ~ "yes",
    grepl("SHA", exp) & experiment_duration > 115 & excel_date == start_date~ "yes",
    grepl("LGA", exp) & experiment_duration > 355 & excel_date == start_date~ "yes",
    grepl("PR", exp) & experiment_duration < 60 & excel_date == start_date~ "yes"),
    valid = replace(valid, is.na(valid), "no")
  ) # 9808 for cohorts C01-C09 (no C06) ## change the minimum times - Olivier (from 1/24 meeting)

# date_time_subject_df_comp %>% dplyr::filter(valid == "yes", labanimalid != 0) %>% group_by(labanimalid, exp) %>% dplyr::filter(n() > 1)
# date_time_subject_df_comp %>% dplyr::filter(valid == "yes", labanimalid != 0, is.na(dbcomment)) %>% group_by(labanimalid, exp) %>% dplyr::filter(n() > 1)

# View(date_time_subject_df_comp  %>% dplyr::filter(valid == "yes") %>% select(cohort, exp, start_date) %>% distinct() %>% spread(exp, start_date))# # files that have different dates 
# View(date_time_subject_df %>% select(cohort, exp, start_date) %>% distinct() %>% group_by(cohort, exp) %>% dplyr::filter(n() > 1))
# ## SECTION OFF R SCRIPT TO DIFFERENTIATE THESE FILES, XX ALSO SECTION OFF BASED ON PR AND FR (DON'T FORGET THE CONSTRAINTS ON THE PR)
# date_time_subject_df %>% 
#   dplyr::filter(!(cohort == "C07"&exp=="SHA07"&start_date==as.Date("2019-01-25")), 
#                 !(cohort == "C07"&exp=="LGA04"&start_date==as.Date("2019-02-11")), 
#                 !(cohort == "C08"&exp=="LGA08"&start_date==as.Date("2019-07-12"))) %>% 
#   select(cohort, exp, start_date) %>% distinct()%>% group_by(cohort, exp) %>% dplyr::filter(n() > 1)
# 
# View(date_time_subject_df %>% 
#   dplyr::filter(!(cohort == "C07"&exp=="SHA07"&start_date==as.Date("2019-01-25")), 
#                 !(cohort == "C07"&exp=="LGA04"&start_date==as.Date("2019-02-11")), 
#                 !(cohort == "C08"&exp=="LGA08"&start_date==as.Date("2019-07-12"))) %>% 
#   select(cohort, exp, start_date) %>% distinct()%>% spread(exp, start_date)
#   )
# do any valid files have rats that are assigned to the wrong boxes (no)
# date_time_subject_df_comp %>% dplyr::filter(subject != 0, valid == "yes") %>% select(subject, box, exp, valid) %>% group_by(subject, box, exp) %>% dplyr::filter(n() > 1) %>% arrange(subject) 

################################
########## SHA #################
################################

###### NEW FILES ##############
## get the filename and subject names 
setwd("~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS")
sha_new_files <- grep(grep(list.files(path = ".", recursive = T, full.names = T), pattern = ".*txt", inv = T, value = T), pattern = ".*SHA", value = T) # 178 files
# get subjects 
sha_subjects_new <- process_subjects_new(sha_new_files) ## sha_subjects_new %>% dplyr::filter(!grepl("^[MF]", labanimalid)) %>% head # will need to use the comp df to name sessions
# get raw data and name raw data with subjects 
rewards_sha_new <- lapply(sha_new_files, read_fread_new, "rewards") %>% unlist(recursive = F)
length(rewards_sha_new)==nrow(sha_subjects_new) # check if they match in length
names(rewards_sha_new) <- sha_subjects_new$labanimalid # 2783 matches! 

rewards_sha_df <- rewards_sha_new %>%
  rbindlist(fill = T, idcol = "labanimalid") %>% 
  subset(bin == "total") %>%
  separate(labanimalid, c("labanimalid", "cohort", "exp", "filename", "start_date", "start_time"), "_") %>% 
  mutate(start_date = lubridate::mdy(start_date),
         start_time = chron::times(start_time)) %>% 
  select(-c("labanimalid", "exp")) %>% # bc the labanimalid and exp are not processed as in the comp df
  merge(., date_time_subject_df_comp %>% subset(valid == "yes"&grepl("SHA", exp)), by = c("cohort", "filename", "start_date", "start_time")) %>% # dropped from 87395 to 50205 files
  distinct()  
   # only 2202


# subset(filename == "C01HSSHA06"&bin == "total"&start_time==chron::times("07:28:22"))
# date_time_subject_df_comp %>% subset(valid == "yes"&exp=="SHA06"&labanimalid=="F10")



# each animal should only have one value for each exp
rewards_sha_df %>% group_by(labanimalid, exp) %>% dplyr::filter(n() > 1) %>% arrange(labanimalid) %>% View()

rewards_sha_totalsession_with_date_df <- rewards_sha_df %>% select(-c(bin, filename)) %>% spread(exp, counts) %>% arrange(labanimalid, start_date, start_time)
rewards_sha_totalsession_df <- rewards_sha_df %>% 
  select(-c("bin", "filename", "start_date", "start_time", "end_time", "excel_date", "experiment_duration")) %>% 
  spread(exp, counts) %>% 
  mutate_at("box", as.numeric) %>% 
  arrange(cohort, box) # cannot spread bc there are duplicate values for the experiment
# see the duplicated values
## XXXXXXXXXXXXXXXXXXX LAUREN 
rewards_sha_df %>% dplyr::filter(bin == "total", labanimalid != "NA") %>% group_by(labanimalid) %>% add_count(file_exp) %>% subset(n != 1) %>% View()
rewards_sha_bybin_df <- rewards_sha_df %>% dplyr::filter(bin != "total")

ratinfo_list_replacements_processed %>% head
## merge(WFU_Olivier_co_test_df[, c("cohort", labanimalnumber", "rfid")])

## exclude the following files - Olivier (1/24 meeting)
### -----
### -----

########### NEW SHA
# 12/16 for rsm: 
# names_sha_rsm <- names_sha %>%
#   rename("labanimalid"="V1") %>%
#   mutate(labanimalid = paste0(str_extract(toupper(labanimalid), "[MF]\\d{1,3}"), "_", str_extract(filename, "C\\d+")),
#          #file_exp = str_extract(toupper(filename), "\\D+\\d+$"))
#   file_exp = sub('.*HS', '', toupper(filename))) %>%
#   group_by(labanimalid) %>%
#   add_count(file_exp) %>%
#   ungroup()
# 
# exps_vs_id <- table(names_sha_rsm$file_exp, factor(names_sha_rsm$labanimalid, levels = unique(gtools::mixedsort(names_sha_rsm$labanimalid)))) %>% t()
# exps_vs_id <- cbind(exps_vs_id, Total = rowSums(exps_vs_id))
# # # exps_vs_id <- rbind(exps_vs_id, Total = colSums(exps_vs_id)) # this line gets rid of the row name values for the excel
# openxlsx::write.xlsx(exps_vs_id, file = "exps_vs_id.xlsx",col.names=TRUE, row.names=TRUE)

###### OLD SHA ##############
## get the filename and subject names 
sha_old_files <- grep(list.files(path = ".", recursive = T, full.names = T), pattern = ".*Old.*SHA", value = T) # 214 files
# get subjects 
sha_subjects_old <- process_subjects_old(sha_old_files) #19 failed to parse warning message #1255 
# extract the sha data
rewards_sha_old <- lapply(sha_old_files, read_fread_old, "rewards") %>% rbindlist() #1255
nrow(sha_subjects_old) == nrow(rewards_sha_old)
# bind to data
rewards_sha_old_df <- rewards_sha_old %>% 
  cbind(., sha_subjects_old[,"labanimalid"]) %>% 
  rename("rewards" = "V1",
         "labanimalid" =  "V2") %>% 
  separate(labanimalid, into = c("labanimalid", "box", "cohort", "exp", "computer", "date", "valid"), sep = "_") %>%
  mutate_if(is.character, .funs = gsub, pattern = "[[:space:]]", replacement = "") %>% 
  mutate(date = lubridate::ymd(date)) %>% 
  subset(valid == "valid") %>% 
  select(-valid)
# only needed if we extract the binned data, but the rewards now only extracts the value after TotalRewards
# rewards_sha_old <- lapply(rewards_sha_old, function(x){
#   Index <- which(x[,1]=="list")
#   if(x[(Index+1),] == 12){
#     x <- x[-(Index+1),]
#   }
#   return(x)
# }) %>% rbindlist() # use indexing to remove the 12 value if it follows list
# rewards_sha_old <- lapply(sha_old_files, read_fread_old, "rewards") %>% rbindlist()
# rewards_sha_old_indices <- grep("list", rewards_sha_old$V1)
# rewards_sha_old_split <- split(rewards_sha_old, cumsum(1:nrow(rewards_sha_old) %in% rewards_sha_old_indices))
# names(rewards_sha_old_split) <- sha_subjects_old$labanimalid ##  1255 matches! subject lines should match the length of list data
# rewards_sha_old_df <- rewards_sha_old_split %>% 
#   rbindlist(idcol = "labanimalid") %>% 
#   rename("counts" = "V1") %>% 
#   group_by(labanimalid) %>% 
#   slice(3:n()) %>% 
#   summarize(tot_rewards = sum(as.numeric(counts), na.rm = T)) %>% 
#   separate(labanimalid, into = c("labanimalid", "file_cohort", "file_exp", "computer", "file_date"), sep = "_")




# old_sha_varname = c("leftresponses", "rightresponses", "rewards")
# for(i in 1:length(old_sha_varname)){
#   sha_old <- lapply(sha_old_files, read_fread_old, old_sha_varname[i]) 
#   sha_old <- lapply(sha_old, function(x){
#     Index <- which(x[,1]=="list")
#     if(x[(Index+1),] == 12){
#       x <- x[-(Index+1),]
#     }
#     return(x)
#   }) %>% rbindlist() 
# }


## try to make into for loop
# list_indices <- grep("list", sha_old$V1)
# sha_old_split <- split(sha_old, cumsum(1:nrow(sha_old) %in% list_indices))
# names(sha_old_split) <- sha_subjects_old$labanimalid ##  1255 matches! subject lines should match the length of list data
# # find . -regex ".*/C0[0-9]/Old/SHA/[^/]*.txt" -exec grep -i 'ratnumber' {} \; | wc -l
# # find . -regextype sed -regex ".*/C0[0-9]/Old/SHA" -exec sh -c 'grep -ira1 "ratnumber" | grep -Eiv "[F|M]" | grep -vi "ratnumber"' sh {} \;
# # find . -regex ".*/C0[0-9]/Old/SHA/[^/]*.txt" -exec sh -c 'grep -ira1 "ratnumber"' sh {} \;

# sha_old_df <- sha_old_split %>% 
#   rbindlist(idcol = "labanimalid") %>%
#   rename("counts" = "V1") %>% 
#   dplyr::filter(counts != "list") %>% 
#   separate(labanimalid, into = c("labanimalid", "file_cohort", "file_exp", "computer", "file_date"), sep = "_")

# check if labanimalids are valid
sha_old_df[str_detect(sha_old_df$labanimalid, "^[MF]\\d+$", negate = T),] %>% View() # no NA labanimalid, but what is F000
sha_old_df %>% subset(labanimalid == "F000") # using the placeholder for problematic cases

# extract iri data 
iri <- lapply(sha_old_files, read_iri_old) %>% rbindlist()
iri %>%
  dplyr::filter(iritime == "list", iricode != "list") %>% 
  nrow() == 0  # should be none/TRUE! since we are using cbind, this ensrues that the df's are "synced"
iri_indices <- grep("list", iri$iritime)
sha_old_iri <- split(iri, cumsum(1:nrow(iri) %in% iri_indices))
length(sha_old_iri) == length(sha_subjects_old$labanimalid)
names(sha_old_iri) <- sha_subjects_old$labanimalid ## subject lines now match the length of list data

old_iri_df <- lapply(sha_old_iri, convert_iri_matrix_to_df) %>% rbindlist(idcol = "labanimalid") %>% 
  separate(labanimalid, c("labanimalid", "file_cohort", "file_exp", "computer",  "filedate"), "_") %>% 
  mutate(filedate = lubridate::ymd(gsub("-\\d", "", filedate))) # bc the file extension has -2 (and maybe other tails) that need to be removed before we can parse into dates
# checking valid id's
old_iri_df[str_detect(old_iri_df$labanimalid, "^[MF]\\d+$", negate = T),] %>% View()
old_iri_df %>% subset(labanimalid == "F000") 

################################
########## LGA #################
################################
lga_new_files <- grep(grep(list.files(path = ".", recursive = T, full.names = T), pattern = ".*txt", inv = T, value = T), pattern = ".*LGA", value = T) # 329 files
lga_new_files <- lga_new_files[-1] # duplicate, 329 -> 328 # confirmed on bash
lga_old_files <- grep(list.files(path = ".", recursive = T, full.names = T), pattern = ".*Old.*LGA", value = T) # 424 files


# get subjects 
lga_subjects_new <- process_subjects_new(lga_new_files)
lga_subjects_old <- process_subjects_old(lga_old_files)


################################ NEW

# process new files subjects that were lost in the group info
lga_subjects_new[which(lga_subjects_new$filename %in% groups_files$filename),]$labanimalid <- groups_files$labanimalid # replace with group info

missingsubjects_filename <- lga_subjects_new_totroubleshoot %>% 
  group_by(filename) %>% 
  summarise(occurence_seq = max(occurence_seq)) %>% 
  mutate(filename = gsub("./", "", filename, fixed = T)) %>% 
  merge(rewards_lga_new_troubleshoot_count) %>% 
  dplyr::filter(occurence_seq != count)  %>% 
  rename("num_subject" = "occurence_seq", "num_rewards_array" = "count") %>% 
  select(filename) %>% unlist() %>% as.character()%>% paste0("./", .)
r <- as.numeric(c("2546", "2561", "2592", "2607", "2622", "2637")) - 1
for(i in 1:length(missingsubjects_filename)){
  newrow <- c(paste0("M753", "_", 
                     str_extract(missingsubjects_filename[i], "C\\d+"), "_", 
                     sub('.*HS', '', toupper(missingsubjects_filename[i])), "_", 
                     sub(".*/.*/.*/", '', missingsubjects_filename[i])), missingsubjects_filename[i])
  lga_subjects_new <- rbind(lga_subjects_new[1:r[i],], newrow, lga_subjects_new[-(1:r[i]),])
}

# XX issue: lga_subjects_new %>% dplyr::filter(!grepl("^[MF]", labanimalid)) %>% head


lga_subjects_new_totroubleshoot <- lga_subjects_new %>% group_by(filename) %>% mutate(occurence_seq = row_number())
lga_subjects_new_totroubleshoot %>% ungroup()  %>% rownames_to_column() %>% subset(filename %in% missingsubjects_filename&occurence_seq == 1) %>% select(rowname) %>% unlist() %>% as.character() 


## get data

rewards_lga_new <- lapply(lga_new_files, read_fread, "rewards") %>% unlist(recursive = F)

rewards_lga_new_troubleshoot <- rewards_lga_new %>% rbindlist() %>% subset(bin == "total")

rewards_lga_new_troubleshoot_count <- data.frame(results = system("grep -r -c \"W:\" */New_medassociates/LGA*", intern = T)) 
rewards_lga_new_troubleshoot_count <- rewards_lga_new_troubleshoot_count %>% 
  separate(results, into = c("filename", "count"), sep = ":") %>% 
  mutate(count = as.numeric(count))
sum(as.numeric(rewards_lga_new_troubleshoot_count$count), na.rm = T) # 4935 

rewards_lga_new_troubleshoot_sub_count <- data.frame(results = system("grep -r -c \"Subject:\" */New_medassociates/LGA*", intern = T)) 
rewards_lga_new_troubleshoot_sub_count <- rewards_lga_new_troubleshoot_sub_count %>% 
  separate(results, into = c("filename", "count"), sep = ":") %>% 
  mutate(count = as.numeric(count))
sum(as.numeric(rewards_lga_new_troubleshoot_sub_count$count), na.rm = T) # 4935 

lga_subjects_new_totroubleshoot %>% group_by(filename) %>% summarise(occurence_seq = max(occurence_seq)) %>% mutate(filename = gsub("./", "", filename, fixed = T)) %>%  merge(rewards_lga_new_troubleshoot_count) %>% dplyr::filter(occurence_seq != count)  %>% rename("num_subject" = "occurence_seq", "num_rewards_array" = "count")

# to get the row numbers 
lga_subjects_new_totroubleshoot %>% group_by(filename) %>% summarise(occurence_seq = max(occurence_seq)) %>% mutate(filename = gsub("./", "", filename, fixed = T)) %>%  merge(rewards_lga_new_troubleshoot_count) %>% rownames_to_column() %>% subset(occurence_seq != count)  %>% rename("num_subject" = "occurence_seq", "num_rewards_array" = "count") %>%  `[[`("rowname") 


# names(rightresponses_sha) <- names_sha_append
names(rewards_lga_new) <- lga_subjects_new$labanimalid
rewards_lga_new_df <- rewards_lga_new %>% 
  rbindlist(fill = T, idcol = "labanimalid") %>% 
  separate(labanimalid, c("labanimalid", "file_cohort", "file_exp", "filename"), "_")
lga_subjects_new %>% 
  separate(labanimalid, c("labanimalid", "file_cohort", "file_exp", "filename"), "_") %>% add_count(file_cohort)



ratinfo_list_deaths_processed 



























# extract the subject id information
lga_subjects <- lapply(olivier_cocaine_files_lga, readsubjects) %>% rbindlist()
lga_subjects_use <- lga_subjects %>% 
  select(V1) %>% 
  unlist() %>% 
  as.vector() %>% 
  paste0(gsub(".*(C\\d+)HS(.*)","\\1\\2", names$filename)) %>% 
  toupper() %>%
  str_extract("F\\d+.*")
lga_subjects_use <- lga_subjects_use[!is.na(lga_subjects_use)]


# rightresponseslga01 <- read_fread(olivier_cocaine_files[[2]], "rightresponses")
olivier_cocaine_files_lga <- grep(grep(list.files(path = ".", recursive = T, full.names = T), pattern = ".*txt", inv = T, value = T), pattern = ".*LGA", value = T) # 329 files
names_lga <- lapply(olivier_cocaine_files_lga, readsubjects) %>% rbindlist()
names_lga_append <- names_lga %>% 
  select(V1) %>% 
  unlist() %>% 
  as.vector() %>% 
  paste0(gsub(".*(C\\d+)HS(.*)","\\1\\2", names_lga$filename)) %>% 
  toupper() %>%
  str_extract("F\\d+.*")
names_lga_append <- names_lga_append[!is.na(names_lga_append)]

rightresponses_lga <- lapply(olivier_cocaine_files_lga, read_fread, "rightresponses") %>% unlist(recursive = F)
names(rightresponses_lga) <- names_lga_append
rightresponses_lga <- rightresponses_lga %>% 
  rbindlist(fill = T, idcol = "labanimalid") %>% 
  mutate(file_cohort = str_extract(labanimalid, "C\\d+"), 
         file_exp = str_extract(labanimalid, "\\D+\\d+$"), 
         labanimalid = str_extract(labanimalid,"^\\D\\d+")) 
# %>% 
## merge(WFU_Olivier_co_test_df[, c("cohort", labanimalnumber", "rfid")])


right_time_responses <- lapply(olivier_cocaine_files_lga, read_fread, definedvars[4]) %>% unlist(recursive = F)
names(right_time_responses) <- names_append
right_time_responses[[3]]

### for loop to process of variables of interest and create objects in the glob env
definedvars <- c("leftresponses", "rightresponses", "rewards", "lefttimestamps", "righttimestamps", "rewardstimestamps")
# for(i in 1:length(definedvars)){
# definedvars_list[i] <- lapply(olivier_cocaine_files, read_fread, definedvars[i]) %>% unlist(recursive = F)
# provide names as well 
# #list2env(definedvars_list, envir = .GlobalEnv)
# }





################################
########## PR ##################
################################

pr_new_files <- grep(list.files(path = ".", recursive = T, full.names = T), pattern = ".*New.*PR/", value = T) # 53 files
pr_old_files <- grep(list.files(path = ".", recursive = T, full.names = T), pattern = ".*Old.*PR/", value = T) # 62 files

pr_subjects_new <- process_subjects_new(pr_new_files) #800 
pr_rewards_new <- lapply(pr_new_files, read_fread, "rewards")  %>% unlist(recursive = F) #785 XX MISSING REWARDS ARRAY

read_pr <- system("grep -iEa1r --no-group-separator  \"(Subject|Start Date|Start Time|H|A|G):\" */New_medassociates/PR* | grep -iE \"(Subject|Start Date|Start Time|A|G| 0):\"", intern = T) # THE REASON WHY THIS DOESN'T WORK ON OLD IS BC THE FORMAT IS DIFFERENT
read_pr <- gsub("\r", "", read_pr)

breakpoint = grep("0:", read_pr, value = T) %>% strsplit("       ") %>% sapply(., '[[', 2) %>% gsub("[[:space:]]", "", .) %>% unlist() %>% as.numeric()
breakpoint = append(breakpoint, rep(NA, 16), after=48)
extra_presses = grep("0:", read_pr, value = T) %>% strsplit("       ") %>% sapply(., '[', 4)%>% gsub("[[:space:]]", "", .) %>% unlist() %>% as.numeric()
extra_presses = append(extra_presses, rep(NA, 16), after=48)

pr_df <- data.frame(subject = gsub(".*Subject: ", "", grep("Subject", read_pr, value = T)),
                    cohort = str_match(grep("Subject", read_pr, value = T), "C\\d{2}") %>% unlist() %>% as.character(),  
                    exp = toupper(sub('.*HS', '', grep("Subject", read_pr, value = T))) %>% gsub(":SUBJECT.*", "", .),
                    start_date = gsub(".*Start Date: ", "", grep("Start Date:", read_pr, value = T)),
                    start_time = gsub(".*Start Time: ", "", grep("Start Time", read_pr, value = T)),
                    left_inactivepresses = gsub(".*A: ", "", grep("A:", read_pr, value = T)) %>% gsub("[[:space:]]", "", .) %>% as.numeric(),
                    right_activepresses = gsub(".*G: ", "", grep("G:", read_pr, value = T)) %>% gsub("[[:space:]]", "", .) %>% as.numeric(),
                    breakpoint = breakpoint,
                    extra_presses = extra_presses,
                    filename = sub(".*/.*/.*/", '', grep("Subject", read_pr, value = T)) %>% gsub(":Subject.*", "", .),
                    directory = str_match(grep("Subject", read_pr, value = T) %>% gsub("-Subject.*", "", .), "New_medassociates|Old") %>% unlist() %>% as.character()
) %>% arrange(cohort, exp) %>% 
  mutate(start_date = lubridate::mdy(format(as.Date(start_date, "%m/%d/%y"), "%m/%d/20%y")),
         start_time = chron::chron(times = start_time)) %>% 
  mutate_if(is.factor, as.character)

#looking for missing 0: arrays 
# test <- system("grep -iEa1r --no-group-separator  \"(Subject|H):\" */New_medassociates/PR* | grep -iE \"(Subject| 0):\"", intern = T) # THE REASON WHY THIS DOESN'T WORK ON OLD IS BC THE FORMAT IS DIFFERENT
# C02HSPR03 trouble file

