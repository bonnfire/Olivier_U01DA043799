## extract raw 
setwd("~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS")

# after cohort 5, there are only new files

## USEFUL FUNCTIONS

# FOR ~NEW~ DIRECTORIES
# #extract names to be assigned for various tables later
process_subjects_new <- function(x){
  
  read_subjects_new <- function(x){
    subjects <- fread(paste0("awk '/Subject/{print NR \"_\" $2}' ", "'", x, "'"),fill = T,header=F)
    subjects$filename <- x
    return(subjects)
  }

  read_meta_subjects_new <- function(x){
    date_time <- fread(paste0("awk '/Start/{print $3}' ", "'", x, "'", " | sed 'N;s/\\r\\n/_/g'"),fill = T,header=F)
    return(date_time)
  }
  date_time <- lapply(x, read_meta_subjects_new) %>% rbindlist() 
  
  names_sha_append <- lapply(x, read_subjects_new) %>% rbindlist() %>% rename("labanimalid"="V1") %>%
    cbind(., date_time) %>% 
    mutate(labanimalid = paste0( str_extract(labanimalid, "\\d+"), "_", 
                                str_extract(toupper(labanimalid), "[MF]\\d{1,3}"), "_",
                                str_extract(filename, "C\\d+"), "_",
                                sub('.*HS', '', toupper(filename)), "_",
                                sub(".*/.*/.*/", '', filename), "_",
                                V1)) %>% # subject id, cohort, experiment, file/location perhaps
  select(-V1)
  
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
    subject_old <- fread(paste0("grep -inEA1 --no-group-separator \"ratnumber|boxnumber\" ", "'", x, "'", "| grep -vE \"Rat|Box\""), header = F)
    subject_old$filename <- x 
    return(subject_old)
  }
  
  subjects_old <- lapply(x, read_subjects_old) %>% rbindlist() %>% 
    separate(V1, into = c("row", "value"), sep = "-", remove = T) %>% 
    mutate(box_id = ifelse((row_number() %% 2) == 0, "labanimalid", "box"))
  
  box <- subjects_old %>% dplyr::filter(box_id == "box") %>% select(value) %>% unlist() %>% as.character()
  labanimalid <- subjects_old %>% dplyr::filter(box_id == "labanimalid") %>% select(value) %>% unlist() %>% as.character()
  filename <- subjects_old %>% dplyr::filter(box_id == "box") %>% select(filename) %>% unlist() %>% as.character()
  row <- subjects_old %>% dplyr::filter(box_id == "box") %>% select(row) %>% unlist() %>% as.numeric()
  
  box_id_bind <- data.frame(box = box, 
                            labanimalid = labanimalid, 
                            filename = filename,
                            row = row) %>% 
    # mutate_all(as.character) %>% 
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
    select(one_of("labanimalid", "filename", "row"))
  # 
  # %>% dplyr::filter(!grepl("^NA", labanimalid))
  
  return(box_id_bind)
}


# FOR ~OLD~ DIRECTORIES
read_fread_old <- function(x, varname){
  
  fread_old_statements <- data.frame(varname = c("leftresponses", "rightresponses", "rewards"),
                                 statement = c("awk '/^BinsInActiveResponses/{flag=1;next}/endl/{flag=0}flag' ",
                                               "awk '/^ResponsesActBins/{flag=1;next}/endl/{flag=0}flag' ",
                                               "awk '/totalRewards/{flag=1;next}/TotalResponses/{flag=0; print NR}flag' ")) 
                                               # "awk '/BinRewards/{flag=1;next}/endl/{flag=0}flag' "))  #### 	 In=L Act=R  Rew=W InTS=U	ActTS=Y  RewTS=V  RewIRI=Z 	
  statement <- fread_old_statements[which(fread_old_statements$varname == varname),]$statement
  rawdata <- fread(paste0(statement, "'", x, "'"), fill = T, header = F)
  rawdata$filename <- x
  
  return(rawdata)
}

### XX 
## repeating the awk statements with conditions in shell script and trying to fread shell script

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

sha_subjects_new <- process_subjects_new(sha_new_files) %>% separate(labanimalid, c("row", "labanimalid"), sep = "_", extra = "merge") %>% 
  arrange(filename, row)
read_rewards_new <- function(x){
  rewards <- fread(paste0("awk '/W:/{flag=1;next}/5:/{flag=0}flag' ", "'", x, "' | awk '/0:/{print $2}'"))
  return(rewards)
}
sha_rewards_new <- lapply(sha_new_files, read_rewards_new) %>% rbindlist() %>% bind_cols(sha_subjects_new) %>% 
  separate(labanimalid, into = c("", "labanimalid", "cohort", "exp", "filename", "date", "time"), sep = "_") %>% 
  mutate(date = lubridate::mdy(date), time = chron::chron(times = time)) %>%  
  left_join(., date_time_subject_df_comp %>% 
              select(cohort, exp, filename, valid, start_date, start_time) %>% 
              rename("date" = "start_date", "time" = "start_time"), 
            by = c("cohort", "exp", "filename", "date", "time")) %>% 
  dplyr::filter(valid == "yes") %>% 
  rename("rewards" = "V1") %>% 
  mutate(time = as.character(time))

###### OLD FILES ##############

sha_subjects_old <- process_subjects_old(sha_old_files)
sha_rewards_old <- lapply(sha_old_files, read_fread_old, "rewards") %>% rbindlist() 

rewards <- sha_rewards_old %>% dplyr::filter(row_number() %% 2 == 1) %>% select(V1) %>% unlist() %>% as.character()
row <- sha_rewards_old %>% dplyr::filter(row_number() %% 2 == 0) %>% select(V1) %>% unlist() %>% as.character()
filename <- sha_rewards_old %>% dplyr::filter(row_number() %% 2 == 1) %>% select(filename) %>% unlist() %>% as.character()

rewards_bind <- data.frame(rewards = rewards, 
                           row = row, 
                           filename = filename) 

sha_rewards_old <- rewards_bind %>% arrange(filename, as.numeric(as.character(row))) %>% bind_cols(sha_subjects_old %>% arrange(filename, row)) %>% 
  select(-c("filename1", "row", "row1")) %>% 
  separate(labanimalid, into = c("labanimalid", "box", "cohort", "exp", "computer", "date", "valid"), sep = "_") %>% 
  mutate(date = lubridate::ymd(date),
         rewards = rewards %>% unlist() %>% as.character() %>% as.numeric(),
         filename = filename %>% unlist() %>% as.character()) %>%  
  dplyr::filter(valid == "valid") 




















################################
########## LGA #################
################################

###### NEW FILES ##############
lga_new_files <- grep(grep(list.files(path = ".", recursive = T, full.names = T), pattern = ".*txt", inv = T, value = T), pattern = ".*LGA", value = T) # 328 files
lga_subjects_new <- process_subjects_new(lga_new_files) %>% separate(labanimalid, c("row", "labanimalid"), sep = "_", extra = "merge") %>% 
  arrange(filename, as.numeric(row)) %>% select(-row)
read_rewards_new <- function(x){
  rewards <- fread(paste0("awk '/W:/{flag=1;next}/5:/{flag=0}flag' ", "'", x, "' | awk '/0:/{print NR \"_\" $2}'"), header = F, fill = T)
  rewards$filename <- x
  return(rewards)
}
lga_rewards_new <- lapply(lga_new_files, read_rewards_new) %>% rbindlist() %>% separate(V1, into = c("row", "rewards"), sep = "_") %>% arrange(filename, as.numeric(row)) %>% select(-row) %>% 
  bind_cols(lga_subjects_new) %>% 
  separate(labanimalid, into = c("labanimalid", "cohort", "exp", "filename", "date", "time"), sep = "_") %>% 
  mutate(date = lubridate::mdy(date), time = chron::chron(times = time), rewards = as.numeric(rewards)) %>%  
  left_join(., date_time_subject_df_comp %>% 
              select(cohort, exp, filename, valid, start_date, start_time) %>% 
              rename("date" = "start_date", "time" = "start_time"), 
            by = c("cohort", "exp", "filename", "date", "time")) %>% 
  dplyr::filter(valid == "yes") %>% 
  mutate(time = as.character(time))

###### OLD FILES ##############
lga_old_files <- grep(list.files(path = ".", recursive = T, full.names = T), pattern = ".*Old.*LGA", value = T) # 424 files
lga_subjects_old <- process_subjects_old(lga_old_files)
lga_rewards_old <- lapply(lga_old_files, read_fread_old, "rewards") %>% rbindlist() 

rewards <- lga_rewards_old %>% dplyr::filter(row_number() %% 2 == 1) %>% select(V1) %>% unlist() %>% as.character()
row <- lga_rewards_old %>% dplyr::filter(row_number() %% 2 == 0) %>% select(V1) %>% unlist() %>% as.character()
filename <- lga_rewards_old %>% dplyr::filter(row_number() %% 2 == 1) %>% select(filename) %>% unlist() %>% as.character()

rewards_bind <- data.frame(rewards = rewards, 
                           row = row, 
                           filename = filename) 

lga_rewards_old <- rewards_bind %>% arrange(filename, as.numeric(as.character(row))) %>% bind_cols(lga_subjects_old %>% arrange(filename, row)) %>% 
  select(-c("filename1", "row", "row1")) %>% 
  separate(labanimalid, into = c("labanimalid", "box", "cohort", "exp", "computer", "date", "valid"), sep = "_") %>% 
  mutate(date = lubridate::ymd(date),
         rewards = rewards %>% unlist() %>% as.character() %>% as.numeric(),
         filename = filename %>% unlist() %>% as.character()) %>%  
  dplyr::filter(valid == "valid") 


## exclude the following files - Olivier (1/24 meeting)
### -----
### -----










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

