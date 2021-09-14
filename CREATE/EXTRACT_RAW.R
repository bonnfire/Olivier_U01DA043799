## extract raw 
setwd("~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS")


# after cohort 5, there are only new files

## USEFUL FUNCTIONS

# FOR ~NEW~ DIRECTORIES
# #extract names to be assigned for various tables later
process_subjects_new <- function(x){
  
  read_subjects_new <- function(x) {
    subjects <-
      fread(
        paste0("awk '/Subject/{print NR \"_\" $2}' ", "'", x, "'"),
        fill = T,
        header = F
      )
    subjects$filename <- x
    return(subjects)
  }
  
  read_meta_subjects_new <- function(x) {
    date_time <-
      fread(
        paste0("awk '/Start/{print $3}' ", "'", x, "'", " | sed 'N;s/\\r\\n/_/g'"), # requires this line bc it has two fields of information
        fill = T,
        header = F
      )
    return(date_time)
  }
  date_time <- lapply(x, read_meta_subjects_new) %>% rbindlist()
  
  read_meta_box_new <- function(x) {
    box_new <-
      fread(
        paste0("awk '/Box/{print $2}' ", "'", x, "'"),
        fill = T,
        header = F
      )
    return(box_new)
  }
  box_new <- lapply(x, read_meta_box_new) %>% rbindlist()
  
  
  names_sha_append <- lapply(x, read_subjects_new) %>% rbindlist() %>% rename("labanimalid"="V1") %>%
    cbind(., date_time) %>% 
    rename("date_time" = "V1") %>% 
    cbind(., box_new) %>% 
    rename("box" = "V1") %>%
    mutate(labanimalid = paste0( str_extract(labanimalid, "\\d+"), "_",
                                 str_extract(toupper(labanimalid), "[MF]\\d{1,3}"), "_",
                                 str_extract(filename, "C\\d+"), "_",
                                 sub('.*HS', '', toupper(filename)), "_",
                                 sub(".*/.*/.*/", '', filename), "_",
                                 date_time, "_",
                                 box)) %>% # subject id, cohort, experiment, file/location perhaps
    select(-c("date_time", "box"))
  
  return(names_sha_append)
  
}


# bc some animals don't have id's under subject, mistake in manual entry; extract subject information from experiment or group? 

# #notice system call, doing it recurisvely for the entire directory; for these cases, the group number is not 0 and the subject is found elsewhere.  # 12/18 only pr
# only calling cohort 8?
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
      
      processeddata_df <- processeddata_df %>%
        mutate(intertrial_time = lead(timestamps) - timestamps,
               bin = as.character(bin))
      
      return(processeddata_df)
    }) 
  }
  
  
  return(processeddata)
}

## 08/12/2020
lapply(sha_new_files[1:2], read_fread_new, "rewardstimestamps")

## 08/02/2020
# sha_new_files <- grep(list.files(path = ".", recursive = T, full.names = T), pattern = ".*C01.*New.*SHA", value = T) 
lapply( grep(list.files(path = ".", recursive = T, full.names = T), pattern = ".*C01.*New.*SHA", value = T),  read_fread_new, "leftresponses")



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
    mutate(date = str_extract(filename, "\\d{8}") %>% lubridate::ymd(),
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
                                                   "awk '/totalRewards/{flag=1;next}/TotalResponses/{flag=0}flag' ")) 
  # "awk '/BinRewards/{flag=1;next}/endl/{flag=0}flag' "))  #### 	 In=L Act=R  Rew=W InTS=U	ActTS=Y  RewTS=V  RewIRI=Z 	
  statement <- fread_old_statements[which(fread_old_statements$varname == varname),]$statement
  rawdata <- fread(paste0(statement, "'", x, "' | nl -s _ | sed \"s/[[:blank:]]//g\""), fill = T, header = F)
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
  
  x <- x %>% 
    mutate(iritime = as.numeric(iritime),
           iritime = iritime/100)
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
  return(x)
}



## 08/12/2020
lapply(sha_old_files[1], read_iri_old) %>% lapply(convert_iri_matrix_to_df)

# FOR ~NEW~ DIRECTORIES
# NON EXPERIMENT SPECIFIC 


## TO VALIDATE ENTRIES  
## extract date and start time/end time to determine valid sessions
read_date_time_subject <- system("grep -a7r --no-group-separator \"Start Date: \" . | grep -E \"(Start Date|End Date|Subject|Box|Start Time|End Time):\"", intern = T)
read_date_time_subject <- gsub("\\r", "", read_date_time_subject)
read_date_time_subject <- read_date_time_subject[!grepl("/LGA/:", read_date_time_subject)] # remove the duplicate file # removed on 1/27

date_time_subject <- data.frame(labanimalid = gsub(".*Subject: ", "", grep("Subject", read_date_time_subject, value = T)) %>% toupper,
                                cohort = str_match(grep("Subject", read_date_time_subject, value = T), "C\\d{2}") %>% unlist() %>% as.character(),  
                                exp = toupper(sub('.*HS', '', grep("Subject", read_date_time_subject, value = T) %>% gsub("-Subject.*", "", .))),
                                start_date = gsub(".*Start Date: ", "", grep("Start Date:", read_date_time_subject, value = T)),
                                end_date = gsub(".*End Date: ", "", grep("End Date:", read_date_time_subject, value = T)),
                                box = gsub(".*Box: ", "", grep("Box", read_date_time_subject, value = T)),
                                start_time = gsub(".*Start Time: ", "", grep("Start Time", read_date_time_subject, value = T)),
                                end_time = gsub(".*End Time: ", "", grep("End Time", read_date_time_subject, value = T)),
                                filename = sub(".*/.*/.*/", '', grep("Subject", read_date_time_subject, value = T)) %>% gsub("-Subject.*", "", .),
                                directory = str_match(grep("Subject", read_date_time_subject, value = T) %>% gsub("-Subject.*", "", .), "New_medassociates|Old") %>% unlist() %>% as.character()
) # 13548 

date_time_subject_mut <- date_time_subject %>% 
  mutate_if(is.factor, as.character) %>% 
  mutate(start = lubridate::mdy_hms(paste0(format(as.Date(start_date, "%m/%d/%y"), "%m/%d/20%y"), start_time)),
         end = lubridate::mdy_hms(paste0(format(as.Date(end_date, "%m/%d/%y"), "%m/%d/20%y"), end_time))) %>% 
  mutate(exp_dur_min = difftime(end, start, units= "mins") %>% as.numeric) %>% 
  select(-matches("date|time")) %>% 
  mutate(labanimalid = replace(labanimalid, labanimalid == "M7678", "M768"),
         labanimalid = replace(labanimalid, labanimalid == "X", "F507"),
         labanimalid = replace(labanimalid, labanimalid == "717", "F717")) %>% # verified by box
  mutate(labanimalid = if_else(grepl("^C0", labanimalid), str_match(labanimalid,"[FM]\\d{1,3}" %>% unlist() %>% as.character()), labanimalid %>% as.character()))

### 
## 3/4/2020 plotting the durations 
date_time_subject_mut %>% 
  mutate(exp_abv = sub("^([[:alpha:]]*).*", "\\1", exp)) %>% 
  ggplot(aes(x = exp_dur_min, color = exp_abv)) + geom_density() + facet_grid( ~ cohort)


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
    arrange(labanimalid, start) 
  return(x)
}) %>% rbindlist(., idcol = "cohort")

# replace rbindlist... with openxlsx::write.xlsx(., "labanimalid_assign_bybox.xlsx") to create the excel sheets that I sent to their lab 
# subject0[[3]] <- NULL
# subject0 %>% openxlsx::write.xlsx(., "labanimalid_assign_bybox.xlsx")

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
  arrange(cohort, start, as.numeric(box)) %>% 
  distinct() %>%  # needed bc otherwise subject0 will "double count" the "reference" rows
  mutate(exp = gsub("-.*", "", exp),
         exp = replace(exp, as.numeric(str_extract(cohort, "\\d+")) > 6&grepl("SHOCK", exp), "SHOCK03"),
         room = ifelse(grepl("[[:alnum:]]+C\\d{2}HS", filename), gsub("C\\d{2}HS.*", "", filename), NA)) # for all cohorts later than cohort 6, they only use one shock value, but we can change it to shock03 so that we can have uniformity 
## XX Note: remove all but SHOCK03 and Pre Shock - Olivier (from meeting 1/24)

# fix these subjects because their labanimalids should not be this long
date_time_subject_df %>% subset(!cohort %in% c("C10", "C11")&nchar(labanimalid)>4)
date_time_subject_df <- date_time_subject_df %>% 
  mutate(labanimalid = replace(labanimalid, labanimalid == "F8256", "F826"),
         labanimalid = replace(labanimalid, labanimalid == "M5556", "M556")) %>% 
  mutate(exp = replace(exp, grepl("PRESHOCK", filename), "PRESHOCK"))

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
  # mutate(date_pr19 = replace(date_pr19, cohort == "cohort5", lubridate::ymd("2018-09-19")),
  mutate(date_sha02 = replace(date_sha02, cohort == "cohort5", lubridate::ymd("2018-07-31"))) %>% select(matches("date|cohort")) %>% distinct() %>% # for record keeping, make sure to make this change on the actual excel! 
  gather(v, value, date_sha01:date_lga23) %>% 
  separate(v, c("date", "exp")) %>% 
  arrange(cohort) %>% 
  select(-date) %>% 
  mutate(cohort = paste0("C", str_pad(gsub("COHORT", "", toupper(cohort)),  2, "left","0")),
         exp = toupper(exp),
         value = lubridate::ymd(value)) %>% 
  rename("excel_date" = "value")

date_time_subject_df_comp <- left_join(date_time_subject_df, cohorts_exp_date, by = c("cohort", "exp")) %>%
  mutate(start_date = as.Date(start),
         start_time = format(start, "%H:%M:%S") %>% chron::chron(times = .)) %>% 
  mutate(valid = case_when(
    grepl("SHOCK", exp) & exp_dur_min > 58 & excel_date == start_date ~ "yes",
    grepl("SHA", exp) & exp_dur_min > 115 & excel_date == start_date~ "yes",
    grepl("LGA", exp) & exp_dur_min > 355 & excel_date == start_date~ "yes",
    grepl("PR", exp) & exp_dur_min > 60 & excel_date == start_date~ "yes"),
    valid = replace(valid, is.na(valid), "no")
  )   # 9808 for cohorts C01-C09 (no C06) ## change the minimum times - Olivier (from 1/24 meeting)
# mutate()# temporarily give a free pass for cohort 9


date_time_subject_df_comp <- date_time_subject_df_comp %>% 
  subset(!(labanimalid == "806"&box =="6"&room=="BSB273C"&exp =="LGA02")) # remove after manual check


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


## age table (at start of session), box, directory, rooom
subjects_exp_age_source <-  rat_info_allcohort_xl_df[, c("cohort", "labanimalid", "rfid")] %>% 
  subset(grepl("^[MF]", labanimalid)) %>%
  distinct() %>% 
  left_join(
    rbind(cohorts_exp_date, 
          olivierxl_df %>% 
            subset(cohort == "C09") %>% 
            select(cohort, matches("date")) %>% 
            gather("exp", "date", -cohort) %>% distinct() %>% 
            mutate(exp = gsub("date_", "", exp) %>% toupper) %>%
            rename("excel_date" = "date"))
    , by = c("cohort")
  ) %>% # XX TEMP get cohort 9 dates
  rename("start_date" = "excel_date") %>% 
  # left_join(date_time_subject_df_comp %>% 
  #             subset(valid == "yes") %>% 
  #             select(labanimalid, cohort, exp, box, room), 
  #           by = c("labanimalid", "cohort", "exp")) %>% ## use rat_info_allcohort_xl_df to get rfid and use the comp to get the start date for exp
  # left_join()
  # mutate(comment_date = "NA") %>%
  # mutate(comment_date = replace(comment, !is.na(start_date), ""))
  left_join(WFU_OlivierCocaine_test_df[, c("rfid", "dob")], by = c("rfid")) %>% # use WFU_OlivierCocaine_test_df to get dob
  left_join(WFU_OlivierOxycodone_naive_test[, c("rfid", "dob")], by = c("rfid")) %>% # extract the dob of the scrubs from subset(is.na(dob)) (XX temp fix)
  mutate_all(as.character) %>% 
  mutate(dob = coalesce(dob.x, dob.y), 
         source_dob = ifelse(!is.na(dob), "wfu", "NA")) %>% # merge the two dob columns
  left_join(rat_info_allcohort_xl_df[, c("labanimalid", "d_o_b")] %>% 
              mutate(d_o_b = ifelse(!grepl("-", d_o_b), openxlsx::convertToDate(d_o_b) %>% as.character, d_o_b)), by = "labanimalid") %>%
  mutate(source_dob = replace(source_dob, source_dob == "NA"&!is.na(d_o_b), "olivier excel")) %>%
  mutate(dob = coalesce(dob, d_o_b)) %>% 
  mutate_at(vars(matches("^(start_date|dob)")), as.Date) %>% 
  mutate(age = difftime(start_date, dob, units = c("days")) %>% as.numeric) %>% # calculate the age
  distinct(cohort, labanimalid, rfid, exp, age, source_dob, source_dob) %>% 
  subset(grepl("SHA|SHOCK|LGA(0[1-9]|1[1-4])|PR", exp)) %>%
  mutate(exp = paste0(gsub("(\\D+)(\\d+)", "\\1_\\2", tolower(exp)), "_age")) %>% 
  mutate(source_date = ifelse(!is.na(age), "olivierexcel", "NA")) %>% 
  group_by(labanimalid) %>% fill(age) %>% 
  mutate(source_date = replace(source_date, source_date == "NA", "fill by previous date")) %>% 
  ungroup() 

subjects_exp_age <- subjects_exp_age_source %>%
  select(-matches("source")) %>% 
  spread(exp, age) # spread to get _age columns




# mutate_at(vars(-matches("labanimalid|cohort")), list(esc = ~.-dob)) # calculate the age 


## box table
# record of the raw data (combined with Brent's knowledge for the rooms)
box_metadata_long <- rat_info_allcohort_xl_df[, c("cohort", "labanimalid", "rfid")] %>% 
  subset(grepl("^[MF]", labanimalid)&!is.na(labanimalid)) %>%
  distinct() %>% 
  left_join(cohorts_exp_date %>% select(-excel_date), by = c("cohort")) %>% # get all exps 
  left_join(date_time_subject_df_comp %>% 
              subset(valid == "yes") %>% 
              select(labanimalid, cohort, exp, box, room), 
            by = c("labanimalid", "cohort", "exp")) %>% ## use rat_info_allcohort_xl_df to get rfid and use the comp to get the box and room 
  full_join(rewards[, c("directory", "cohort", "labanimalid", "exp", "room", "box", "computer")], by = c("cohort", "labanimalid", "exp")) %>% # get the old directory animals since the comp object only gives new directory animals
  mutate_all(as.character) %>% 
  mutate(box = coalesce(box.x, box.y),
         room = coalesce(room.y, room.x), # chose this order bc room y has less missingness
         room_computer = coalesce(room, computer)) %>% 
  select(-matches("\\.[xy]$|^computer|^room$")) %>% 
  rowwise() %>% 
  mutate(room_computer = replace(room_computer, is.na(room_computer)&grepl("F\\d?([1-9]|0[1-9]|1[0-6]|90)$", labanimalid)&parse_number(cohort)<7, "MED1110"),
         room_computer = replace(room_computer, is.na(room_computer)&grepl("F\\d?(1[7-9]|2[0-2])$", labanimalid)&parse_number(cohort)<7, "K1"),
         room_computer = replace(room_computer, is.na(room_computer)&grepl("F\\d?(2[3-8]$)", labanimalid)&parse_number(cohort)<7, "K2"),
         room_computer = replace(room_computer, is.na(room_computer)&grepl("(F\\d?(29|30))|(M\\d?(5[1-4]))$", labanimalid)&parse_number(cohort)<7, "K3"),
         room_computer = replace(room_computer, is.na(room_computer)&grepl("M\\d?(5[5-9]|6[0-b8])$", labanimalid)&parse_number(cohort)<7, "MED1113"),
         room_computer = replace(room_computer, is.na(room_computer)&grepl("M\\d?(69|7[0-4])$", labanimalid)&parse_number(cohort)<7, "Q2"),
         room_computer = replace(room_computer, is.na(room_computer)&grepl("(M\\d?(7[5-9]|80))$", labanimalid)&parse_number(cohort)<7, "Q3"),
         room_computer = replace(room_computer, grepl("SHOCK", exp)&parse_number(cohort)>5, "BSB273B"),
         room_computer = replace(room_computer, grepl("SHOCK", exp)&parse_number(cohort)<=5, "MED1110")) %>% #using Brent's comments 09/23/2020 
  ungroup() %>% 
  subset(!is.na(labanimalid)|grepl("[MF]", labanimalid))


## qc 
box_metadata_long <- box_metadata_long %>% 
  subset(!grepl("[MF]0\\d", labanimalid)) %>% # remove cases like F03
  subset(!grepl("[MF]\\d+-\\d", labanimalid)) # remove cases like M156-1 and F79-1

## create dataframe with all filled for NA 
# record of the raw data (combined with Brent's knowledge for the rooms) ON TOP OF the knowledge theoretically, each animal should be run on the same boxes throughout the exps
box_metadata_long_fill_2 <- box_metadata_long %>% # for those animals that don't have raw data # exlude animals for all na
  naniar::replace_with_na_all(~.x %in% c("NA", "<NA>", "N/A")) %>% 
  mutate(phenotype_box = ifelse(!is.na(box), paste0(room_computer, "-", box), NA)) %>% 
  subset(!grepl("SHOCK", exp)) %>% 
  group_by(labanimalid) %>%
  dplyr::filter(!all(is.na(box))) %>% 
  arrange(desc(!is.na(phenotype_box))) %>% 
  # fill(phenotype_box) %>%
  mutate(phenotype_box = zoo::na.locf(phenotype_box)) %>%
  ungroup() %>% 
  rbind(box_metadata_long %>% # for those animals that don't have raw data
          naniar::replace_with_na_all(~.x %in% c("NA", "<NA>", "N/A")) %>% 
          mutate(phenotype_box = ifelse(!is.na(box), paste0(room_computer, "-", box), NA)) %>% 
          subset(grepl("SHOCK", exp)) %>% 
          group_by(labanimalid) %>%
          dplyr::filter(!all(is.na(box))) %>% 
          arrange(desc(!is.na(phenotype_box))) %>% 
          # fill(phenotype_box) %>%
          mutate(phenotype_box = zoo::na.locf(phenotype_box)) %>%
          ungroup()) %>% 
  rbind(date_time_subject_df_comp %>% 
          subset(cohort == "C09") %>%
          mutate(labanimalid = replace(labanimalid, labanimalid == "F951", "M951"),
                 labanimalid = replace(labanimalid, labanimalid == "F952", "M952")) %>% 
          mutate(phenotype_box = paste0(room, "-", box)) %>% 
          rename("room_computer" = "room") %>% 
          mutate(directory = "new_sha") %>% 
          subset(!(labanimalid == "F925"&phenotype_box == "BSB273E-9")) %>% # animal died before SHA06, should not have LGA sessions
          left_join(rat_info_allcohort_xl_df[, c("labanimalid", "rfid")], by = "labanimalid") %>% # give df rfid column
          distinct(cohort, labanimalid, rfid, exp, directory, box, room_computer, phenotype_box)
  ) ## XX TEMP FOR COHORT 9


# box_metadata_long_fill <- box_metadata_long %>% # for those animals that don't have raw data
#   naniar::replace_with_na_all(~.x %in% c("NA", "<NA>", "N/A")) %>% 
#   mutate(phenotype_box = ifelse(!is.na(box), paste0(room_computer, "-", box), NA)) %>% 
#   subset(!grepl("SHOCK", exp)) %>% 
#   group_by(labanimalid) %>%
#   fill(phenotype_box) %>%
#   # mutate(phenotype_box = zoo::na.locf(phenotype_box, fromLast = T)) %>%
#   ungroup() %>% 
#   rbind(box_metadata_long %>% # for those animals that don't have raw data
#           naniar::replace_with_na_all(~.x %in% c("NA", "<NA>", "N/A")) %>% 
#           mutate(phenotype_box = ifelse(!is.na(box), paste0(room_computer, "-", box), NA)) %>% 
#           subset(grepl("SHOCK", exp)) %>% 
#           group_by(labanimalid) %>%
#           fill(phenotype_box) %>%
#           # mutate(phenotype_box = zoo::na.locf(phenotype_box, fromLast = T)) %>%
#           ungroup())

# by cohort and sex for all non shock exps
box_metadata_long_fill <- box_metadata_long_fill %>% mutate(phenotype_box = replace(phenotype_box, cohort == "C02"&grepl("F(80|10[6-9])", labanimalid)&!grepl("SHOCK", exp), "No phenotype data"))
box_metadata_long_fill <- box_metadata_long_fill %>% mutate(phenotype_box = replace(phenotype_box, cohort == "C02"&grepl("F90", labanimalid)&!grepl("SHOCK", exp), "MED1110-11"))
box_metadata_long_fill %>% subset(cohort == "C02"&grepl("F", labanimalid)&!grepl("SHOCK", exp)) %>% group_by(labanimalid) %>% mutate(phenotype_box = zoo::na.locf(phenotype_box, fromLast=T)) %>% ungroup()%>% distinct(labanimalid, phenotype_box) %>% 
  add_count(labanimalid) %>% subset(n!= 1) %>% View()
# box_metadata_long_fill %>% subset(cohort == "C02"&grepl("F80", labanimalid)&!grepl("SHOCK", exp))%>% select(labanimalid, exp, phenotype_box) %>% spread("exp", "phenotype_box")
# box_metadata_long_fill %>% subset(cohort == "C02"&grepl("F9[36]", labanimalid)&!grepl("SHOCK", exp))%>% subset(labanimalid == "F93"&phenotype_box != "MED1110-14"| labanimalid == "F96"&phenotype_box != "MED1110-15")
box_metadata_long_fill %>% subset(cohort == "C02"&grepl("M", labanimalid)&!grepl("SHOCK", exp)) %>% group_by(labanimalid) %>% mutate(phenotype_box = zoo::na.locf(phenotype_box, fromLast=T)) %>% ungroup()%>% distinct(labanimalid, phenotype_box) %>% add_count(labanimalid) %>% subset(n!= 1) %>% View()
# box_metadata_long_fill %>% subset(cohort == "C02"&grepl("M17[35]", labanimalid)&!grepl("SHOCK", exp))%>% subset(labanimalid == "M173"&phenotype_box != "Q3-1"| labanimalid == "M175"&phenotype_box != "Q3-3")
box_metadata_long_fill %>% subset(cohort == "C03"&grepl("M", labanimalid)&!grepl("SHOCK", exp))%>% distinct(labanimalid, phenotype_box) %>% add_count(labanimalid) %>% subset(n== 1&is.na(phenotype_box)|n!=1) %>% select(labanimalid) %>% mutate(labanimalid = gsub("M", "", labanimalid)) %>% unlist() %>% paste0(collapse = ", ")
box_metadata_long_fill %>% subset(cohort == "C04"&grepl("M", labanimalid)&!grepl("SHOCK", exp))%>% distinct(labanimalid, phenotype_box) %>% add_count(labanimalid) %>% subset(n==2&phenotype_box %>% is.na) %>% select(labanimalid) %>% mutate(labanimalid = gsub("M", "", labanimalid)) %>% unlist() %>% paste0(collapse = ", ")

# by cohort and sex for all shock exps
# once you filter out the n==1&is.na(phenotype_box), use group_by(labanimalid) %>% mutate(phenotype_box = zoo::na.locf(phenotype_box, fromLast=T)) %>% ungroup()
box_metadata_long_fill %>% subset(cohort == "C01"&grepl("F", labanimalid)&grepl("SHOCK", exp)) %>% distinct(labanimalid, phenotype_box) %>% add_count(labanimalid) %>% subset(n== 1&is.na(phenotype_box)|n!=1) %>% View()


# clean up the duplicates for non shock
box_metadata_long_fill %>% 
  subset(!grepl("SHOCK", exp)) %>% 
  distinct(labanimalid, room_box) %>% 
  add_count(labanimalid) %>% 
  subset(n!=1)

# make sure that shock has value


# breakdown by cohorts 
# boxes for cohort 1
setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE")
c01_boxes <- rewards %>% # already subsetted to "valid" or "yes" in valid
  subset(cohort == "C01") %>% 
  select(cohort, labanimalid, exp, room_box) %>% 
  subset(exp %in% c("SHA08", "SHA09", "SHA10", "SHOCK03", "PR01", "PR02", "PR03", "LGA11", "LGA12", "LGA13", "LGA14")) %>% 
  spread(key = "exp", value = "room_box") %>% 
  mutate(labanimalid_num = parse_number(labanimalid),
         sex = str_extract(labanimalid, "[MF]")) %>% 
  arrange(cohort, sex, labanimalid_num) %>% select(-c("labanimalid_num", "sex")) %>% 
  
  
  openxlsx::write.xlsx(c01_boxes, file = "cocaine_c01_boxes_toqc.xlsx")

setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE")
c02_boxes <- rewards %>% # already subsetted to "valid" or "yes" in valid
  subset(cohort == "C02") %>% 
  mutate(room = ifelse(grepl("[[:alnum:]]+C\\d{2}HS", filename), gsub("([.]/.*/.*/.*/)?(.*)C\\d{2}HS.*", "\\2", filename), NA)) %>% 
  select(cohort, labanimalid, exp, box, room) %>% 
  naniar::replace_with_na_all( ~ .x %in% c("<NA>")) %>% 
  rowwise() %>% 
  mutate(room_box = ifelse(!is.na(room), paste0(room, "-", box), box)) %>%
  select(-c("box", "room")) %>% 
  subset(exp %in% c("SHA08", "SHA09", "SHA10", "SHOCK03", "PR01", "PR02", "PR03", "LGA11", "LGA12", "LGA13", "LGA14")) %>% 
  ungroup() %>% 
  spread(key = "exp", value = "room_box") %>% 
  mutate(labanimalid_num = parse_number(labanimalid),
         sex = str_extract(labanimalid, "[MF]")) %>% 
  arrange(cohort, sex, labanimalid_num) %>% select(-c("labanimalid_num", "sex"))
openxlsx::write.xlsx(c02_boxes, file = "cocaine_c02_boxes_toqc.xlsx")


c03_boxes <- rewards %>% # already subsetted to "valid" or "yes" in valid
  subset(cohort == "C03") %>% 
  mutate(room = ifelse(grepl("[[:alnum:]]+C\\d{2}HS", filename), gsub("([.]/.*/.*/.*/)?(.*)C\\d{2}HS.*", "\\2", filename), NA)) %>% 
  select(cohort, labanimalid, exp, box, room) %>% 
  naniar::replace_with_na_all( ~ .x %in% c("<NA>")) %>% 
  rowwise() %>% 
  mutate(room_box = ifelse(!is.na(room), paste0(room, "-", box), box)) %>%
  select(-c("box", "room")) %>% 
  subset(exp %in% c("SHA08", "SHA09", "SHA10", "SHOCK03", "PR01", "PR02", "PR03", "LGA11", "LGA12", "LGA13", "LGA14")) %>% 
  ungroup() %>% 
  spread(key = "exp", value = "room_box") %>% 
  mutate(labanimalid_num = parse_number(labanimalid),
         sex = str_extract(labanimalid, "[MF]")) %>% 
  arrange(cohort, sex, labanimalid_num) %>% select(-c("labanimalid_num", "sex"))
openxlsx::write.xlsx(c03_boxes, file = "cocaine_c03_boxes_toqc.xlsx")




# XX 09/17/2020 fix this eventually (maybe go back to the comp object and fix from there) box_metadata_long %>% distinct() %>% get_dupes(labanimalid, exp)
# box_metadata_wide <- box_metadata_long %>% 
#   subset(grepl("SHA08|SHOCK03|LGA11|PR\\d", exp)) %>%
#   mutate(exp = paste0(gsub("(\\D+)(\\d+)", "\\1_\\2", tolower(exp)), "_box")) %>% 
#   naniar::replace_with_na_all(~.x %in% c("<NA>")) %>% 
#   mutate(room = replace(room, is.na(room), "MED1110")) %>% 
#   mutate(box = ifelse(!is.na(room)&!is.na(box), paste0(box, "-", room), box)) %>% 
#   select(-room) %>%
#   distinct() %>% 
#   mutate(box = replace(box, is.na(box), "NA")) %>% 
#   spread(exp, box)



box_metadata_wide_2 <- box_metadata_long_fill_2 %>% 
  subset(grepl("SHA08|SHOCK03|LGA11|PR\\d", exp)) %>%
  distinct(cohort, labanimalid, rfid, exp, phenotype_box) %>% 
  mutate(exp = paste0(gsub("(\\D+)(\\d+)", "\\1_\\2", tolower(exp)), "_box")) %>%
  spread(exp, phenotype_box)


























################################  
########## SHA #################
################################

setwd("~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS")

###### NEW FILES ##############
sha_new_files_01_11 <- grep(grep(list.files(path = "~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS", recursive = T, full.names = T), pattern = ".*txt", inv = T, value = T), pattern = ".*SHA", value = T) %>% grep("C(0[1-9]|1[01])", ., value = T) # 254 files

# check cohorts and experiments represented 
sha_new_files_01_11 %>% as.data.frame() %>% 
  mutate(cohort = str_extract(., "C\\d+"),
         exp = str_extract(., "SHA\\d+")) %>% 
  select(cohort, exp) %>% table(exclude = NULL)

# check which ones are NA
sha_new_files_01_11 %>% as.data.frame() %>% 
  mutate(cohort = str_extract(., "C\\d+"),
         exp = str_extract(., "SHA\\d+")) %>% 
  subset(is.na(exp))

# extract subject, box, and rewards, active, inactive (rai)
sha_rai <- lapply(sha_new_files_01_11, function(x){
  setwd("~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS")
  internal_filename = fread(paste0("awk '/MSN:/{print $1 $2 $3}' '", x, "'"), header = F, fill = T)
  subject = fread(paste0("awk '/Subject:/{print $2}' '", x, "'"), header = F, fill = T)
  
  box = fread(paste0("awk '/Box:/{print $2}' '", x, "'"), header = F, fill = T)
  
  date = fread(paste0("awk '/^Start Date:/{print $3}' '", x, "'"), header = F, fill = T)
  enddate = fread(paste0("awk '/^End Date:/{print $3}' '", x, "'"), header = F, fill = T)
  
  time = fread(paste0("awk '/^Start Time:/{print $3}' '", x, "'"), header = F, fill = T)
  endtime = fread(paste0("awk '/^End Time:/{print $3}' '", x, "'"), header = F, fill = T)
  
  session_duration <- cbind(date = date, enddate = enddate) %>%
    cbind(time = time) %>% cbind(endtime = endtime) 
  names(session_duration) <- c("date", "enddate", "time", "endtime")
  session_duration <- session_duration %>% 
    mutate(start_datetime = lubridate::mdy_hms(paste(format(as.Date(date, "%m/%d/%y"), "%m/%d/20%y"), time)),
           end_datetime = lubridate::mdy_hms(paste(format(as.Date(enddate, "%m/%d/%y"), "%m/%d/20%y"), endtime)),
           session_duration = difftime(end_datetime, start_datetime, units = "mins") %>% as.numeric() %>% round(0))
  
  metadata = cbind(subject = subject, box = box) %>% 
    cbind(startdate = date) %>% 
    cbind(session_duration = session_duration$session_duration) %>% 
    cbind(internal_filename)
  
  names(metadata) <- c("subject", "box", "startdate", "session_duration", "internal_filename")
  
  if(grepl("GWAS 2hNEW", unique(internal_filename$V1))){
    rewards = fread(paste0("awk '/^B: /{print $2}' ", "'", x, "'"), header = F, fill = T)
    active = fread(paste0("awk '/^G: /{print $2}' ", "'", x, "'"), header = F, fill = T)
    inactive = fread(paste0("awk '/^A: /{print $2}' ", "'", x, "'"), header = F, fill = T)
  }
  else{
    rewards = fread(paste0("awk '/^W:/{flag=1;next}/5:/{flag=0}flag' ", "'",  x, "' | awk '/0:/{print $2}'"), header = F, fill = T)
    active =  fread(paste0("awk '/^R:/{flag=1;next}/5:/{flag=0}flag' ", "'",  x, "' | awk '/0:/{print $2}'"), header = F, fill = T)
    inactive = fread(paste0("awk '/^L:/{flag=1;next}/5:/{flag=0}flag' ", "'",  x, "' | awk '/0:/{print $2}'"), header = F, fill = T)
  }
  
  data = list("rewards" = rewards,
              "active" = active,
              "inactive" = inactive  
  ) %>% do.call(cbind, .)
  
  data <- cbind(metadata, data) %>% 
    rename_all(~ stringr::str_replace_all(., '[.](V1)?', '')) %>% 
    mutate(filename = x )
  
  return(data)
})

sha_rai_df <- sha_rai %>% 
  rbindlist(fill = T) %>% 
  mutate(cohort = str_match(filename, "C\\d+") %>% as.character,
         subject = str_extract(toupper(subject), "[MF]\\d+"),
         filename = gsub(".*/SHA/", "", filename),
         exp = gsub(".*HS", "", toupper(filename)) %>% gsub("-\\d+$", "", .),
         room = ifelse(grepl("[[:alnum:]]+C\\d{2}HS", filename), gsub("C\\d{2}HS.*", "", filename) %>% gsub(".*SHA/", "", .), NA),
         box = as.character(box)) %>% 
  distinct()

# check msn
sha_rai_df <- sha_rai_df %>% 
  subset(internal_filename != "MSN:OXYGWAS12h")

# fix NA subjects and duplicated values
sha_rai_df %>% get_dupes(subject, exp) %>% View()

sha_rai_df <- sha_rai_df %>% 
  subset(!(session_duration %in% c("0", "1"))) 

# NA subjects 
sha_rai_df <- sha_rai_df %>% 
  mutate(subject = ifelse(!grepl("[MF]\\d+", subject), NA, subject)) %>% 
  group_by(cohort, box, room) %>% 
  fill(subject, .direction = "downup") %>% 
  ungroup()

# check again 
sha_rai_df %>% get_dupes(subject, exp) %>% View()

# replace NA subjects for animals that died from the start of SHA, don't have any values to fill
sha_rai_df <- sha_rai_df %>% 
  rowwise() %>% 
  mutate(subject = ifelse(grepl("C04HSSHA\\d+", filename)&is.na(subject), paste0("F40", box), subject)) %>% 
  ungroup()

# duplicated values
#### figure this one out 
sha_rai_df <- sha_rai_df %>% 
  add_count(subject, exp) %>% 
  subset(!(n > 1&rewards==0&active==0&inactive==0)) %>%  # remove if bad session
  select(-n)

sha_rai_df %>% get_dupes(subject, exp) %>% View()

sha_rai_df <- sha_rai_df %>% 
  add_count(subject, exp) %>% 
  subset(!(cohort == "C01" & session_duration < 200&n!=1)) %>% 
  mutate(exp = replace(exp, n!=1&cohort == "C07"&startdate == "01/25/19", "SHA05")) %>% 
  select(-n) 

sha_rai_df %>% get_dupes(subject, exp) %>% View()

sha_rai_df <- sha_rai_df %>% 
  add_count(subject, exp) %>% 
  subset(!(n > 1&rewards==0&active==0&inactive==0)) %>%  # remove if bad session
  select(-n) 

sha_rai_df %>% get_dupes(subject, exp) %>% View()

sha_rai_df <- sha_rai_df %>% 
  subset(!(subject == "M759"&session_duration == 120&exp =="SHA05"))

sha_rai_df %>% get_dupes(subject, exp) %>% View()



###### OLD FILES ##############
# label data with... 
sha_old_files <- grep(list.files(path = ".", recursive = T, full.names = T), pattern = ".*Old.*SHA", value = T) # 214 files

# extract data...
sha_rai_old <- lapply(sha_old_files, function(x){
  setwd("~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS")
  internal_filename = fread(paste0("awk '/^ProgramName/{flag=1;next}/ProgramDate/{flag=0}flag' ", "'", x, "' | sed \"s/[[:blank:]]//g\""), fill = T, header = F)
  
  subject = fread(paste0("awk '/^RatNumber/{flag=1;next}/ProgramName/{flag=0}flag' ", "'", x, "' | sed \"s/[[:blank:]]//g\""), fill = T, header = F)
  
  box = fread(paste0("awk '/^BoxNumber/{flag=1;next}/Sessionlength/{flag=0}flag' ", "'", x, "' | sed \"s/[[:blank:]]//g\""), fill = T, header = F)
  
  metadata = cbind(subject = subject, box = box)
  
  # checked that all msn/program name are "SAOneLever"
  rewards = fread(paste0("awk '/^totalRewards/{flag=1;next}/TotalResponses/{flag=0}flag' ", "'", x, "' | sed \"s/[[:blank:]]//g\""), fill = T, header = F)
  to_responses = fread(paste0("awk '/^TotalTOResponses/{flag=1;next}/TotalRspInAct/{flag=0}flag' ", "'", x, "' | sed \"s/[[:blank:]]//g\""), fill = T, header = F)
  inactive = fread(paste0("awk '/^TotalRspInAct/{flag=1;next}/TotalTORspInAct/{flag=0}flag' ", "'", x, "' | sed \"s/[[:blank:]]//g\""), fill = T, header = F)
  
  active = rewards + to_responses
  
  data = list("rewards" = rewards,
              "active" = active,
              "inactive" = inactive
  ) %>% do.call(cbind, .)
  
  
  internal_filename <- cbind(internal_filename = internal_filename, cbind(metadata, data)) %>% 
    rename_all(~ stringr::str_replace_all(., '[.](V1)?', '')) %>% 
    mutate(filename = x )
  
  return(internal_filename)
})

sha_rai_old_df <- sha_rai_old %>% 
  rbindlist(fill = T) %>% 
  mutate(cohort = str_match(filename, "C\\d+") %>% as.character,
         filename = gsub(".*/SHA/", "", filename),
         exp = gsub(".*HS(.*)-.*", "\\1", filename) %>% gsub("-\\d+", "", .),
         room = ifelse(grepl("[[:alnum:]]+C\\d{2}HS", filename), gsub("C\\d{2}HS.*", "", filename) %>% gsub(".*SHA/", "", .), NA),
         box = as.character(box)) %>% 
  mutate(subject = toupper(subject) %>% gsub(".*([MF].*)", "\\1", .))

# qc with...
sha_rai_old_df %>% get_dupes(subject, exp)
sha_rai_old_df %>% subset(!grepl("[MF]\\d+", subject))

# deal with the missing subjects...
sha_rai_old_df <- sha_rai_old_df %>% 
  mutate(subject = ifelse(!grepl("[MF]\\d+", subject), NA, subject)) %>% 
  group_by(cohort, room, box) %>% 
  fill(subject, .direction = "downup") %>% 
  ungroup()

sha_rai_old_df %>% get_dupes(subject, exp)

# fix repeats 
sha_rai_old_df <- sha_rai_old_df %>% 
  add_count(subject, exp) %>% 
  subset(!(n!=1&rewards==0&active==0&inactive==0)) %>% 
  select(-n)

sha_rai_old_df %>% get_dupes(subject, exp)


################################
########## SHA ITI #############
################################

sha_sessions08_10_new <- grep(list.files(path = "~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS", recursive = T, full.names = T), pattern = "SHA(0[89]|10)", value = T) %>% grep("New",., value = T) #  files



########## ITI #################

sha_iti_timebin_new <- lapply(sha_sessions08_10_new, function(x){
  setwd("~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS")
  internal_filename = fread(paste0("awk '/MSN:/{print $1 $2 $3}' '", x, "'"), header = F, fill = T)
  
  subject = fread(paste0("awk '/Subject:/{print $2}' '", x, "'"), header = F, fill = T)
  
  box = fread(paste0("awk '/Box:/{print $2}' '", x, "'"), header = F, fill = T)
  
  date = fread(paste0("awk '/^Start Date:/{print $3}' '", x, "'"), header = F, fill = T)
  enddate = fread(paste0("awk '/^End Date:/{print $3}' '", x, "'"), header = F, fill = T)
  
  time = fread(paste0("awk '/^Start Time:/{print $3}' '", x, "'"), header = F, fill = T)
  endtime = fread(paste0("awk '/^End Time:/{print $3}' '", x, "'"), header = F, fill = T)
  
  session_duration <- cbind(date = date, enddate = enddate) %>%
    cbind(time = time) %>% cbind(endtime = endtime) 
  names(session_duration) <- c("date", "enddate", "time", "endtime")
  session_duration <- session_duration %>% 
    mutate(start_datetime = lubridate::mdy_hms(paste(format(as.Date(date, "%m/%d/%y"), "%m/%d/20%y"), time)),
           end_datetime = lubridate::mdy_hms(paste(format(as.Date(enddate, "%m/%d/%y"), "%m/%d/20%y"), endtime)),
           session_duration = difftime(end_datetime, start_datetime, units = "mins") %>% as.numeric() %>% round(0))
  
  metadata = cbind(subject = subject, box = box) %>% 
    cbind(startdate = date) %>% 
    cbind(session_duration = session_duration$session_duration) %>% 
    cbind(internal_filename)
  
  names(metadata) <- c("subject", "box", "startdate", "session_duration", "internal_filename")
  
  
  
  # metadata = cbind(subject = subject, box = box)
  
  rewards_first10min =  fread(paste0("awk '/^W:/{flag=1;next}/5:/{flag=0}flag' ", "'",  x, "' | awk '/0:/{print $3 + $4}'"), header = F, fill = T)
  rewards_last60min = fread(paste0("awk '/^W:/{flag=1;next}/^Y:/{flag=0}flag' ", "'", x, "'| awk '{ if( $1==\"10:\" || $1==\"15:\" || $1==\"20:\" ) print $2+$3+$4+$5+$6}'"), header = F, fill = T) %>%
    zoo::rollapply(., 3, sum, by = 3)
  
  # to get the median and sd of the iti
  iti = fread(paste0("awk '/^V:/{flag=1;next}/^W:/{flag=0}flag' ", "'", x, "'"), header = F, fill = T) 
  
  data_indices <- grep("^0:$", iti$V1) # start of each subject
  split_data <- split(iti, cumsum(1:nrow(iti) %in% data_indices))
  
  # all timestamps and intertrial times
  iti_by_subject <- lapply(split_data, function(x){
    indexremoved <- x[,-1]
    nonzerorows <- indexremoved[rowSums(indexremoved) > 0, ] # remove excessively trailing 0's 
    processeddata_df <- data.frame(timestamps = as.vector(t(data.matrix(nonzerorows)))) # transpose to get by row
    if(any(processeddata_df$timestamps > 7500)){
      processeddata_df %<>% 
        mutate(bin = cut(timestamps, breaks=seq(from = 0, length.out = 73, by = 300), right = T, labels = seq(from = 1, to = 72, by =1))) %<>% ## for long access 
        dplyr::filter(timestamps != 0)
    }
    else{
      processeddata_df %<>% 
        mutate(bin = cut(timestamps, breaks=seq(from = 0, length.out = 25, by = 300), right = T, labels = seq(from = 1, to = 24, by =1))) %<>% ## for short access 
        dplyr::filter(timestamps != 0) 
    }
    
    processeddata_df <- processeddata_df %>%
      mutate(intertrial_time = lead(timestamps) - timestamps,
             bin = as.character(bin))
    
    return(processeddata_df)
  }) 
  
  # get the rewards within the first hour 
  rewards_first1hr <- iti_by_subject %>% 
    lapply(function(x){
      x %>% 
        subset(bin <= 12) %>% 
        nrow
    }) %>% unlist %>% as.data.frame()
  
  
  # iti median, mean, and sd 
  iti_by_subject_agg <- iti_by_subject %>% 
    lapply(function(x){
      if(nrow(x) > 2){
        x %>% 
          summarize(sha_median_iti = median(intertrial_time, na.rm = T),
                    sha_mean_iti = mean(intertrial_time, na.rm = T),
                    sha_sd_iti = sd(intertrial_time, na.rm = T)) %>% 
          mutate(sha_sd_iti_na = sha_sd_iti)  
      }
      else if(nrow(x) %in% c(1, 2)){
        x %>% 
          summarize(sha_median_iti = median(intertrial_time, na.rm = T),
                    sha_mean_iti = mean(intertrial_time, na.rm = T),
                    sha_sd_iti = sd(intertrial_time, na.rm = T),
                    sha_sd_iti_na = 0)
      }
      else if(nrow(x) == 0){
        x %>% 
          summarize(sha_median_iti = median(intertrial_time, na.rm = T),
                    sha_mean_iti = mean(intertrial_time, na.rm = T),
                    sha_sd_iti = sd(intertrial_time, na.rm = T)) %>% 
          mutate(sha_sd_iti_na = sha_sd_iti)
      }
      
    }) %>% 
    rbindlist()
  
  
  
  data = list("sha_rewards_first10min" = rewards_first10min,
              "sha_rewards_last60min" = rewards_last60min,
              iti_by_subject_agg,
              "sha_rewards_first1hr" = rewards_first1hr
  ) %>% do.call(cbind, .)
  
  data <- cbind(metadata, data) %>% 
    rename_all(~ stringr::str_replace_all(., '[.](V1)?', '')) %>% 
    mutate(filename = x )
  
  return(data)
})

sha_iti_new_df <- sha_iti_timebin_new %>% 
  rbindlist(fill = T) %>% 
  mutate(cohort = str_match(filename, "C\\d+") %>% as.character,
         filename = gsub(".*/SHA/", "", filename),
         exp = gsub(".*HS", "", filename) %>% gsub("-\\d+$", "", .),
         room = ifelse(grepl("[[:alnum:]]+C\\d{2}HS", filename), gsub("C\\d{2}HS.*", "", filename) %>% gsub(".*SHA/", "", .), NA),
         box = as.character(box)) 




################### 
### OLD DIRECTORY


sha_sessions08_10_old <- grep(list.files(path = "~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS", recursive = T, full.names = T), pattern = "SHA(0[89]|10)", value = T) %>% grep("Old",., value = T) #  files

# extract data...
sha_iti_old <- lapply(sha_sessions08_10_old, function(x){
  setwd("~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS")
  internal_filename = fread(paste0("awk '/^ProgramName/{flag=1;next}/ProgramDate/{flag=0}flag' ", "'", x, "' | sed \"s/[[:blank:]]//g\""), fill = T, header = F)
  
  subject = fread(paste0("awk '/^RatNumber/{flag=1;next}/ProgramName/{flag=0}flag' ", "'", x, "' | sed \"s/[[:blank:]]//g\""), fill = T, header = F)
  
  box = fread(paste0("awk '/^BoxNumber/{flag=1;next}/Sessionlength/{flag=0}flag' ", "'", x, "' | sed \"s/[[:blank:]]//g\""), fill = T, header = F)
  
  metadata = cbind(subject = subject, box = box)
  
  # checked that all msn/program name are "SAOneLever"
  rewards_bins = fread(paste0("awk '/^BinRewards\\r/{flag=1;next}/endl/{flag=0}flag' ", "'", x, "'"), fill = T, header = F) # 10 min bins

  data_indices <- grep("^list$", rewards_bins$V1) # start of each subject
  split_data <- split(rewards_bins, cumsum(1:nrow(rewards_bins) %in% data_indices))
  
  # all reward bins
  subject_traits <- lapply(split_data, function(x){
    
    if(nrow(x) != 0){ # if the data frame/vector isn't empty 
    indexremoved <- x[-c(1:2),] # remove excess info
    
    data = data.frame(sha_rewards_first10min = indexremoved[1,] %>% unlist() %>% as.numeric,
                      sha_rewards_first1hr = indexremoved[1:6,] %>% unlist() %>% as.numeric %>% sum(na.rm = T),
                      sha_rewards_last60min = indexremoved[7:12] %>% unlist() %>% as.numeric() %>% sum(na.rm = T))
    
    iti <- 10/as.numeric(unlist(indexremoved)) * 60
    iti[!is.finite(iti)] <- NA
    
    data$sha_median_iti = median(iti, na.rm = T)# from Brent's email -- since we do not have proper time stamp data for all the old comp files, the ITI should be calculated using (bin time)/(bin rewards). So for a 10 min bin with 5 rewards, the “Bin ITI” or average ITI for that bin would be 2 min. With that, we can still use the old comp data for ITI and we will at least be able to have a reasonable stdev value based off the average ITI per bin. 
    data$sha_mean_iti = mean(iti, na.rm = T)
    data$sha_sd_iti = sd(iti, na.rm = T)
    data$sha_sd_iti_na = ifelse(indexremoved %>% unlist %>% as.numeric %>% sum %in% c(1, 2), 0, data$sha_sd_iti)
    
    }
    else{
      data = data.frame(sha_rewards_first10min = NA,
                        sha_rewards_first1hr = NA, 
                        sha_rewards_last60min = NA,
                        sha_median_iti = NA,
                        sha_mean_iti = NA, 
                        sha_sd_iti = NA,
                        sha_sd_iti_na = NA)
    }
    
    return(data)
  }) %>% rbindlist()
  
  internal_filename <- cbind(internal_filename = internal_filename, cbind(metadata, subject_traits)) %>% 
    rename_all(~ stringr::str_replace_all(., '[.](V1)?', '')) %>% 
    mutate(filename = x )
  
  return(internal_filename)
})

sha_iti_old_df <- sha_iti_old %>% 
  rbindlist(fill = T) %>% 
  mutate(cohort = str_match(filename, "C\\d+") %>% as.character,
         filename = gsub(".*/SHA/", "", filename),
         exp = gsub(".*HS(.*)-.*", "\\1", filename) %>% gsub("-\\d+", "", .),
         room = ifelse(grepl("[[:alnum:]]+C\\d{2}HS", filename), gsub("C\\d{2}HS.*", "", filename) %>% gsub(".*SHA/", "", .), NA),
         box = as.character(box)) %>% 
  mutate(subject = toupper(subject) %>% gsub(".*([MF].*)", "\\1", .))

# qc with...
sha_iti_old_df %>% get_dupes(subject, exp)
sha_iti_old_df %>% subset(!grepl("[MF]\\d+", subject))

# fix repeats 
sha_iti_old_df <- sha_iti_old_df %>% 
  add_count(subject, exp) %>% 
  subset(!(n!=1&sha_rewards_first10min==0&sha_rewards_last60min==0&sha_median_iti==0&sha_sd_iti==0)) %>% 
  select(-n)

sha_iti_old_df %>% get_dupes(subject, exp)



### join the ITI dataframes
sha_iti_new_df %>% names
sha_iti_old_df %>% names

sha_gwasiti_df <- bind_rows(sha_iti_new_df, sha_iti_old_df)

sha_gwasiti_df %>% get_dupes(subject, exp)

sha_gwasiti_df <- sha_gwasiti_df %>% 
  add_count(subject, exp) %>% 
  subset(!(subject == 0|n != 1 & session_duration < 5)) %>% 
  select(-n)

sha_gwasiti_df %>% get_dupes(subject, exp)

# take means by labanimalid
sha_gwasiti_df_agg <- sha_gwasiti_df %>% 
  select(subject, matches("^sha")) %>% 
  group_by(subject) %>% 
  summarise_all(~mean(., na.rm = T)) %>% 
  ungroup()

## keep an eye on the C13 ID's, C1317, C1319 subset(!grepl("[MF]", subject))





################################
########## SHA (01-03) ITI #############
################################

sha_sessions01_03_new <- grep(list.files(path = "~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS", recursive = T, full.names = T), pattern = "SHA(0[123])", value = T) %>% grep("New",., value = T) #  files



########## ITI #################

sha_iti_timebin_new_2 <- lapply(sha_sessions01_03_new, function(x){
  setwd("~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS")
  internal_filename = fread(paste0("awk '/MSN:/{print $1 $2 $3}' '", x, "'"), header = F, fill = T)
  
  subject = fread(paste0("awk '/Subject:/{print $2}' '", x, "'"), header = F, fill = T)
  
  box = fread(paste0("awk '/Box:/{print $2}' '", x, "'"), header = F, fill = T)
  
  date = fread(paste0("awk '/^Start Date:/{print $3}' '", x, "'"), header = F, fill = T)
  enddate = fread(paste0("awk '/^End Date:/{print $3}' '", x, "'"), header = F, fill = T)
  
  time = fread(paste0("awk '/^Start Time:/{print $3}' '", x, "'"), header = F, fill = T)
  endtime = fread(paste0("awk '/^End Time:/{print $3}' '", x, "'"), header = F, fill = T)
  
  session_duration <- cbind(date = date, enddate = enddate) %>%
    cbind(time = time) %>% cbind(endtime = endtime) 
  names(session_duration) <- c("date", "enddate", "time", "endtime")
  session_duration <- session_duration %>% 
    mutate(start_datetime = lubridate::mdy_hms(paste(format(as.Date(date, "%m/%d/%y"), "%m/%d/20%y"), time)),
           end_datetime = lubridate::mdy_hms(paste(format(as.Date(enddate, "%m/%d/%y"), "%m/%d/20%y"), endtime)),
           session_duration = difftime(end_datetime, start_datetime, units = "mins") %>% as.numeric() %>% round(0))
  
  metadata = cbind(subject = subject, box = box) %>% 
    cbind(startdate = date) %>% 
    cbind(session_duration = session_duration$session_duration) %>% 
    cbind(internal_filename)
  
  names(metadata) <- c("subject", "box", "startdate", "session_duration", "internal_filename")
  
  
  
  # metadata = cbind(subject = subject, box = box)
  
  rewards_first10min =  fread(paste0("awk '/^W:/{flag=1;next}/5:/{flag=0}flag' ", "'",  x, "' | awk '/0:/{print $3 + $4}'"), header = F, fill = T)
  rewards_last60min = fread(paste0("awk '/^W:/{flag=1;next}/^Y:/{flag=0}flag' ", "'", x, "'| awk '{ if( $1==\"10:\" || $1==\"15:\" || $1==\"20:\" ) print $2+$3+$4+$5+$6}'"), header = F, fill = T) %>%
    zoo::rollapply(., 3, sum, by = 3)
  
  # to get the median and sd of the iti
  iti = fread(paste0("awk '/^V:/{flag=1;next}/^W:/{flag=0}flag' ", "'", x, "'"), header = F, fill = T) 
  
  data_indices <- grep("^0:$", iti$V1) # start of each subject
  split_data <- split(iti, cumsum(1:nrow(iti) %in% data_indices))
  
  # all timestamps and intertrial times
  iti_by_subject <- lapply(split_data, function(x){
    indexremoved <- x[,-1]
    nonzerorows <- indexremoved[rowSums(indexremoved) > 0, ] # remove excessively trailing 0's 
    processeddata_df <- data.frame(timestamps = as.vector(t(data.matrix(nonzerorows)))) # transpose to get by row
    if(any(processeddata_df$timestamps > 7500)){
      processeddata_df %<>% 
        mutate(bin = cut(timestamps, breaks=seq(from = 0, length.out = 73, by = 300), right = T, labels = seq(from = 1, to = 72, by =1))) %<>% ## for long access 
        dplyr::filter(timestamps != 0)
    }
    else{
      processeddata_df %<>% 
        mutate(bin = cut(timestamps, breaks=seq(from = 0, length.out = 25, by = 300), right = T, labels = seq(from = 1, to = 24, by =1))) %<>% ## for short access 
        dplyr::filter(timestamps != 0) 
    }
    
    processeddata_df <- processeddata_df %>%
      mutate(intertrial_time = lead(timestamps) - timestamps,
             bin = as.character(bin))
    
    return(processeddata_df)
  }) 
  
  # get the rewards within the first hour 
  rewards_first1hr <- iti_by_subject %>% 
    lapply(function(x){
      x %>% 
        subset(bin <= 12) %>% 
        nrow
    }) %>% unlist %>% as.data.frame()
  
  
  # iti median, mean, and sd 
  iti_by_subject_agg <- iti_by_subject %>% 
    lapply(function(x){
      if(nrow(x) > 2){
        x %>% 
          summarize(sha_median_iti = median(intertrial_time, na.rm = T),
                    sha_mean_iti = mean(intertrial_time, na.rm = T),
                    sha_sd_iti = sd(intertrial_time, na.rm = T)) %>% 
          mutate(sha_sd_iti_na = sha_sd_iti)  
      }
      else if(nrow(x) %in% c(1, 2)){
        x %>% 
          summarize(sha_median_iti = median(intertrial_time, na.rm = T),
                    sha_mean_iti = mean(intertrial_time, na.rm = T),
                    sha_sd_iti = sd(intertrial_time, na.rm = T),
                    sha_sd_iti_na = 0)
      }
      else if(nrow(x) == 0){
        x %>% 
          summarize(sha_median_iti = median(intertrial_time, na.rm = T),
                    sha_mean_iti = mean(intertrial_time, na.rm = T),
                    sha_sd_iti = sd(intertrial_time, na.rm = T)) %>% 
          mutate(sha_sd_iti_na = sha_sd_iti)
      }
      
    }) %>% 
    rbindlist()
  
  
  
  data = list("sha_rewards_first10min" = rewards_first10min,
              "sha_rewards_last60min" = rewards_last60min,
              iti_by_subject_agg,
              "sha_rewards_first1hr" = rewards_first1hr
  ) %>% do.call(cbind, .)
  
  data <- cbind(metadata, data) %>% 
    rename_all(~ stringr::str_replace_all(., '[.](V1)?', '')) %>% 
    mutate(filename = x )
  
  return(data)
})

sha_iti_new_df_2 <- sha_iti_timebin_new_2 %>% 
  rbindlist(fill = T) %>% 
  mutate(cohort = str_match(filename, "C\\d+") %>% as.character,
         filename = gsub(".*/SHA/", "", filename),
         exp = gsub(".*HS", "", filename) %>% gsub("-\\d+$", "", .),
         room = ifelse(grepl("[[:alnum:]]+C\\d{2}HS", filename), gsub("C\\d{2}HS.*", "", filename) %>% gsub(".*SHA/", "", .), NA),
         box = as.character(box)) 




################### 
### OLD DIRECTORY


sha_sessions01_03_old <- grep(list.files(path = "~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS", recursive = T, full.names = T), pattern = "SHA(0[123])", value = T) %>% grep("Old",., value = T) #  files

# extract data...
sha_iti_old_2 <- lapply(sha_sessions01_03_old, function(x){
  setwd("~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS")
  internal_filename = fread(paste0("awk '/^ProgramName/{flag=1;next}/ProgramDate/{flag=0}flag' ", "'", x, "' | sed \"s/[[:blank:]]//g\""), fill = T, header = F)
  
  subject = fread(paste0("awk '/^RatNumber/{flag=1;next}/ProgramName/{flag=0}flag' ", "'", x, "' | sed \"s/[[:blank:]]//g\""), fill = T, header = F)
  
  box = fread(paste0("awk '/^BoxNumber/{flag=1;next}/Sessionlength/{flag=0}flag' ", "'", x, "' | sed \"s/[[:blank:]]//g\""), fill = T, header = F)
  
  metadata = cbind(subject = subject, box = box)
  
  # checked that all msn/program name are "SAOneLever"
  rewards_bins = fread(paste0("awk '/^BinRewards\\r/{flag=1;next}/endl/{flag=0}flag' ", "'", x, "'"), fill = T, header = F) # 10 min bins
  
  data_indices <- grep("^list$", rewards_bins$V1) # start of each subject
  split_data <- split(rewards_bins, cumsum(1:nrow(rewards_bins) %in% data_indices))
  
  # all reward bins
  subject_traits <- lapply(split_data, function(x){
    
    if(nrow(x) != 0){ # if the data frame/vector isn't empty 
      indexremoved <- x[-c(1:2),] # remove excess info
      
      data = data.frame(sha_rewards_first10min = indexremoved[1,] %>% unlist() %>% as.numeric,
                        sha_rewards_first1hr = indexremoved[1:6,] %>% unlist() %>% as.numeric %>% sum(na.rm = T),
                        sha_rewards_last60min = indexremoved[7:12] %>% unlist() %>% as.numeric() %>% sum(na.rm = T))
      
      iti <- 10/as.numeric(unlist(indexremoved)) * 60
      iti[!is.finite(iti)] <- NA
      
      data$sha_median_iti = median(iti, na.rm = T)# from Brent's email -- since we do not have proper time stamp data for all the old comp files, the ITI should be calculated using (bin time)/(bin rewards). So for a 10 min bin with 5 rewards, the “Bin ITI” or average ITI for that bin would be 2 min. With that, we can still use the old comp data for ITI and we will at least be able to have a reasonable stdev value based off the average ITI per bin. 
      data$sha_mean_iti = mean(iti, na.rm = T)
      data$sha_sd_iti = sd(iti, na.rm = T)
      data$sha_sd_iti_na = ifelse(indexremoved %>% unlist %>% as.numeric %>% sum %in% c(1, 2), 0, data$sha_sd_iti)
      
    }
    else{
      data = data.frame(sha_rewards_first10min = NA,
                        sha_rewards_first1hr = NA, 
                        sha_rewards_last60min = NA,
                        sha_median_iti = NA,
                        sha_mean_iti = NA, 
                        sha_sd_iti = NA,
                        sha_sd_iti_na = NA)
    }
    
    return(data)
  }) %>% rbindlist()
  
  internal_filename <- cbind(internal_filename = internal_filename, cbind(metadata, subject_traits)) %>% 
    rename_all(~ stringr::str_replace_all(., '[.](V1)?', '')) %>% 
    mutate(filename = x )
  
  return(internal_filename)
})

sha_iti_old_df_2 <- sha_iti_old_2 %>% 
  rbindlist(fill = T) %>% 
  mutate(cohort = str_match(filename, "C\\d+") %>% as.character,
         filename = gsub(".*/SHA/", "", filename),
         exp = gsub(".*HS(.*)-.*", "\\1", filename) %>% gsub("-\\d+", "", .),
         room = ifelse(grepl("[[:alnum:]]+C\\d{2}HS", filename), gsub("C\\d{2}HS.*", "", filename) %>% gsub(".*SHA/", "", .), NA),
         box = as.character(box)) %>% 
  mutate(subject = toupper(subject) %>% gsub(".*([MF].*)", "\\1", .))

# qc with...
sha_iti_old_df_2 %>% get_dupes(subject, exp)
sha_iti_old_df_2 %>% subset(!grepl("[MF]\\d+", subject))

# fix repeats 
sha_iti_old_df_2 <- sha_iti_old_df_2 %>% 
  add_count(subject, exp) %>% 
  subset(!(n!=1&sha_rewards_first10min==0&sha_rewards_last60min==0&sha_median_iti==0&sha_sd_iti==0)) %>% 
  select(-n)

sha_iti_old_df %>% get_dupes(subject, exp)



### join the ITI dataframes
sha_iti_new_df_2 %>% names
sha_iti_old_df_2 %>% names

sha_gwasiti_df_2 <- bind_rows(sha_iti_new_df_2, sha_iti_new_df_2)

sha_gwasiti_df_2 <- sha_gwasiti_df_2 %>% 
  mutate(cohort = str_match(filename, "C\\d+") %>% as.character,
         filename = gsub(".*/SHA/", "", filename),
         exp = gsub(".*HS(.*)-.*", "\\1", filename) %>% gsub("-\\d+", "", .),
         room = ifelse(grepl("[[:alnum:]]+C\\d{2}HS", filename), gsub("C\\d{2}HS.*", "", filename) %>% gsub(".*SHA/", "", .), NA),
         box = as.character(box)) %>% 
  mutate(subject = toupper(subject) %>% gsub(".*([MF].*)", "\\1", .)) %>% 
  distinct

  
sha_gwasiti_df_2 %>% get_dupes(subject, exp)

sha_gwasiti_df_2 <- sha_gwasiti_df_2 %>% 
  add_count(subject, exp) %>% 
  subset(!(subject == 0|n != 1 & session_duration < 5)) %>% 
  select(-n)

sha_gwasiti_df_2 %>% get_dupes(subject, exp)

# take means by labanimalid
sha_gwasiti_df_agg_2 <- sha_gwasiti_df_2 %>% 
  select(subject, matches("^sha")) %>% 
  group_by(subject) %>% 
  summarise_all(~mean(., na.rm = T)) %>% 
  ungroup() %>% 
  rename_all(~ stringr::str_replace_all(., '(.*)', '\\1_2')) %>% 
  rename_all(~ stringr::str_replace_all(., '_2_2', '_2')) %>% 
  rename("subject" = "subject_2")

## keep an eye on the C13 ID's, C1317, C1319 subset(!grepl("[MF]", subject))




################################
########## LGA #################
################################


###### NEW FILES ##############
# label data with... 
lga_new_files_01_11 <- grep(grep(list.files(path = "~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS", recursive = T, full.names = T), pattern = ".*txt", inv = T, value = T), pattern = ".*LGA", value = T) %>% grep("C(0[1-9]|1[0-1])", ., value = T) # 464 files

# check cohorts and experiments represented 
lga_new_files_01_11 %>% as.data.frame() %>% 
  mutate(cohort = str_extract(., "C\\d+"),
         exp = str_extract(gsub("-\\d+$", "", .), "LGA(\\d+)$")) %>% 
  select(cohort, exp) %>% table(exclude = NULL)

# check which ones are NA
lga_new_files_01_11 %>% as.data.frame() %>% 
  mutate(cohort = str_extract(., "C\\d+"),
         exp = str_extract(gsub("-\\d+$", "", .), "LGA(\\d+)$")) %>% 
  subset(is.na(exp))


# extract subject, box, and rewards, active, inactive (rai)
lga_rai <- lapply(lga_new_files_01_11, function(x){
  setwd("~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS")
  internal_filename = fread(paste0("awk '/MSN:/{print $1 $2 $3}' '", x, "'"), header = F, fill = T)
  
  subject = fread(paste0("awk '/Subject:/{print $2}' '", x, "'"), header = F, fill = T)
  
  box = fread(paste0("awk '/Box:/{print $2}' '", x, "'"), header = F, fill = T)
  
  date = fread(paste0("awk '/^Start Date:/{print $3}' '", x, "'"), header = F, fill = T)
  enddate = fread(paste0("awk '/^End Date:/{print $3}' '", x, "'"), header = F, fill = T)
  
  time = fread(paste0("awk '/^Start Time:/{print $3}' '", x, "'"), header = F, fill = T)
  endtime = fread(paste0("awk '/^End Time:/{print $3}' '", x, "'"), header = F, fill = T)
  
  session_duration <- cbind(date = date, enddate = enddate) %>%
    cbind(time = time) %>% cbind(endtime = endtime) 
  names(session_duration) <- c("date", "enddate", "time", "endtime")
  session_duration <- session_duration %>% 
    mutate(start_datetime = lubridate::mdy_hms(paste(format(as.Date(date, "%m/%d/%y"), "%m/%d/20%y"), time)),
           end_datetime = lubridate::mdy_hms(paste(format(as.Date(enddate, "%m/%d/%y"), "%m/%d/20%y"), endtime)),
           session_duration = difftime(end_datetime, start_datetime, units = "mins") %>% as.numeric() %>% round(0))
  
  metadata = cbind(subject = subject, box = box) %>% 
    cbind(startdate = date) %>% 
    cbind(session_duration = session_duration$session_duration) %>% 
    cbind(internal_filename)
  
  names(metadata) <- c("subject", "box", "startdate", "session_duration", "internal_filename")
  
  
  # metadata = cbind(subject = subject, box = box)
  
  rewards = fread(paste0("awk '/^B: /{print $2}' ", "'", x, "'"), header = F, fill = T)
  active = fread(paste0("awk '/^G: /{print $2}' ", "'", x, "'"), header = F, fill = T)
  inactive = fread(paste0("awk '/^A: /{print $2}' ", "'", x, "'"), header = F, fill = T)
  
  data = list("rewards" = rewards,
              "active" = active,
              "inactive" = inactive  
  ) %>% do.call(cbind, .)
  
  data <- cbind(metadata, data) %>% 
    rename_all(~ stringr::str_replace_all(., '[.](V1)?', '')) %>% 
    mutate(filename = x )
  
  return(data)
})

lga_rai_df <- lga_rai %>% 
  rbindlist(fill = T) %>% 
  mutate(cohort = str_match(filename, "C\\d+") %>% as.character,
         filename = gsub(".*/LGA/", "", filename),
         exp = gsub(".*HS", "", filename) %>% gsub("-\\d+$", "", .),
         room = ifelse(grepl("[[:alnum:]]+C\\d{2}HS", filename), gsub("C\\d{2}HS.*", "", filename) %>% gsub(".*LGA/", "", .), NA),
         box = as.character(box)) 

# qc with...
lga_rai_df %>% get_dupes(subject, exp, cohort) 

# deal with the missing subjects...
lga_rai_df <- lga_rai_df %>% 
  mutate(subject = ifelse(!grepl("[MF]", subject), NA, subject)) %>% 
  group_by(cohort, room, box) %>% 
  fill(subject, .direction = "downup") %>% 
  ungroup()

# remove duplicates if measurements are all 0
lga_rai_df <- lga_rai_df %>% 
  add_count(subject, exp, cohort) %>% 
  subset(!(n!=1&rewards==0&active==0&inactive==0)) %>% 
  select(-n)

lga_rai_df %>% get_dupes(subject, exp, cohort) 

# correct duplicates because of error filenames

lga_rai_df <- lga_rai_df %>% 
  mutate(filename = replace(filename, cohort == "C07"&startdate == "02/11/19", "MED1110C07HSLGA05"),
         exp  = replace(exp, cohort == "C07"&startdate == "02/11/19", "LGA05")) %>% 
  subset(!(filename == "BSB273DC11HSLGA02"&startdate == "03/05/20"))




###### OLD FILES ##############
lga_old_files <- grep(list.files(path = "~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS", recursive = T, full.names = T), pattern = ".*Old.*LGA", value = T) # 424 files

# label data with... 

lga_rai_old <- lapply(lga_old_files, function(x){
  setwd("~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS")
  internal_filename = fread(paste0("awk '/^ProgramName/{flag=1;next}/ProgramDate/{flag=0}flag' ", "'", x, "' | sed \"s/[[:blank:]]//g\""), fill = T, header = F)
  
  subject = fread(paste0("awk '/^RatNumber/{flag=1;next}/ProgramName/{flag=0}flag' ", "'", x, "' | sed \"s/[[:blank:]]//g\""), fill = T, header = F)
  
  box = fread(paste0("awk '/^BoxNumber/{flag=1;next}/Sessionlength/{flag=0}flag' ", "'", x, "' | sed \"s/[[:blank:]]//g\""), fill = T, header = F)
  
  metadata = cbind(subject = subject, box = box)
  
  # checked that all msn/program name are "SAOneLever"
  rewards = fread(paste0("awk '/^BinRewards/{flag=1;next}/TotalResponses/{flag=0}flag' ", "'", x, "' | sed \"s/[[:blank:]]//g\""), fill = T, header = F)
  to_responses = fread(paste0("awk '/^TotalTOResponses/{flag=1;next}/TotalRspInAct/{flag=0}flag' ", "'", x, "' | sed \"s/[[:blank:]]//g\""), fill = T, header = F)
  inactive = fread(paste0("awk '/^TotalRspInAct/{flag=1;next}/TotalTORspInAct/{flag=0}flag' ", "'", x, "' | sed \"s/[[:blank:]]//g\""), fill = T, header = F)
  
  active = rewards + to_responses
  
  data = list("rewards" = rewards,
              "active" = active,
              "inactive" = inactive
  ) %>% do.call(cbind, .)
  
  
  internal_filename <- cbind(internal_filename = internal_filename, cbind(metadata, data)) %>% 
    rename_all(~ stringr::str_replace_all(., '[.](V1)?', '')) %>% 
    mutate(filename = x )
  
  return(internal_filename)
})

lga_rai_old_df <- lga_rai_old %>% 
  rbindlist(fill = T) %>% 
  mutate(cohort = str_match(filename, "C\\d+") %>% as.character,
         filename = gsub(".*/LGA/", "", filename),
         exp = gsub(".*HS(.*)-.*", "\\1", filename) %>% gsub("-\\d+", "", .),
         room = ifelse(grepl("[[:alnum:]]+C\\d{2}HS", filename), gsub("C\\d{2}HS.*", "", filename) %>% gsub(".*LGA/", "", .), NA),
         box = as.character(box)) %>% 
  mutate(subject = toupper(subject))

# qc with...
lga_rai_old_df %>% get_dupes(subject, exp)
lga_rai_old_df %>% subset(!grepl("[MF]\\d+", subject))

# assign wrong subjects 
lga_rai_old_df <- lga_rai_old_df %>% 
  mutate(subject = ifelse(filename == "Q3C03HSLGA01-20180221.txt", paste0("M", subject), 
                          ifelse(filename == "K3C04HSLGA12-20180516.txt"&box == "4", "M459",
                                 ifelse(filename == "K2C04HSLGA15-20180524.txt"&box == "2", "M451", subject))))

# deal with the missing subjects...
lga_rai_old_df <- lga_rai_old_df %>% 
  mutate(subject = ifelse(!grepl("[MF]\\d+", subject), NA, subject)) %>% 
  group_by(cohort, room, box) %>% 
  fill(subject, .direction = "downup") %>% 
  ungroup()

lga_rai_old_df %>% get_dupes(subject, exp)

# fix repeats 
lga_rai_old_df <- lga_rai_old_df %>% 
  add_count(subject, exp, cohort) %>% 
  subset(!(n!=1&rewards==0&active==0&inactive==0)) %>% 
  select(-n)

lga_rai_old_df %>% get_dupes(subject, exp)

lga_rai_old_df <- lga_rai_old_df %>% 
  add_count(subject, exp) %>% 
  rowwise() %>% 
  mutate(subject = ifelse(n!=1&parse_number(cohort)!= substr(as.character(parse_number(subject)), 1, 1)&parse_number(cohort)>3, gsub("\\d(\\d+)$", "\\1" %>% paste0(as.character(parse_number(cohort)),.), subject), subject)) %>% 
  select(-n) %>% 
  ungroup()

lga_rai_old_df <- lga_rai_old_df %>% 
  subset(!(subject == "M460"&rewards == 1&active==1))

lga_rai_old_df %>% get_dupes(subject, exp)



################################
########## LGA ITI #############
################################


# NEW DIRECTORY 
lga_sessions12_14_new <- grep(list.files(path = "~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS", recursive = T, full.names = T), pattern = "LGA(1[2-4])", value = T) %>% grep("New",., value = T) #  files



########## ITI #################

lga_iti_timebin <- lapply(lga_sessions12_14_new, function(x){
  setwd("~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS")
  internal_filename = fread(paste0("awk '/MSN:/{print $1 $2 $3}' '", x, "'"), header = F, fill = T)
  
  subject = fread(paste0("awk '/Subject:/{print $2}' '", x, "'"), header = F, fill = T)
  
  box = fread(paste0("awk '/Box:/{print $2}' '", x, "'"), header = F, fill = T)
  
  date = fread(paste0("awk '/^Start Date:/{print $3}' '", x, "'"), header = F, fill = T)
  enddate = fread(paste0("awk '/^End Date:/{print $3}' '", x, "'"), header = F, fill = T)
  
  time = fread(paste0("awk '/^Start Time:/{print $3}' '", x, "'"), header = F, fill = T)
  endtime = fread(paste0("awk '/^End Time:/{print $3}' '", x, "'"), header = F, fill = T)
  
  session_duration <- cbind(date = date, enddate = enddate) %>%
    cbind(time = time) %>% cbind(endtime = endtime) 
  names(session_duration) <- c("date", "enddate", "time", "endtime")
  session_duration <- session_duration %>% 
    mutate(start_datetime = lubridate::mdy_hms(paste(format(as.Date(date, "%m/%d/%y"), "%m/%d/20%y"), time)),
           end_datetime = lubridate::mdy_hms(paste(format(as.Date(enddate, "%m/%d/%y"), "%m/%d/20%y"), endtime)),
           session_duration = difftime(end_datetime, start_datetime, units = "mins") %>% as.numeric() %>% round(0))
  
  metadata = cbind(subject = subject, box = box) %>% 
    cbind(startdate = date) %>% 
    cbind(session_duration = session_duration$session_duration) %>% 
    cbind(internal_filename)
  
  names(metadata) <- c("subject", "box", "startdate", "session_duration", "internal_filename")
  
  
  
  # metadata = cbind(subject = subject, box = box)
  
  rewards_first10min =  fread(paste0("awk '/^W:/{flag=1;next}/5:/{flag=0}flag' ", "'",  x, "' | awk '/0:/{print $3 + $4}'"), header = F, fill = T)
  rewards_last60min = fread(paste0("awk '/^W:/{flag=1;next}/^Y:/{flag=0}flag' ", "'", x, "'| awk '{ if( $1==\"60:\" || $1==\"65:\" || $1==\"70:\" ) print $2+$3+$4+$5+$6}'"), header = F, fill = T) %>%
    zoo::rollapply(., 3, sum, by = 3)
  
  # to get the median and sd of the iti
  iti = fread(paste0("awk '/^V:/{flag=1;next}/^W:/{flag=0}flag' ", "'", x, "'"), header = F, fill = T) 
  
  data_indices <- grep("^0:$", iti$V1) # start of each subject
  split_data <- split(iti, cumsum(1:nrow(iti) %in% data_indices))
  
  # all timestamps and intertrial times
  iti_by_subject <- lapply(split_data, function(x){
    indexremoved <- x[,-1]
    nonzerorows <- indexremoved[rowSums(indexremoved) > 0, ] # remove excessively trailing 0's 
    processeddata_df <- data.frame(timestamps = as.vector(t(data.matrix(nonzerorows)))) # transpose to get by row
    if(any(processeddata_df$timestamps > 7500)){
      processeddata_df %<>% 
        mutate(bin = cut(timestamps, breaks=seq(from = 0, length.out = 73, by = 300), right = T, labels = seq(from = 1, to = 72, by =1))) %<>% ## for long access 
        dplyr::filter(timestamps != 0)
    }
    else{
      processeddata_df %<>% 
        mutate(bin = cut(timestamps, breaks=seq(from = 0, length.out = 25, by = 300), right = T, labels = seq(from = 1, to = 24, by =1))) %<>% ## for short access 
        dplyr::filter(timestamps != 0) 
    }
    
    processeddata_df <- processeddata_df %>%
      mutate(intertrial_time = lead(timestamps) - timestamps,
             bin = as.character(bin))
    
    return(processeddata_df)
  }) 
  
  # get the rewards within the first hour 
  rewards_first1hr <- iti_by_subject %>% 
    lapply(function(x){
      x %>% 
        subset(bin <= 12) %>% 
        nrow
    }) %>% unlist %>% as.data.frame()
  
  # iti median, mean, and sd 
  iti_by_subject_agg <- iti_by_subject %>% 
    lapply(function(x){
      if(nrow(x) > 2){
        x %>% 
          summarize(lga_median_iti = median(intertrial_time, na.rm = T),
                    lga_mean_iti = mean(intertrial_time, na.rm = T),
                    lga_sd_iti = sd(intertrial_time, na.rm = T)) %>% 
          mutate(lga_sd_iti_na = lga_sd_iti)  
      }
      else if(nrow(x) %in% c(1, 2)){
        x %>% 
          summarize(lga_median_iti = median(intertrial_time, na.rm = T),
                    lga_mean_iti = mean(intertrial_time, na.rm = T),
                    lga_sd_iti = sd(intertrial_time, na.rm = T),
                    lga_sd_iti_na = 0)
      }
      else if(nrow(x) == 0){
        x %>% 
          summarize(lga_median_iti = median(intertrial_time, na.rm = T),
                    lga_mean_iti = mean(intertrial_time, na.rm = T), 
                    lga_sd_iti = sd(intertrial_time, na.rm = T)) %>% 
          mutate(lga_sd_iti_na = lga_sd_iti)
      }
      
    }) %>% 
    rbindlist()
  
  
  
  data = list("lga_rewards_first10min" = rewards_first10min,
              "lga_rewards_last60min" = rewards_last60min,
              iti_by_subject_agg,
              "lga_rewards_first1hr" = rewards_first1hr
  ) %>% do.call(cbind, .)
  
  data <- cbind(metadata, data) %>% 
    rename_all(~ stringr::str_replace_all(., '[.](V1)?', '')) %>% 
    mutate(filename = x )
  
  return(data)
})

lga_iti_new_df <- lga_iti_timebin %>% 
  rbindlist(fill = T) %>% 
  mutate(cohort = str_match(filename, "C\\d+") %>% as.character,
         filename = gsub(".*/LGA/", "", filename),
         exp = gsub(".*HS", "", filename) %>% gsub("-\\d+$", "", .),
         room = ifelse(grepl("[[:alnum:]]+C\\d{2}HS", filename), gsub("C\\d{2}HS.*", "", filename) %>% gsub(".*LGA/", "", .), NA),
         box = as.character(box)) 




# OLD DIRECTORY 
lga_sessions12_14_old <- grep(list.files(path = "~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS", recursive = T, full.names = T), pattern = "LGA(1[2-4])", value = T) %>% grep("Old",., value = T) #  files


# extract data...
lga_iti_old <- lapply(lga_sessions12_14_old, function(x){
  setwd("~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS")
  internal_filename = fread(paste0("awk '/^ProgramName/{flag=1;next}/ProgramDate/{flag=0}flag' ", "'", x, "' | sed \"s/[[:blank:]]//g\""), fill = T, header = F)
  
  subject = fread(paste0("awk '/^RatNumber/{flag=1;next}/ProgramName/{flag=0}flag' ", "'", x, "' | sed \"s/[[:blank:]]//g\""), fill = T, header = F)
  
  box = fread(paste0("awk '/^BoxNumber/{flag=1;next}/Sessionlength/{flag=0}flag' ", "'", x, "' | sed \"s/[[:blank:]]//g\""), fill = T, header = F)
  
  metadata = cbind(subject = subject, box = box)
  
  # checked that all msn/program name are "SAOneLever"
  rewards_bins = fread(paste0("awk '/^BinRewards\\r/{flag=1;next}/endl/{flag=0}flag' ", "'", x, "'"), fill = T, header = F) # 10 min bins
  
  data_indices <- grep("^list$", rewards_bins$V1) # start of each subject
  split_data <- split(rewards_bins, cumsum(1:nrow(rewards_bins) %in% data_indices))
  
  # all reward bins
  subject_traits <- lapply(split_data, function(x){
    
    if(nrow(x) != 0){
      indexremoved <- x[-c(1:2),] # remove excess
      
      data = data.frame(lga_rewards_first10min = indexremoved[1] %>% unlist() %>% as.numeric,
                        lga_rewards_first1hr = indexremoved[1:6] %>% unlist() %>% as.numeric %>% sum(na.rm = T),
                        lga_rewards_last60min = indexremoved[31:36] %>% unlist() %>% as.numeric() %>% sum(na.rm = T))
      
      iti <- 10/as.numeric(unlist(indexremoved)) * 60
      iti[!is.finite(iti)] <- NA
      
      data$lga_median_iti = median(iti, na.rm = T)# from Brent's email -- since we do not have proper time stamp data for all the old comp files, the ITI should be calculated using (bin time)/(bin rewards). So for a 10 min bin with 5 rewards, the “Bin ITI” or average ITI for that bin would be 2 min. With that, we can still use the old comp data for ITI and we will at least be able to have a reasonable stdev value based off the average ITI per bin. 
      data$lga_mean_iti = mean(iti, na.rm = T)
      data$lga_sd_iti = sd(iti, na.rm = T)
      data$lga_sd_iti_na = ifelse(indexremoved %>% unlist %>% as.numeric %>% sum %in% c(1, 2), 0, data$lga_sd_iti)
      
    }
    else{
      data = data.frame(lga_rewards_first10min = NA,
                        lga_rewards_first1hr = NA,
                        lga_rewards_last60min = NA,
                        lga_median_iti = NA,
                        lga_mean_iti = NA,
                        lga_sd_iti = NA,
                        lga_sd_iti_na = NA)
    }
    
    return(data)
  }) %>% rbindlist()
  
  internal_filename <- cbind(internal_filename = internal_filename, cbind(metadata, subject_traits)) %>% 
    rename_all(~ stringr::str_replace_all(., '[.](V1)?', '')) %>% 
    mutate(filename = x )
  
  return(internal_filename)
})

lga_iti_old_df <- lga_iti_old %>% 
  rbindlist(fill = T) %>% 
  mutate(cohort = str_match(filename, "C\\d+") %>% as.character,
         filename = gsub(".*/lga/", "", filename),
         exp = gsub(".*HS(.*)-.*", "\\1", filename) %>% gsub("-\\d+", "", .),
         room = ifelse(grepl("[[:alnum:]]+C\\d{2}HS", filename), gsub("C\\d{2}HS.*", "", filename) %>% gsub(".*lga/", "", .), NA),
         box = as.character(box)) %>% 
  mutate(subject = toupper(subject) %>% gsub(".*([MF].*)", "\\1", .))

# qc with...
lga_iti_old_df %>% get_dupes(subject, exp)
lga_iti_old_df %>% subset(!grepl("[MF]\\d+", subject))

# fix repeats 
lga_iti_old_df <- lga_iti_old_df %>% 
  add_count(subject, exp) %>% 
  subset(!(n!=1&lga_rewards_first10min==0&lga_rewards_last60min==0&lga_median_iti==0&lga_sd_iti==0)) %>% 
  subset(subject != "0") %>% 
  select(-n)

lga_iti_old_df %>% get_dupes(subject, exp)



### join the ITI dataframes
lga_iti_new_df %>% names
lga_iti_old_df %>% names

lga_gwasiti_df <- bind_rows(lga_iti_new_df, lga_iti_old_df)

lga_gwasiti_df %>% get_dupes(subject, exp)

lga_gwasiti_df <- lga_gwasiti_df %>% 
  add_count(subject, exp) %>% 
  subset(!(subject == 0|n != 1 & session_duration < 5)) %>% 
  select(-n)

lga_gwasiti_df %>% get_dupes(subject, exp)

lga_gwasiti_df <- lga_gwasiti_df %>% 
  select(-box, -session_duration) %>% # inconsistent box info
  distinct()

lga_gwasiti_df %>% get_dupes(subject, exp)


# take means by labanimalid
lga_gwasiti_df_agg <- lga_gwasiti_df %>% 
  select(subject, matches("^lga")) %>% 
  group_by(subject) %>% 
  summarise_all(~mean(., na.rm = T)) %>% 
  ungroup()






################################
########## PR ##################
################################

###### NEW FILES ##############

# label data with... 
pr_new_files_01_11 <- grep(grep(list.files(path = "~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS", recursive = T, full.names = T), pattern = ".*txt", inv = T, value = T), pattern = ".*PR/", value = T) %>% grep("C(0[1-9]|1[0-1])", ., value = T)

# extract subject, box, and rewards, active, inactive, pr 
pr_rai <- lapply(pr_new_files_01_11, function(x){
  setwd("~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS")
  internal_filename = fread(paste0("awk '/MSN:/{print $1 $2 $3}' '", x, "'"), header = F, fill = T)
  
  subject = fread(paste0("awk '/Subject:/{print $2}' '", x, "'"), header = F, fill = T)
  
  box = fread(paste0("awk '/Box:/{print $2}' '", x, "'"), header = F, fill = T)
  
  metadata = cbind(subject = subject, box = box)
  
  if(grepl("GWAS", unique(internal_filename$V1))){
    rewards = fread(paste0("awk '/^B: /{print $2}' ", "'", x, "'"), header = F, fill = T)
    active = fread(paste0("awk '/^G: /{print $2}' ", "'", x, "'"), header = F, fill = T)
    inactive = fread(paste0("awk '/^A: /{print $2}' ", "'", x, "'"), header = F, fill = T)
    pr <- fread(paste0("awk '/^S: /{print $2}' ", "'", x, "'"), header = F, fill = T)
  }
  else{
    rewards = fread(paste0("awk '/^C: /{print $2}' ", "'", x, "'"), header = F, fill = T)
    active = fread(paste0("awk '/^B: /{print $2}' ", "'", x, "'"), header = F, fill = T)
    inactive = fread(paste0("awk '/^A: /{print $2}' ", "'", x, "'"), header = F, fill = T)
    pr <- fread(paste0("awk '/^L: /{print $2}' ", "'", x, "'"), header = F, fill = T)
  }
  
  data = list("rewards" = rewards,
              "active" = active,
              "inactive" = inactive, 
              "pr" = pr
  ) %>% do.call(cbind, .)
  
  
  internal_filename <- cbind(msn = internal_filename, cbind(metadata, data)) %>% 
    rename_all(~ stringr::str_replace_all(., '[.](V1)?', '')) %>% 
    mutate(filename = x )
  
  return(internal_filename)
})

pr_rai_df <- pr_rai %>% 
  rbindlist(fill = T) %>% 
  mutate(cohort = str_match(filename, "C\\d+") %>% as.character,
         filename = gsub(".*/PR/", "", filename),
         exp = gsub(".*HS", "", filename),
         room = ifelse(grepl("[[:alnum:]]+C\\d{2}HS", filename), gsub("C\\d{2}HS.*", "", filename) %>% gsub(".*PR/", "", .), NA),
         box = as.character(box)) 


# qc with...
pr_rai_df <- pr_rai_df %>% subset(subject != "717")
pr_rai_df %>% count(subject, exp, cohort) %>% subset(n!=1)
pr_rai_df %>% get_dupes(subject, exp, cohort)
pr_rai_df %>% distinct() %>% add_count(subject, exp, cohort) %>% subset(n!=1)

# deal with the missing subjects...
pr_rai_df <- pr_rai_df %>% 
  mutate(subject = ifelse(!grepl("[MF]", subject), NA, subject)) %>% 
  group_by(cohort, room, box) %>% 
  fill(subject, .direction = "downup") %>% 
  ungroup()

pr_rai_df <- pr_rai_df %>% subset(grepl("[MF]", subject))

pr_rai_df %>% get_dupes(subject, exp)



###### OLD FILES ##############

pr_old_files <- grep(list.files(path = ".", recursive = T, full.names = T), pattern = ".*Old.*PR/", value = T) # 62 files

# extract subject, box, and rewards, active, inactive, pr 







################################
########## SHOCK ###############
################################

###### NEW FILES ##############
shock_new_files_c01_11 <- grep(list.files(path = "~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS", recursive = T, full.names = T), pattern = ".*New.*SHOCK/", value = T) %>% grep("C(0[1-9]|1[0-1])", ., value = T) # 54 files

# check cohorts and experiments represented 
shock_new_files_c01_11 %>% as.data.frame() %>% 
  mutate(cohort = str_extract(., "C\\d+"),
         exp = str_extract(gsub("-\\d+$", "", .), "(PRE)?SHOCK(\\d+)?$")) %>% 
  select(cohort, exp) %>% table(exclude = NULL)

# check which ones are NA
shock_new_files_c01_11 %>% as.data.frame() %>% 
  mutate(cohort = str_extract(., "C\\d+"),
         exp = str_extract(gsub("-\\d+$", "", .), "(PRE)?SHOCK(\\d+)?$"))%>% 
  subset(is.na(exp))

# extract subject, box, and rewards, active, inactive (rai)
shock_rai <- lapply(shock_new_files_c01_11, function(x){
  setwd("~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS")
  internal_filename = fread(paste0("awk '/MSN:/{print $1 $2 $3}' '", x, "'"), header = F, fill = T)
  subject = fread(paste0("awk '/Subject:/{print $2}' '", x, "'"), header = F, fill = T)
  
  box = fread(paste0("awk '/Box:/{print $2}' '", x, "'"), header = F, fill = T)
  
  date = fread(paste0("awk '/^Start Date:/{print $3}' '", x, "'"), header = F, fill = T)
  enddate = fread(paste0("awk '/^End Date:/{print $3}' '", x, "'"), header = F, fill = T)
  
  time = fread(paste0("awk '/^Start Time:/{print $3}' '", x, "'"), header = F, fill = T)
  endtime = fread(paste0("awk '/^End Time:/{print $3}' '", x, "'"), header = F, fill = T)
  
  session_duration <- cbind(date = date, enddate = enddate) %>%
    cbind(time = time) %>% cbind(endtime = endtime) 
  names(session_duration) <- c("date", "enddate", "time", "endtime")
  session_duration <- session_duration %>% 
    mutate(start_datetime = lubridate::mdy_hms(paste(format(as.Date(date, "%m/%d/%y"), "%m/%d/20%y"), time)),
           end_datetime = lubridate::mdy_hms(paste(format(as.Date(enddate, "%m/%d/%y"), "%m/%d/20%y"), endtime)),
           session_duration = difftime(end_datetime, start_datetime, units = "mins") %>% as.numeric() %>% round(0))
  
  metadata = cbind(subject = subject, box = box) %>% 
    cbind(startdate = date) %>% 
    cbind(session_duration = session_duration$session_duration) %>% 
    cbind(internal_filename)
  
  names(metadata) <- c("subject", "box", "startdate", "session_duration", "internal_filename")
  
  rewards = fread(paste0("awk '/^B: /{print $2}' ", "'", x, "'"), header = F, fill = T)
  active = fread(paste0("awk '/^G: /{print $2}' ", "'", x, "'"), header = F, fill = T)
  inactive = fread(paste0("awk '/^A: /{print $2}' ", "'", x, "'"), header = F, fill = T)
  firstshock = fread(paste0("awk '/^M:/{flag=1;next}/5:/{flag=0}flag' ", "'", x, "' | awk '/0:/{print $2}'"), header = F, fill = T)
  
  data = list("rewards" = rewards,
              "active" = active,
              "inactive" = inactive,
              "firstshock" = firstshock
  ) %>% do.call(cbind, .)
  
  data <- cbind(metadata, data) %>%
    rename_all(~ stringr::str_replace_all(., '[.](V1)?', '')) %>%
    mutate(filename = x ) %>% 
    mutate(rewards_w_postshock1 = ifelse(rewards != 0&firstshock !=0, rewards-firstshock+1, NA))
  
  return(data)
})

shock_rai_df <- shock_rai %>% 
  rbindlist(fill = T) %>% 
  mutate(cohort = str_match(filename, "C\\d+") %>% as.character,
         subject = str_extract(toupper(subject), "[MF]\\d+"),
         filename = gsub(".*/SHOCK/", "", filename),
         exp = gsub(".*HS", "", filename) %>% gsub("-\\d+$", "", .),
         exp = ifelse(grepl("PRE", exp), str_extract(exp, "PRESHOCK"), 
                      ifelse(parse_number(cohort) > 6&!grepl("PRE", exp), "SHOCK03", exp)), #clean up exp
         room = ifelse(grepl("[[:alnum:]]+C\\d{2}HS", filename), gsub("C\\d{2}HS.*", "", filename) %>% gsub(".*SHOCK/", "", .), NA),
         box = as.character(box)) %>% 
  distinct()



# qc with...
shock_rai_df %>% get_dupes(subject, exp, cohort) 

# deal with the missing subjects...
shock_rai_df <- shock_rai_df %>% 
  group_by(cohort, room, box) %>% 
  fill(subject, .direction = "downup") %>% 
  ungroup()

shock_rai_df <- shock_rai_df %>% 
  add_count(subject, exp, cohort) %>% 
  subset(!(n!=1&rewards==0&active==0&inactive==0)) %>% 
  select(-n)

shock_rai_df %>% get_dupes(subject, exp, cohort) 

# fix duplicates 
shock_rai_df <- shock_rai_df %>% 
  mutate(subject = ifelse(subject == "F7"&box == "7"&exp == "SHOCK03"&rewards == 34, "F16", subject))


###### OLD FILES ##############
## no old files for SHOCK


###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################

## new code setup -- run by cohort
# 
# setwd("~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS")
# # olivier_directory_names <- grep("C\\d+/.*/", list.dirs(), value = T)
# 
# # make function for the entire cohort
# # if new directories, if old directories
# # and specify the exp
# read_exp_by_cohort <- function(x){
#   # read the new dir (sha)
#   # pattern = 
#   sha_files <- grep(list.files(path = ".", recursive = T, full.names = T), pattern = ".*C01.*New.*SHA", value = T) 
#   
#   # read the old dir (sha)
#   # pattern = 
#     
# }
# 
# 
#   
#   
# sha_new_files <- grep(grep(list.files(path = ".", recursive = T, full.names = T), pattern = ".*txt", inv = T, value = T), pattern = ".*SHA", value = T) # 254 files
# sha_subjects_new <- process_subjects_new(sha_new_files) %>% separate(labanimalid, c("row", "labanimalid"), sep = "_", extra = "merge") %>% 
#   arrange(filename, as.numeric(row)) %>% select(-c(row, filename))
# read_rewards_new <- function(x){
#   rewards <- fread(paste0("awk '/W:/{flag=1;next}/5:/{flag=0}flag' ", "'", x, "' | awk '/0:/{print NR \"_\" $2}'"), header = F, fill = T)
#   rewards$filename <- x
#   return(rewards)
# }
# sha_rewards_new <-  lapply(sha_new_files, read_rewards_new) %>% rbindlist() %>% separate(V1, into = c("row", "rewards"), sep = "_") %>% arrange(filename, as.numeric(row)) %>% select(-row) %>% 
#   bind_cols(sha_subjects_new) %>% 
#   separate(labanimalid, into = c("labanimalid", "cohort", "exp", "filename", "date", "time", "box"), sep = "_") %>% 
#   mutate(date = lubridate::mdy(date), time = chron::chron(times = time)) %>%  
#   left_join(., date_time_subject_df_comp %>% 
#               select(cohort, exp, filename, valid, start_date, start_time, exp_dur_min) %>% 
#               rename("date" = "start_date", "time" = "start_time"), 
#             by = c("cohort", "exp", "filename", "date", "time")) 
# 
# ## reconsider since the data may not be updated from the date_time_subject_df_comp object
# sha_rewards_new_valid <- sha_rewards_new %>% 
#   dplyr::filter(valid == "yes") %>%
#   # dplyr::filter(valid == "yes"|is.na(valid)) %>% ## XX ALLOW INTO THE CODE ONCE WE HAVE THE EXCEL SHEETS FOR C10 AND C11
#   mutate(time = as.character(time)) %>%
#   dplyr::filter(!filename %in% c("C01HSSHA06", "MED1113C07HSSHA05", "MED1114C07HSSHA08")) %>% # update records several lines down from meeting to show other team's confirmation 
#   distinct() # fixes duplicates in filenames %in% c("MED1113C07HSSHA06", "MED1110C05HSSHA08", "MED1110C05HSSHA09") ### there are no dupes for dplyr::filter(!grepl("[MF]\\d+", labanimalid)) 
# 
# ## notes 
# ## exclude files (from meeting)
# # c("C01HSSHA06", "MED1113C07HSSHA05", "MED1114C07HSSHA08")
# ## exclude cases (from meeting )
# # c("F720") for SHA03 bc both files with her data seem incorrect (MED1112C07HSSHA03 and MED1112C07HSSHA03-2)
# # MED1113C07HSSHA07 is actually LGA data (code that validates the date is filtering out these cases, and in the file, sha07 data and pr data follows)
# 
# # deal with the missing subjects...
# # join and update "df" by reference, i.e. without copy 
# 
# ## ADDED _valid 5/20 -- remove once unneeded 
# setDT(sha_rewards_new_valid)             # convert to data.table without copy
# sha_rewards_new_valid[setDT(sha_rewards_new_valid %>% dplyr::filter(!grepl("[MF]", labanimalid)) %>% # this captures all "NA" cases as checked with mutate_at(vars(labanimalid), na_if, "NA") %>% dplyr::filter(is.na(labanimalid))
#                               left_join(., date_time_subject_df_comp %>% 
#                                           select(labanimalid, cohort, exp, filename, start_date, start_time, exp_dur_min) %>% 
#                                           rename("date" = "start_date", "time" = "start_time") %>% 
#                                           mutate(time = as.character(time)), 
#                                         by = c("cohort", "exp", "filename", "date", "time")) ), 
#                       on = c("rewards", "cohort", "exp", "filename", "date", "time", "valid"), labanimalid := labanimalid.y] # don't want to make another missing object
# setDF(sha_rewards_new_valid)
# sha_rewards_new_valid %<>% 
#   mutate_at(vars(rewards), as.numeric)
# ## case: deal with mislabelled subject?
# sha_rewards_new_valid %>% count(labanimalid, cohort,exp) %>% subset(n != 1)
# sha_rewards_new_valid %<>% mutate(labanimalid = replace(labanimalid, exp=="SHA01"&time=="09:24:16", "M768")) ## extracted as M7678 from file, but verified to have the same box (box 16)
# 
# 
# 
# ###### OLD FILES ##############
# # label data with... 
# sha_subjects_old <- process_subjects_old(sha_old_files)
# # extract data...
# sha_rewards_old <- lapply(sha_old_files, read_fread_old, "rewards") %>% rbindlist() %>% separate(V1, into = c("row", "rewards"), sep = "_") %>% arrange(filename, as.numeric(row)) %>% select(-row) %>% 
#   bind_cols(sha_subjects_old %>% arrange(filename, as.numeric(row)) %>% select(-c("row", "filename"))) %>% 
#   separate(labanimalid, into = c("labanimalid", "box", "cohort", "exp", "computer", "date", "valid"), sep = "_") %>% 
#   mutate(date = lubridate::ymd(date),
#          rewards = rewards %>% as.numeric()) %>% 
#   dplyr::filter(valid == "valid") # no need for distinct() bc it is not an issue here
# 
# # deal with the missing subjects...
# sha_rewards_old %>% dplyr::filter(!grepl("[MF]", labanimalid)) %>% dim
# # will remove these cases bc these files have 7 subjects and both misssing subjects have another "session" (matched box)
# sha_rewards_old %<>% dplyr::filter(grepl("[MF]", labanimalid)) 
# 
# ## case: deal with mislabelled subject?
# sha_rewards_old %>% add_count(labanimalid, cohort,exp) %>% subset(n != 1)
# sha_rewards_old %<>% add_count(labanimalid, cohort,exp) %<>% dplyr::filter(n == 1|(n==2&rewards!=0)) %<>% select(-n)
# 
# sha_rewards_old %>% get_dupes(labanimalid, cohort,exp)
# 
# 
# 
# 



### Extract timeout presses

#####
## lga
#####

lga_c01_11_files <- list.files(path = "~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS", recursive = T, full.names = T) %>% 
  grep("LGA", ., value = T)


lga_c01_11_timeout <- lapply(lga_c01_11_files, function(x){
  if(grepl("Old", x)){ # process if directory is old
    df <- fread(paste0("awk '/^RatNumber/{flag=1;next}/^ProgramName/{flag=0}flag' ", "'", x, "'"), fill = T, header = F) 
    df$filename = x
    df$box = fread(paste0("awk '/^BoxNumber/{flag=1;next}/^Sessionlength/{flag=0}flag' ", "'", x, "'"), fill = T, header = F)  
    df$to_active_presses = fread(paste0("awk '/^TotalTOResponses/{flag=1;next}/^TotalRspInAct/{flag=0}flag' ", "'", x, "'"), fill = T, header = F)
  }
  
  if(grepl("New", x)){ # process if directory is new
    df <- fread(paste0("awk '/Subject/{print NR \"_\" $2}' ", "'", x, "'"), fill = T, header = F)
    df$filename = x 
    df$rewards = fread(paste0("awk '/B:/{print $2}' ", "'", x, "'"), fill = T, header = F)
    df$presses = fread(paste0("awk '/G:/{print $2}' ", "'", x, "'"), fill = T, header = F)
    df$box = fread(paste0("awk '/Box:/{print $2}' ", "'", x, "'"), fill = T, header = F) 
    df$to_active_presses =  df$presses -  df$rewards
    df <- df[, c("V1", "filename", "box", "to_active_presses")]
  }
  
  return(df)
  
})

# use to try for old  lga_c01_11_files[c(19:22, 538:541)]
# use to try for new  lga_c01_11_files[c(1:4, 885:888)]
# use to try for old and new  lga_c01_11_files[c(1:4, 885:888, 19:22, 538:541)]

lga_c01_11_timeout_df <- lga_c01_11_timeout %>% rbindlist() %>% 
  rename("labanimalid" = "V1") %>% 
  mutate(labanimalid = str_extract(toupper(labanimalid), "[MF]\\d+"),
         session = str_extract(toupper(filename), "LGA\\d+"),
         cohort = str_extract(toupper(filename), "/C\\d+/") %>% gsub("/", "", .),
         sex = str_extract(toupper(labanimalid), "[MF]"),
         room = ifelse(grepl("[[:alnum:]]+C\\d{2}HS", filename), gsub("C\\d{2}HS.*", "", filename) %>% gsub(".*LGA/", "", .), NA)
  ) %>% 
  select(cohort, labanimalid, sex, session, box, room, to_active_presses, filename) 

lga_c01_11_timeout_trials1_14 <- lga_c01_11_timeout_df %>% 
  subset(parse_number(session) < 15) %>% 
  mutate(box = as.character(box))

# fix these 
lga_c01_11_timeout_trials1_14 %>% get_dupes(labanimalid, session)
# fix the na subject with boxes 
lga_c01_11_timeout_trials1_14_join <- lga_c01_11_timeout_trials1_14 %>% 
  left_join(box_metadata_long %>% 
              select(cohort, exp, box, room, labanimalid) %>% 
              rename("session" = "exp"), by = c("cohort", "room", "box", "session")) %>% 
  mutate(labanimalid = coalesce(labanimalid.x, labanimalid.y)) %>% 
  select(-labanimalid.x, -labanimalid.y) 
# lga_c01_11_timeout_trials1_14_join %>% get_dupes(labanimalid, session) %>% naniar::vis_miss
lga_c01_11_timeout_trials1_14_distinct <- lga_c01_11_timeout_trials1_14_join %>% 
  select(cohort, labanimalid, sex, session, box, room, to_active_presses) %>% 
  distinct() %>% 
  left_join(lga_c01_11_timeout_trials1_14_join %>% 
              select(cohort, labanimalid, sex, session, box, room, to_active_presses) %>% 
              distinct() %>% 
              get_dupes(labanimalid, session) %>% 
              group_by(labanimalid, session, .drop = F) %>% 
              mutate(nzero = sum(to_active_presses)) %>% 
              subset(to_active_presses != nzero&to_active_presses == 0) %>% ungroup() %>% mutate(drop = "drop") %>% 
              select(labanimalid, cohort, session, to_active_presses, drop), by = c("labanimalid", "cohort", "session", "to_active_presses")) %>% 
  subset(is.na(drop))
lga_c01_11_timeout_trials1_14_brent <- lga_c01_11_timeout_trials1_14_distinct %>% 
  left_join(., lga_c01_11_timeout_trials1_14_distinct %>% get_dupes(labanimalid, session) %>% select(labanimalid, cohort, session, to_active_presses, dupe_count), by = c("labanimalid", "cohort", "session", "to_active_presses")) %>% 
  rename("need_check"="dupe_count") %>% 
  mutate(need_check = ifelse(!is.na(need_check), "check", NA)) %>% 
  select(cohort, labanimalid, sex, session, box, room, to_active_presses, need_check) %>% 
  mutate(labanimalid_num = parse_number(labanimalid)) %>% 
  arrange(cohort, sex, labanimalid_num, session) %>% 
  select(-labanimalid_num)


write.xlsx(lga_c01_11_timeout_trials1_14_brent, "~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE/cocaine_lga_to_presses.xlsx")

# use brent's file to correct 
lga_c01_11_timeout_trials1_14_brent_decision <- lga_c01_11_timeout_trials1_14_brent %>% 
  cbind(openxlsx::read.xlsx("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE/cocaine_lga_to_presses_BB.xlsx")) %>% 
  clean_names() %>% 
  mutate_all(as.character) %>% 
  mutate(cohort = coalesce(cohort, cohort_2),
         labanimalid = coalesce(labanimalid, labanimalid_2),
         sex = coalesce(sex, sex_2), 
         session = coalesce(session, session_2), 
         box = coalesce(box, box_2), 
         room = coalesce(room, room_2), 
         to_active_presses = coalesce(to_active_presses, to_active_presses_2)) %>% 
  select(-matches("_2$")) %>% 
  subset(is.na(decision)|decision=="KEEP") %>% 
  rowwise() %>% 
  mutate(labanimalid = replace(labanimalid, cohort == "C07"&grepl("M3", labanimalid), gsub("M3", "M7", labanimalid))) %>% 
  mutate(session = replace(session, grepl("typo LGA", notes), str_extract(notes,"LGA\\d+")),
         labanimalid = replace(labanimalid, grepl("typo [FM]\\d+", notes), str_extract(notes, "[FM]\\d+"))) %>% 
  # mutate(notes = replace(notes, grepl("blank replaced with corr(r)?ect animal ID, toss repeats", notes)|grepl("typo", notes), NA)) %>%  # after fixing
  ungroup() 
# %>% 
# select(-decision, -notes) # after fixing

lga_c01_11_timeout_trials1_14_brent_decision %>% get_dupes(labanimalid, session) %>% View()  
# provide graphs for olivier
lga_c01_11_timeout_trials1_14_brent_decision %>% mutate(to_active_presses = as.numeric(to_active_presses)) %>% 
  ggplot(aes(x = to_active_presses)) + geom_histogram(bins = 1000) + facet_wrap(~ cohort)
lga_c01_11_timeout_trials1_14_brent_decision %>% mutate(to_active_presses = as.numeric(to_active_presses)) %>% 
  ggplot(aes(x = to_active_presses)) + geom_density() 
lga_c01_11_timeout_trials1_14_brent_decision %>% mutate(to_active_presses = as.numeric(to_active_presses)) %>% 
  ggplot(aes(x = cohort, y = to_active_presses)) + geom_boxplot()

# lga_c01_11_timeout_trials1_14_brent_decision <- lga_c01_11_timeout_trials1_14_brent_decision %>% 
#   mutate(labanimalid)

# mutate()
# mutate_at(vars(-decision, -notes, -matches("_2")), ~ . = coalesce(., select(., !! matches("_2"))))

lga_c01_11_timeout_trials1_14_brent_decision_2 <- lga_c01_11_timeout_trials1_14_brent_decision %>% 
  select(cohort, labanimalid, sex, session, box, room, to_active_presses) %>% 
  distinct() %>% 
  subset(!(labanimalid == "M459"&session == "LGA12"&to_active_presses == 0)) %>% 
  subset(!(labanimalid == "M751"&session == "LGA03"&to_active_presses == 3))

# spread, mean, and deviation 
lga_c01_11_timeout_trials1_14_spread <- lga_c01_11_timeout_trials1_14_brent_decision_2 %>% 
  spread(session, to_active_presses) %>% mutate_at(vars(matches("^LGA")), as.numeric) %>% mutate(mean = rowMeans(select(., matches("^LGA")), na.rm = T)) %>% mutate_at(vars(matches("^LGA")), list(dev = ~round((.-mean)^2))) %>% 
  mutate(labanimalid_num = parse_number(labanimalid)) %>% 
  arrange(cohort, sex, labanimalid_num) %>% 
  select(-labanimalid_num)

openxlsx::write.xlsx(lga_c01_11_timeout_trials1_14_spread, "~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE/cocaine_lga_to_presses_wide_dev.xlsx")

lga_c01_11_timeout_trials1_14_final <- openxlsx::read.xlsx("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE/cocaine_lga_to_presses_wide_dev_BB.xlsx") %>% 
  clean_names() %>% 
  mutate_all(as.character) %>% 
  select(-matches("mean$|_dev$")) %>% 
  gather("session", "to_active_presses", -cohort, -labanimalid, -sex, -box, -room, -omit) %>% 
  mutate(session = toupper(session)) %>% 
  rowwise() %>% 
  mutate(omit = gsub("; ", "|", omit), 
         to_active_presses = replace(to_active_presses, grepl(omit, session), NA)) %>% 
  ungroup() %>%  
  select(-omit) %>%
  spread(session, to_active_presses) %>% 
  mutate(labanimalid_num = parse_number(labanimalid)) %>% 
  arrange(cohort, sex, labanimalid_num) %>% 
  select(-labanimalid_num) 

openxlsx::write.xlsx(lga_c01_11_timeout_trials1_14_final, "~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE/cocaine_lga_to_presses_final.xlsx")



#####
## sha
#####
sha_c01_11_files <- list.files(path = "~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS", recursive = T, full.names = T) %>% 
  grep("SHA", ., value = T)


sha_c01_11_timeout <- lapply(sha_c01_11_files, function(x){
  if(grepl("Old", x)){ # process if directory is old
    df <- fread(paste0("awk '/^RatNumber/{flag=1;next}/^ProgramName/{flag=0}flag' ", "'", x, "'"), fill = T, header = F) 
    df$filename = x
    df$box = fread(paste0("awk '/^BoxNumber/{flag=1;next}/^Sessionlength/{flag=0}flag' ", "'", x, "'"), fill = T, header = F)  
    df$to_active_presses = fread(paste0("awk '/^TotalTOResponses/{flag=1;next}/^TotalRspInAct/{flag=0}flag' ", "'", x, "'"), fill = T, header = F)
  }
  
  if(grepl("New", x)){ # process if directory is new
    df <- fread(paste0("awk '/Subject/{print NR \"_\" $2}' ", "'", x, "'"), fill = T, header = F)
    df$filename = x 
    df$rewards = fread(paste0("awk '/W:/{flag=1;next}/5:/{flag=0}flag' ", "'", x, "' | awk '/0:/{print $2}'"), header = F, fill = T)
    df$presses = fread(paste0("awk '/R:/{flag=1;next}/5:/{flag=0}flag' ", "'", x, "' | awk '/0:/{print $2}'"), header = F, fill = T)
    df$box = fread(paste0("awk '/Box:/{print $2}' ", "'", x, "'"), fill = T, header = F) 
    df$to_active_presses =  df$presses -  df$rewards
    df <- df[, c("V1", "filename", "box", "to_active_presses")]
  }
  
  return(df)
  
})

# use to test old  sha_c01_11_files[c(10:13, 268:271)]
# use to test new  sha_c01_11_files[c(1:4, 465:468)]
# use to test old and new  sha_c01_11_files[c(1:4, 465:468, 10:13, 268:271))]

sha_c01_11_timeout_df <- sha_c01_11_timeout %>% rbindlist() %>% 
  rename("labanimalid" = "V1") %>% 
  mutate(labanimalid = str_extract(toupper(labanimalid), "[MF]\\d+"),
         session = str_extract(toupper(filename), "SHA\\d+"),
         cohort = str_extract(toupper(filename), "/C\\d+/") %>% gsub("/", "", .),
         sex = str_extract(toupper(labanimalid), "[MF]"),
         room = ifelse(grepl("[[:alnum:]]+C\\d{2}HS", filename), gsub("C\\d{2}HS.*", "", filename) %>% gsub(".*SHA/", "", .), NA)
  ) %>% 
  mutate(box = as.character(box)) %>% 
  select(cohort, labanimalid, sex, session, room, box, to_active_presses, filename)


# fix these cases
sha_c01_11_timeout_df %>% get_dupes(labanimalid, session)
# fix the na subject with boxes 
sha_c01_11_timeout_join <- sha_c01_11_timeout_df %>% 
  left_join(box_metadata_long %>% 
              select(cohort, exp, box, room, labanimalid) %>% 
              rename("session" = "exp"), by = c("cohort", "room", "box", "session")) %>% 
  mutate(labanimalid = coalesce(labanimalid.x, labanimalid.y)) %>% 
  select(-labanimalid.x, -labanimalid.y) 

sha_c01_11_timeout_distinct <- sha_c01_11_timeout_join %>% 
  select(cohort, labanimalid, sex, session, box, room, to_active_presses) %>% 
  distinct() %>% 
  left_join(sha_c01_11_timeout_join %>% 
              select(cohort, labanimalid, sex, session, box, room, to_active_presses) %>% 
              distinct() %>% 
              get_dupes(labanimalid, session) %>% 
              group_by(labanimalid, session, .drop = F) %>% 
              mutate(nzero = sum(to_active_presses)) %>% 
              subset(to_active_presses != nzero&to_active_presses == 0) %>% ungroup() %>% mutate(drop = "drop") %>% 
              select(labanimalid, cohort, session, to_active_presses, drop), by = c("labanimalid", "cohort", "session", "to_active_presses")) %>% 
  subset(is.na(drop))
sha_c01_11_timeout_brent <- sha_c01_11_timeout_distinct %>% 
  left_join(., sha_c01_11_timeout_distinct %>% get_dupes(labanimalid, session) %>% select(labanimalid, cohort, session, to_active_presses, dupe_count), by = c("labanimalid", "cohort", "session", "to_active_presses")) %>% 
  rename("need_check"="dupe_count") %>% 
  mutate(need_check = ifelse(!is.na(need_check), "check", NA)) %>% 
  select(cohort, labanimalid, sex, session, box, room, to_active_presses, need_check) %>% 
  mutate(labanimalid_num = parse_number(labanimalid)) %>% 
  arrange(cohort, sex, labanimalid_num, session) %>% 
  select(-labanimalid_num) %>% 
  mutate(need_check = replace(need_check, is.na(labanimalid), "check"))

openxlsx::write.xlsx(sha_c01_11_timeout_brent, "~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE/cocaine_sha_to_presses.xlsx")

# use brent's file to correct 
sha_c01_11_timeout_brent_decision <- sha_c01_11_timeout_brent %>% 
  cbind(openxlsx::read.xlsx("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE/cocaine_sha_to_presses_BB.xlsx")) %>% 
  clean_names() %>% 
  mutate_all(as.character) %>% 
  mutate(cohort = coalesce(cohort, cohort_2),
         labanimalid = coalesce(labanimalid, labanimalid_2),
         sex = coalesce(sex, sex_2), 
         session = coalesce(session, session_2), 
         box = coalesce(box, box_2), 
         room = coalesce(room, room_2), 
         to_active_presses = coalesce(to_active_presses, to_active_presses_2)) %>% 
  select(-matches("_2$")) %>% 
  subset(is.na(decision)|decision=="KEEP") %>% 
  rowwise() %>% 
  mutate(session = replace(session, grepl("typo SHA", notes), str_extract(notes,"SHA\\d+")),
         labanimalid = replace(labanimalid, grepl("typo [FM]\\d+", notes), str_extract(notes, "[FM]\\d+"))) %>% 
  # mutate(notes = replace(notes, grepl("blank replaced with corr(r)?ect animal ID, toss repeats", notes)|grepl("typo", notes), NA)) %>%  # after fixing
  ungroup() 
# %>% 
# select(-decision, -notes) # after fixing

sha_c01_11_timeout_brent_decision %>% get_dupes(labanimalid, session) %>% View()  


sha_c01_11_timeout_brent_decision %>% get_dupes(labanimalid, session) %>% View()  

sha_c01_11_timeout_brent_decision_2 <- sha_c01_11_timeout_brent_decision %>% 
  select(cohort, labanimalid, sex, session, box, room, to_active_presses) %>% 
  distinct() %>% 
  subset(!(labanimalid == "F2"&session == "SHA06"&to_active_presses == 1)) %>% 
  subset(!(labanimalid == "F817"&session == "SHA04"&to_active_presses == 22)) %>% 
  subset(!(labanimalid == "F98"&session == "SHA01"&to_active_presses == 0)) %>% 
  subset(!(labanimalid == "M354"&session == "SHA10"&box == "6"))

# write the data in wide form and calculate mean and deviations 
sha_c01_11_timeout_brent_decision_dev <- sha_c01_11_timeout_brent_decision_2 %>% 
  spread(session, to_active_presses) %>% mutate_at(vars(matches("^SHA")), as.numeric) %>% mutate(mean = rowMeans(select(., matches("^SHA")), na.rm = T)) %>% mutate_at(vars(matches("^SHA")), list(dev = ~round((.-mean)^2))) %>% 
  mutate(labanimalid_num = parse_number(labanimalid)) %>% 
  arrange(cohort, sex, labanimalid_num) %>% 
  select(-labanimalid_num)

openxlsx::write.xlsx(sha_c01_11_timeout_brent_decision_dev, "~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE/cocaine_sha_to_presses_wide_dev.xlsx")


# provide graphs for olivier
sha_c01_11_timeout_brent_decision %>% mutate(to_active_presses = as.numeric(to_active_presses)) %>% 
  ggplot(aes(x = to_active_presses)) + geom_histogram(bins = 1000) + facet_wrap(~ cohort)
sha_c01_11_timeout_brent_decision %>% mutate(to_active_presses = as.numeric(to_active_presses)) %>% 
  ggplot(aes(x = to_active_presses)) + geom_density() 
sha_c01_11_timeout_brent_decision %>% mutate(to_active_presses = as.numeric(to_active_presses)) %>% 
  ggplot(aes(x = cohort, y = to_active_presses)) + geom_boxplot()


sha_c01_11_timeout_final <- openxlsx::read.xlsx("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE/cocaine_sha_to_presses_wide_dev_BB.xlsx") %>% 
  clean_names() %>% 
  mutate_all(as.character) %>% 
  select(-matches("mean$|_dev$")) %>% 
  gather("session", "to_active_presses", -cohort, -labanimalid, -sex, -box, -room, -omit) %>% 
  mutate(session = toupper(session)) %>% 
  rowwise() %>% 
  mutate(omit = gsub("; ", "|", omit), 
         to_active_presses = replace(to_active_presses, grepl(omit, session), NA)) %>% 
  ungroup() %>%  
  select(-omit) %>%
  spread(session, to_active_presses) %>% 
  mutate(labanimalid_num = parse_number(labanimalid)) %>% 
  arrange(cohort, sex, labanimalid_num) %>% 
  select(-labanimalid_num) 

openxlsx::write.xlsx(sha_c01_11_timeout_final, "~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE/cocaine_sha_to_presses_final.xlsx")









#####
## shock
#####
shock_c01_11_files <- list.files(path = "~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS", recursive = T, full.names = T) %>% 
  grep("SHOCK", ., value = T) %>% 
  grep("PRESHOCK", ., value = T, invert = T) %>% 
  grep("SHOCK0[12]", ., value = T, invert = T)


shock_c01_11_timeout <- lapply(shock_c01_11_files, function(x){
  
  df <- fread(paste0("awk '/Subject/{print NR \"_\" $2}' ", "'", x, "'"), fill = T, header = F)
  df$filename = x 
  df$rewards = fread(paste0("awk '/B:/{print $2}' ", "'", x, "'"), fill = T, header = F)
  df$presses = fread(paste0("awk '/G:/{print $2}' ", "'", x, "'"), fill = T, header = F)
  df$box = fread(paste0("awk '/Box:/{print $2}' ", "'", x, "'"), fill = T, header = F) 
  df$to_active_presses =  df$presses -  df$rewards
  df <- df[, c("V1", "filename", "box", "to_active_presses")]
  return(df)
  
})

# small list, no need to subset 

shock_c01_11_timeout_df <- shock_c01_11_timeout %>% rbindlist() %>% 
  rename("labanimalid" = "V1") %>% 
  mutate(labanimalid = str_extract(toupper(labanimalid), "[MF]\\d+"),
         session = "SHOCK03",
         cohort = str_extract(toupper(filename), "/C\\d+/") %>% gsub("/", "", .),
         sex = str_extract(toupper(labanimalid), "[MF]"),
         room = ifelse(grepl("[[:alnum:]]+C\\d{2}HS", filename), gsub("C\\d{2}HS.*", "", filename) %>% gsub(".*SHOCK/", "", .), NA)
  ) %>% 
  mutate(box = as.character(box)) %>% 
  select(cohort, labanimalid, sex, session, box, room, to_active_presses, filename)

# fix these cases
shock_c01_11_timeout_df %>% get_dupes(labanimalid)
# fix the na subject with boxes 
shock_c01_11_timeout_df_join <- shock_c01_11_timeout_df %>% 
  left_join(box_metadata_long %>% 
              select(cohort, exp, box, room, labanimalid) %>% 
              rename("session" = "exp"), by = c("cohort", "room", "box", "session")) %>% 
  mutate(labanimalid = coalesce(labanimalid.x, labanimalid.y)) %>% 
  select(-labanimalid.x, -labanimalid.y) 
shock_c01_11_timeout_distinct <- shock_c01_11_timeout_df_join %>% 
  select(cohort, labanimalid, sex, session, box, room, to_active_presses) %>% 
  distinct() %>% 
  left_join(shock_c01_11_timeout_df_join %>% 
              select(cohort, labanimalid, sex, session, box, room, to_active_presses) %>% 
              distinct() %>% 
              get_dupes(labanimalid, session) %>% 
              group_by(labanimalid, session, .drop = F) %>% 
              mutate(nzero = sum(to_active_presses)) %>% 
              subset(to_active_presses != nzero&to_active_presses == 0) %>% ungroup() %>% mutate(drop = "drop") %>% 
              select(labanimalid, cohort, session, to_active_presses, drop), by = c("labanimalid", "cohort", "session", "to_active_presses")) %>% 
  subset(is.na(drop))
shock_c01_11_timeout_brent <- shock_c01_11_timeout_distinct %>% 
  left_join(., shock_c01_11_timeout_distinct %>% get_dupes(labanimalid, session) %>% select(labanimalid, cohort, session, to_active_presses, dupe_count), by = c("labanimalid", "cohort", "session", "to_active_presses")) %>% 
  rename("need_check"="dupe_count") %>% 
  mutate(need_check = ifelse(!is.na(need_check), "check", NA)) %>% 
  select(cohort, labanimalid, sex, session, box, room, to_active_presses, need_check) %>% 
  mutate(labanimalid_num = parse_number(labanimalid)) %>% 
  arrange(cohort, sex, labanimalid_num, session) %>% 
  select(-labanimalid_num) %>% 
  distinct()

openxlsx::write.xlsx(shock_c01_11_timeout_brent, "~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE/cocaine_shock_to_presses.xlsx")


# use brent's file to correct 
shock_c01_11_timeout_brent_decision <- shock_c01_11_timeout_brent %>% 
  cbind(openxlsx::read.xlsx("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE/cocaine_shock_to_presses_BB.xlsx")) %>% 
  clean_names() %>% 
  mutate_all(as.character) %>% 
  mutate(cohort = coalesce(cohort, cohort_2),
         labanimalid = coalesce(labanimalid, labanimalid_2),
         sex = coalesce(sex, sex_2), 
         session = coalesce(session, session_2), 
         box = coalesce(box, box_2), 
         room = coalesce(room, room_2), 
         to_active_presses = coalesce(to_active_presses, to_active_presses_2)) %>% 
  select(-matches("_2$")) %>% 
  subset(is.na(decision)|decision=="KEEP") %>% 
  rowwise() %>% 
  mutate(session = replace(session, grepl("typo SHOCK", notes), str_extract(notes,"SHOCK\\d+")),
         labanimalid = replace(labanimalid, grepl("typo [FM]\\d+", notes), str_extract(notes, "[FM]\\d+"))) %>% 
  # mutate(notes = replace(notes, grepl("blank replaced with corr(r)?ect animal ID, toss repeats", notes)|grepl("typo", notes), NA)) %>%  # after fixing
  ungroup() 
# %>% 
# select(-decision, -notes) # after fixing

shock_c01_11_timeout_brent_decision %>% get_dupes(labanimalid, session) %>% View()  

shock_c01_11_timeout_brent_decision_2 <- shock_c01_11_timeout_brent_decision %>% 
  select(cohort, labanimalid, sex, session, box, room, to_active_presses) %>% 
  distinct()

shock_c01_11_timeout_brent_decision_2 %>% get_dupes(labanimalid, session) %>% View()  

# provide graphs for olivier
shock_c01_11_timeout_brent_decision %>% mutate(to_active_presses = as.numeric(to_active_presses)) %>% 
  ggplot(aes(x = to_active_presses)) + geom_histogram(bins = 1000) + facet_wrap(~ cohort)
shock_c01_11_timeout_brent_decision %>% mutate(to_active_presses = as.numeric(to_active_presses)) %>% 
  ggplot(aes(x = to_active_presses)) + geom_density() 
shock_c01_11_timeout_brent_decision %>% mutate(to_active_presses = as.numeric(to_active_presses)) %>% 
  ggplot(aes(x = cohort, y = to_active_presses)) + geom_boxplot()


# write the data in wide form and calculate mean and deviations 
shock_c01_11_timeout_brent_decision_2 <- shock_c01_11_timeout_brent_decision_2 %>% 
  spread(session, to_active_presses) %>% mutate_at(vars(matches("^SHA")), as.numeric) %>% mutate(mean = rowMeans(select(., matches("^SHA")))) %>% mutate_at(vars(matches("^SHA")), list(dev = ~round((.-mean)^2)))
openxlsx::write.xlsx(sha_c01_11_timeout_brent_decision_dev, "~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE/cocaine_lga_to_presses.xlsx")

