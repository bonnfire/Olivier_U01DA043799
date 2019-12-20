## extract raw 
setwd("~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS")

# after cohort 5, there are only new files

# extract the new files 

## USEFUL FUNCTIONS
# #extract names to be assigned for various tables later
readsubjects <- function(x){
  subjects <- fread(paste0("awk '/Subject/{print $2}' ", "'", x, "'"),fill = T,header=F)
  subjects$filename <- x
  return(subjects)
}

# #notice system call, doing it recurisvely for the entire directory; for these cases, the group number is not 0 and the subject is found elsewhere.  # 12/18 only lga
groups_files <- system("grep -ir \"Group:\" | grep -v \"Group: 0\"", intern = TRUE) %>% 
  gsub("\r", "", .) %>% 
  as.data.frame() %>% 
  rename("filename" = ".") %>% 
  separate(filename, c("filename", "labanimalid"), sep = ":", extra = "merge") %>%  # only split specified number of times
  mutate(labanimalid = paste0(str_match(labanimalid, "[FM]\\d{1,3}"), "_", str_extract(filename, "C\\d+"), "_", sub('.*HS', '', toupper(filename)), "_", sub(".*/.*/.*/", '', groups_files$filename) ), 
         comment = "labanimalid extracted from group in raw files") %>% 
  select(-filename)

readboxes <- function(x){
  boxes <- fread(paste0("awk '/Box/{print $2}' ", "'", x, "'"),fill = T,header=F)
  boxes$filename <- x
  return(boxes)
}

# to differentiate the bad sessions from the good ones; find the typical amount of time spent on a session and filter
readdate_time <- function(x){
  
}

#### ONE FUNCTIONS EDITION #####
# create the dataframe with vector in it
# check that the column we are removing ends with 5 or 0 and then remove

# write a df containing the fread statements and function for extracting the different dfs
read_fread <- function(x, varname){
  
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
  # return(split_data)
  
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
  
  # names(processeddata) <- grep("C01LGA01", names_append, value = T)
  
  return(processeddata)
}
## know how many subjects to expect in each filename
## find . -name "SHA" -exec grep -ira1 "NumberOfSubjects" {} +

read_subject_old <- function(x){
    subject_old <- fread(paste0("grep -iEo \"[F|M][0-9]+\" ", "'", x, "'"), header = F)
  subject_old$filename <- x 
  return(subject_old)
}


read_fread_old <- function(x, varname){
  fread_old_statements <- data.frame(varname = c("leftresponses", "rightresponses", "rewards"),
                                 statement = c("awk '/^BinsInActiveResponses/{flag=1;next}/endl/{flag=0}flag' ",
                                               "awk '/^ResponsesActBins/{flag=1;next}/endl/{flag=0}flag' ",
                                               "awk '/BinRewards/{flag=1;next}/endl/{flag=0}flag' "))  #### 	 In=L Act=R  Rew=W InTS=U	ActTS=Y  RewTS=V  RewIRI=Z 	
  statement <- fread_statements[which(fread_statements$varname == varname),]$statement
  rawdata <- fread(paste0(statement, "'", x, "'"), fill = T, header = F)
  rawdata$filename <- x
  return(rawdata)
}

## OLD 
## read_fread_old
# function(x){
#   # binrewards <- fread(paste0("awk '/BinRewards/{flag=1;next}/ResponsesActBins/{flag=0}flag' ", "'", x, "'", " | grep -v \"[list|endl]\""))
#   binrewards <- fread(paste0("awk '/BinRewards/{flag=1;next}/ResponsesActBins/{flag=0}flag' ", "'", x, "'", "| grep -v \"endl\""), header = F)
#   binrewards$filename <- x
#   return(binrewards)
# }

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
  



# join_wfu_oli_cocaine <- function(x){
#   
# }
# olivier_cocaine_files <- grep(grep(list.files(path = ".", recursive = T, full.names = T), pattern = ".*txt", inv = T, value = T), pattern = ".*C01..*LGA", value = T) # filter to the new files for lga in cohort 1

## SECTION OFF R SCRIPT TO DIFFERENTIATE THESE FILES, XX ALSO SECTION OFF BASED ON PR AND FR (DON'T FORGET THE CONSTRAINTS ON THE PR)

################################
########## SHA #################
################################

########### NEW SHA
olivier_cocaine_files_sha <- grep(grep(list.files(path = ".", recursive = T, full.names = T), pattern = ".*txt", inv = T, value = T), pattern = ".*SHA", value = T) # 178 files
names_sha <- lapply(olivier_cocaine_files_sha, readsubjects) %>% rbindlist()

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

# use for actual names vector
names_sha_append <- names_sha %>% 
  rename("labanimalid"="V1") %>% 
  mutate(labanimalid = paste0(str_extract(toupper(labanimalid), "[MF]\\d{1,3}"), "_", str_extract(filename, "C\\d+"), "_", sub('.*HS', '', toupper(filename)), "_", sub(".*/.*/.*/", '', filename)))

# merge to this dataset to acquire more metadata
sha_boxes <- lapply(olivier_cocaine_files_sha, readboxes) %>% 
  rbindlist() %>% 
  rename("box" = "V1")

rewards_sha <- lapply(olivier_cocaine_files_sha, read_fread, "rewards") %>% unlist(recursive = F)
# names(rightresponses_sha) <- names_sha_append
names(rewards_sha) <- names_sha_append$labanimalid
rewards_sha_df <- rewards_sha %>% 
  rbindlist(fill = T, idcol = "labanimalid") %>% 
  separate(labanimalid, c("labanimalid", "file_cohort", "file_exp", "filename"), "_")

# double check all labanimalids are valid
rewards_sha_df[str_detect(rewards_sha_df$labanimalid, "^[MF]\\d+$", negate = T),] %>% View() # NA labanimalid
rewards_sha_df[str_detect(rewards_sha_df$labanimalid, "^[MF]\\d+$", negate = T),] %>% dplyr::filter(bin == "total") %>% View()

## merge(WFU_Olivier_co_test_df[, c("cohort", labanimalnumber", "rfid")])

########### OLD SHA

# get the files of interest
setwd("~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS/")
sha_old_files <- grep(list.files(path = ".", recursive = T, full.names = T), pattern = ".*Old.*SHA", value = T)

# extract the subject id information
sha_subject_old <- lapply(sha_old_files, read_subject_old) %>% rbindlist()
sha_subject_old_use <- sha_subject_old %>% 
  rename("labanimalid" = "V1") %>% 
  mutate(labanimalid = paste0(str_match(toupper(labanimalid), "[FM]\\d{1,3}"), "_", 
                              str_extract(filename, "C\\d+"), "_", 
                              sub("-.*", "", sub(".*HS([^.]+)[-].*", "\\1", toupper(filename))), "_", 
                              sub("C.*", "", sub(".*/.*/.*/.*/", "", filename)), "_", 
                              str_extract(filename, "\\d{8}(-\\d+)?"))) # subject id, cohort, experiment, computer, date
# # handle the problematic cases 
sha_problem_files <- c("./C01/Old/SHA/K2C01HSSHA06-20170807.txt", "./C02/Old/SHA/Q3C02HSSHA01-20171017.txt")
sha_problem_subjects <- lapply(sha_problem_files, function(x){
  subject_old <- fread(paste0("grep -iA1 \"ratnumber\" ", "'", x, "'"), header = F)
  subject_old$filename <- x 
  return(subject_old)
}) # to find the order that the problematic id falls w/in the files

sha_problem_subjects_df <- sha_subject_old_use[which(sha_subject_old_use$filename %in% sha_problem_files),] %>% 
  tibble::rownames_to_column() %>% 
  group_by(filename) %>% 
  do(head(.,1)) %>% # happens to be the first case for both; used sha_problem_subjects to find this
  ungroup() %>% 
  select(rowname, filename) %>% 
  mutate(rowname = as.numeric(rowname))
sha_subject_old_use_test <- do.call("rbind", list(sha_subject_old_use[1:sha_problem_subjects_df$rowname[1]-1, ], 
                                                  data.frame(labanimalid="NA_C01_SHA06_K2_20170807", filename="./C01/Old/SHA/K2C01HSSHA06-20170807.txt"), 
                                                  sha_subject_old_use[sha_problem_subjects_df$rowname[1]+1:sha_problem_subjects_df$rowname[2]-1, ], 
                                                  data.frame(labanimalid="NA_C02_SHA01_K1_20170117", filename="./C02/Old/SHA/Q3C02HSSHA01-20171017.txt"),
                                                  sha_subject_old_use[sha_problem_subjects_df$rowname[2]+1:nrow(sha_subject_old_use), ])) %>% 
  dplyr::filter(!is.na(labanimalid))
rownames(sha_subject_old_use_test) <- seq(length=nrow(sha_subject_old_use_test)) # fix missing id cases


# approach 2 for handling problematic cases
sha_subjects <- lapply(sha_old_files, function(x){
  subject_old <- fread(paste0("grep -iA1 \"ratnumber\" ", "'", x, "'"), header = F)
  subject_old$filename <- x 
  return(subject_old)
}) # to find the order that the problematic id falls w/in the files

sha_subjects_old_use <- sha_subjects %>% rbindlist() %>% 
  rename("labanimalid" = "V1") %>% 
  mutate(labanimalid=replace(labanimalid, labanimalid=="999", "F000"), # create placeholder for the problematic cases
    labanimalid = paste0(str_match(toupper(labanimalid), "[FM]\\d{1,3}"), "_", 
                              str_extract(filename, "C\\d+"), "_", 
                              sub("-.*", "", sub(".*HS([^.]+)[-].*", "\\1", toupper(filename))), "_", 
                              sub("C.*", "", sub(".*/.*/.*/.*/", "", filename)), "_", 
                              str_extract(filename, "\\d{8}(-\\d+)?"))) %>%  # subject id, cohort, experiment, computer, date  
  dplyr::filter(!grepl("^NA", labanimalid))



# extract the subject sha data
sha_old <- lapply(sha_old_files, read_fread_old) 
sha_old <- lapply(sha_old, function(x){
  Index <- which(x[,1]=="list")
  if(x[(Index+1),] == 12){
    x <- x[-(Index+1),]
  }
  return(x)
}) %>% rbindlist() # use indexing to remove the 12 value if it follows list

list_indices <- grep("list", sha_old$V1)
sha_old_split <- split(sha_old, cumsum(1:nrow(sha_old) %in% list_indices))
names(sha_old_split) <- sha_subject_old_use$labanimalid ## XX NOW TWO SUBJECT LNES SHORTER THAN THE DATA # but this code shows that there should be 1255 'ratnumber' values so reconsider how to extract
# find . -regex ".*/C0[0-9]/Old/SHA/[^/]*.txt" -exec grep -i 'ratnumber' {} \; | wc -l
# find . -regextype sed -regex ".*/C0[0-9]/Old/SHA" -exec sh -c 'grep -ira1 "ratnumber" | grep -Eiv "[F|M]" | grep -vi "ratnumber"' sh {} \;
# find . -regex ".*/C0[0-9]/Old/SHA/[^/]*.txt" -exec sh -c 'grep -ira1 "ratnumber"' sh {} \;

sha_old_df <- sha_old_split %>% 
  rbindlist(idcol = "labanimalid") %>%
  rename("counts" = "V1") %>% 
  dplyr::filter(counts != "list") %>% 
  separate(labanimalid, into = c("labanimalid", "file_cohort", "file_exp", "computer", "file_date"), sep = "_")

# extract iri data 
iri <- lapply(sha_old_files, read_iri_old) %>% rbindlist()
iri %>%
  dplyr::filter(iritime == "list", iricode != "list") %>% 
  nrow() # should be none! since we are using cbind, make sure that the df's are "synced"
iri_indices <- grep("list", iri$iritime)
sha_old_iri <- split(iri, cumsum(1:nrow(iri) %in% iri_indices))
names(sha_old_iri) <- sha_subject_old_use$labanimalid ## still wrong because of the data are longer than the number of names we can assign

cohort1_old_iri_df <- lapply(cohort1_old_iri, convert_iri_matrix_to_df) %>% rbindlist(idcol = "labanimalid")


################################
########## LGA #################
################################
olivier_cocaine_files_lga <- grep(grep(list.files(path = ".", recursive = T, full.names = T), pattern = ".*txt", inv = T, value = T), pattern = ".*LGA", value = T) # 329 files
names_lga <- lapply(olivier_cocaine_files_lga, readsubjects) %>% rbindlist()
names_append <- names %>% 
  select(V1) %>% 
  unlist() %>% 
  as.vector() %>% 
  paste0(gsub(".*(C\\d+)HS(.*)","\\1\\2", names$filename)) %>% 
  toupper() %>%
  str_extract("F\\d+.*")
names_append <- names_append[!is.na(names_append)]


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
