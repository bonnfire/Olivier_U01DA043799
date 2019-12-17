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


join_wfu_oli_cocaine <- function(x){
  
}
# olivier_cocaine_files <- grep(grep(list.files(path = ".", recursive = T, full.names = T), pattern = ".*txt", inv = T, value = T), pattern = ".*C01..*LGA", value = T) # filter to the new files for lga in cohort 1

## SECTION OFF R SCRIPT TO DIFFERENTIATE THESE FILES, XX ALSO SECTION OFF BASED ON PR AND FR (DON'T FORGET THE CONSTRAINTS ON THE PR)

################################
########## SHA #################
################################
olivier_cocaine_files_sha <- grep(grep(list.files(path = ".", recursive = T, full.names = T), pattern = ".*txt", inv = T, value = T), pattern = ".*SHA", value = T) # 178 files
names_sha <- lapply(olivier_cocaine_files_sha, readsubjects) %>% rbindlist()

# 12/16 for rsm: 
names_sha_rsm <- names_sha %>% 
  rename("labanimalid"="V1") %>% 
  mutate(labanimalid = paste0(str_extract(toupper(labanimalid), "[MF]\\d{1,3}"), "_", str_extract(filename, "C\\d+")), 
         #file_exp = str_extract(toupper(filename), "\\D+\\d+$"))
  file_exp = sub('.*HS', '', toupper(filename))) %>% 
  group_by(labanimalid) %>% 
  add_count(file_exp) %>% 
  ungroup() 

exps_vs_id <- table(names_sha_rsm$file_exp, factor(names_sha_rsm$labanimalid, levels = unique(gtools::mixedsort(names_sha_rsm$labanimalid)))) %>% t()
exps_vs_id <- cbind(exps_vs_id, Total = rowSums(exps_vs_id))
# exps_vs_id <- rbind(exps_vs_id, Total = colSums(exps_vs_id))
openxlsx::write.xlsx(exps_vs_id, file = "exps_vs_id_2.xlsx",col.names=TRUE, row.names=TRUE)

#openxlsx::write.xlsx(exps_vs_id, file = "exps_vs_id.xlsx",col.names=TRUE, row.names=TRUE)

# names_sha_append <- names_sha %>% 
#   select(V1) %>% 
#   unlist() %>% 
#   as.vector() %>% 
#   paste0(gsub(".*(C\\d+)HS(.*)","\\1\\2", names_sha$filename)) %>% 
#   toupper() %>%
#   str_extract("[FM]\\d+.*")
# names_sha_append <- names_sha_append[!is.na(names_sha_append)] #with na, 2783; 2716 without na

names_sha_append <- names_sha %>% 
  rename("labanimalid"="V1") %>% 
  mutate(labanimalid = paste0(str_extract(toupper(labanimalid), "[MF]\\d{1,3}"), "_", str_extract(filename, "C\\d+"), "_", sub('.*HS', '', toupper(filename))))
                              

rewards_sha <- lapply(olivier_cocaine_files_sha, read_fread, "rewards") %>% unlist(recursive = F)
# names(rightresponses_sha) <- names_sha_append
names(rewards_sha) <- toupper(names_sha_append$labanimalid)
rewards_sha_df <- rewards_sha %>% 
  rbindlist(fill = T, idcol = "labanimalid") %>% 
  mutate(file_cohort = str_extract(labanimalid, "C\\d+"), 
         file_exp = str_extract(labanimalid, "SHA\\d+$"), 
         labanimalid = str_extract(labanimalid,"^\\D\\d+"))
  ## merge(WFU_Olivier_co_test_df[, c("cohort", labanimalnumber", "rfid")])
  
  
  
  
  
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
%>% 
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
