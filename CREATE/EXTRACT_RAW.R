## extract raw 
setwd("~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS")

# after cohort 5, there are only new files

# extract the new files 

olivier_cocaine_files <- grep(grep(list.files(path = ".", recursive = T, full.names = T), pattern = ".*txt", inv = T, value = T), pattern = ".*C01..*LGA", value = T)
# filter to the new files for lga in cohort 1

#extract names to be assigned for various tables later
readsubjects <- function(x){
  subjects <- fread(paste0("awk '/Subject/{print $2}' ", "'", x, "'"),fill = T,header=F)
  subjects$filename <- x
  return(subjects)
}
names <- lapply(olivier_cocaine_files, readsubjects) %>% rbindlist()
names_append <- names$V1 %>% paste0(gsub(".*(C\\d+)HS(.*)","\\1\\2", names$filename))

#extract the rewards 
readrewards <- function(x){
  rewards <- fread(paste0("awk '/W:/{flag=1;next}/Y:/{flag=0}flag' ", "'", x, "'"), fill = T)
  indices_reward <- grep("^0:$", rewards$V1)
  split_rewards <- split(rewards, cumsum(1:nrow(rewards) %in% indices_rewards))
  return(split_rewards)
} 
rewards_cohort1 <- lapply(olivier_cocaine_files, readrewards) %>% unlist(recursive = F) # needs one level of unlisting bc otherwise it is a list of lists
names(rewards_cohort1) <- names_append

#extract timestamps for rewards 
readreward_times <- function(x){
  rewards_times <- fread(paste0("awk '/V:/{flag=1;next}/W:/{flag=0}flag' ", "'", x, "'"), fill = T)
  indices_rewards_times <- grep("^0:$", rewards_times$V1)
  split_rewards_times <- split(rewards_times, cumsum(1:nrow(rewards_times) %in% indices_rewards_times))
  return(split_rewards_times)
}
rewards_times_cohort1 <- lapply(olivier_cocaine_files, readreward_times) %>% unlist(recursive = F) 
names(rewards_times_cohort1) <- names_append

#extract timestamps for right responses
readrightresponses_time <- function(x){
  rightresponses <- fread(paste0("awk '/Y:/{flag=1;next}/^$/{flag=0}flag' ", "'", x, "'"), fill = T) # matching to empty lines
  indices_right_resp_time <- grep("^0:$", rightresponses$V1)
  split_right_resp_time <- split(rightresponses, cumsum(1:nrow(rightresponses) %in% indices_right_resp_time))
  return(split_right_resp_time)
}
rightresponses_time <- lapply(olivier_cocaine_files, readrightresponses_time) %>% unlist(recursive = F) 
names(rightresponses_time) <- names_append

#create the dataframe with vector in it
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
    processeddata_df <- data.frame(timestamps = as.vector(t(data.matrix(indexremoved)))) # transpose to get by row
    return(processeddata_df)
    })
  }
  else{
    processeddata <- lapply(split_data, function(x){
      indexremoved <- x[,-1]
      nonzerorows <- indexremoved[rowSums(indexremoved) > 0, ]
      processeddata_df <- data.frame(timestamps = as.vector(t(data.matrix(nonzerorows)))) # transpose to get by row
      return(processeddata_df)
      })
  }

  names(processeddata) <- grep("C01LGA01", names_append, value = T)

  return(processeddata)
  }

rightresponseslga01 <- read_fread(olivier_cocaine_files[[2]], "rightresponses")

test <- rightresponses_time_split[[1]][,-1]
test <- test[rowSums(test[, -1] > 0) != 0, ]
data.frame(timestamps = as.vector(matrix(test, byrow = T))) 




