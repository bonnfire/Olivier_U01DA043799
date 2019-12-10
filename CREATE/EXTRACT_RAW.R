## extract raw 
setwd("~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS")

# after cohort 5, there are only new files

# extract the new files 

olivier_cocaine_files <- grep(grep(list.files(path = ".", recursive = T, full.names = T), pattern = ".*txt", inv = T, value = T), pattern = ".*C01..*LGA", value = T)
# filter to the new files for lga in cohort 1

#extract and unlist the rewards 
readrewards <- function(x){
  rewards <- fread(paste0("awk '/W:/{flag=1;next}/Y:/{flag=0}flag' ", "'", x, "'"), fill = T)
  return(rewards)
}
rewards_cohort1 <- lapply(olivier_cocaine_files, readrewards)
rewards_cohort1_split <- lapply(rewards_cohort1, function(x){
  indices_rewards <- grep("^0:$", x$V1)
  split_x <- split(x, cumsum(1:nrow(x) %in% indices_rewards))
  return(split_x)
}) %>% unlist(recursive = F) 

#extract and assign subject names to rewards
readsubjects <- function(x){
subjects <- fread(paste0("awk '/Subject/{print $2}' ", "'", x, "'"),fill = T,header=F)
subjects$filename <- x
return(subjects)
}
names <- lapply(olivier_cocaine_files, readsubjects) %>% rbindlist()
names_append <- names$V1 %>% paste0(gsub(".*(C\\d+)HS(.*)","\\1\\2", names$filename))

names(rewards_cohort1_split) <- names_append

#extract timestamps for rewards 
readreward_times <- function(x){
  rewards_times <- fread(paste0("awk '/V:/{flag=1;next}/W:/{flag=0}flag' ", "'", x, "'"), fill = T)
  indices_rewards_times <- grep("^0:$", rewards_times$V1)
  split_rewards_times <- split(rewards_times, cumsum(1:nrow(rewards_times) %in% indices_rewards_times))
  return(split_rewards_times)
}

rewards_times_cohort1 <- lapply(olivier_cocaine_files, readreward_times) %>% unlist(recursive = F) 
names(rewards_times_cohort1) <- names_cohort1



readrightresponses <- function(x){
  subjects <- fread(paste0("awk '/Subject/{print $2}' ", "'", x, "'"),fill = T,header=F)
  return(subjects)
}

readrightresponses_time <- function(x){
  rightresponses <- fread(paste0("awk '/Y:/{flag=1;next}/^$/{flag=0}flag' ", "'", x, "'"), fill = T) # matching to empty lines
  return(rightresponses)
}
rightresponses_time <- lapply(olivier_cocaine_files[2], readrightresponses_time)
rightresponses_time_indices <- grep("^0:$", rightresponses_time[[1]]$V1)
rightresponses_time_split <- split(rightresponses_time[[1]], cumsum(1:nrow(rightresponses_time[[1]]) %in% rightresponses_time_indices))
names(rightresponses_time_split) <- lapply(olivier_cocaine_files[2], readsubjects) %>% rbindlist() %>% unlist() %>% as.character()

test <- rightresponses_time_split[[1]][,-1]
test <- test[rowSums(test[, -1] > 0) != 0, ]
data.frame(timestamps = as.vector(data.matrix(test))) 