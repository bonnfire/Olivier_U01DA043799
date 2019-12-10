## extract raw 
setwd("~/Dropbox (Palmer Lab)/GWAS (1)/Cocaine/Cocaine GWAS")

# after cohort 5, there are only new files

# extract the new files 

olivier_cocaine_files <- grep(grep(list.files(path = ".", recursive = T, full.names = T), pattern = ".*txt", inv = T, value = T), pattern = ".*C01..*LGA", value = T)
# filter to the new files for lga in cohort 1

readrewards <- function(x){
  rewards <- fread(paste0("awk '/W:/{flag=1;next}/Y:/{flag=0}flag' ", "'", x, "'"), fill = T)
  return(rewards)
}
rewardstest <- lapply(olivier_cocaine_files[2], readrewards)
indices <- grep("^0:$", rewardstest[[1]]$V1)
rewardstest_split <- split(rewardstest[[1]], cumsum(1:nrow(rewardstest[[1]]) %in% indices))


readsubjects <- function(x){
subjects <- fread(paste0("awk '/Subject/{print $2}' ", "'", x, "'"),fill = T,header=F)
return(subjects)
}


names(rewardstest_split) <- lapply(olivier_cocaine_files[2], readsubjects) %>% rbindlist() %>% unlist() %>% as.character()
lapply(olivier_cocaine_files[2], readsubjects)


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
gdata::unmatrix(rightresponses_time_split[[1]])

test <- rightresponses_time_split[[1]][,-1]
test <- test[rowSums(test[, -1] > 0) != 0, ]
data.frame(timestamps = as.vector(data.matrix(test))) %>% head