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
subjects <- fread(paste0("awk '/Subject/{print $2}' ", "'", x, "'"), fill = T)
return(subjects)
}


names(rewardstest_split) <- lapply(olivier_cocaine_files[2], readsubjects) %>% rbindlist() %>% unlist() %>% as.character()