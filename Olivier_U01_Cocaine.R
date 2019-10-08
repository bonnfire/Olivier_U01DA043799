setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/Olivier_George_U01DA043799 (Cocaine)")

Olivier_U01_Cocaine_original <- read_excel()

# QC: sex vs shipment data sex  
## process WFU data 
WFU_Olivier_co_sexQC <- WFU_Olivier_co
WFU_Olivier_co_sexQC[[1]] <- data.frame(rfid = WFU_Olivier_co[[1]]$`Transponder ID` %>% as.character(),
                                        sex = WFU_Olivier_co[[1]]$`Sex`)
WFU_Olivier_co_sexQC[[2]] <- data.frame(rfid = WFU_Olivier_co[[2]]$`Transponder ID` %>% as.character(),
                                        sex = WFU_Olivier_co[[2]]$`Sex`)
WFU_Olivier_co_sexQC[[3]] <- data.frame(rfid = WFU_Olivier_co[[3]]$`Transponder ID` %>% as.character(),
                                        sex = WFU_Olivier_co[[3]]$`Sex`)
WFU_Olivier_co_sexQC[[4]] <- data.frame(rfid = WFU_Olivier_co[[4]]$`Transponder ID` %>% as.character(),
                                        sex = WFU_Olivier_co[[4]]$`Sex`)
WFU_Olivier_co_sexQC[[5]] <- data.frame(rfid = WFU_Olivier_co[[5]][6],
                                        sex = WFU_Olivier_co[[5]][5])
colnames(WFU_Olivier_co_sexQC[[5]]) <- c("rfid", "sex")
WFU_Olivier_co_sexQC[[5]] <- WFU_Olivier_co_sexQC[[5]][grep("^9", WFU_Olivier_co_sexQC[[5]]$rfid), ]
WFU_Olivier_co_sexQC[[6]] <- data.frame(rfid = WFU_Olivier_co[[6]][6],
                                        sex = WFU_Olivier_co[[6]][5])
colnames(WFU_Olivier_co_sexQC[[6]]) <- c("rfid", "sex")
WFU_Olivier_co_sexQC[[6]] <- WFU_Olivier_co_sexQC[[6]][grep("^9", WFU_Olivier_co_sexQC[[6]]$rfid), ]
WFU_Olivier_co_sexQC[[7]] <- data.frame(rfid = WFU_Olivier_co[[7]][6],
                                        sex = WFU_Olivier_co[[7]][5])
colnames(WFU_Olivier_co_sexQC[[7]]) <- c("rfid", "sex")
WFU_Olivier_co_sexQC[[7]] <- WFU_Olivier_co_sexQC[[7]][grep("^9", WFU_Olivier_co_sexQC[[7]]$rfid), ]
WFU_Olivier_co_sexQC[[8]] <- data.frame(rfid = WFU_Olivier_co[[8]][6],
                                        sex = WFU_Olivier_co[[8]][5])
colnames(WFU_Olivier_co_sexQC[[8]]) <- c("rfid", "sex")
WFU_Olivier_co_sexQC[[8]] <- WFU_Olivier_co_sexQC[[8]][grep("^9", WFU_Olivier_co_sexQC[[8]]$rfid), ]
WFU_Olivier_co_sexQC[[9]] <- data.frame(rfid = WFU_Olivier_co[[9]][6],
                                        sex = WFU_Olivier_co[[9]][5])
colnames(WFU_Olivier_co_sexQC[[9]]) <- c("rfid", "sex")
WFU_Olivier_co_sexQC[[9]] <- WFU_Olivier_co_sexQC[[9]][grep("^9", WFU_Olivier_co_sexQC[[9]]$rfid), ]

WFU_Olivier_co_sexQC <- bind_rows(WFU_Olivier_co_sexQC, .id = "shipmentdate")

## get this function to work
WFU_Olivier_co_sexQC[[7]] <- WFU_Olivier_co_sexQC[[7]][grep("^9", WFU_Olivier_co_sexQC[[7]]$rfid), ]
WFU_Olivier_co_sexQC[5:9] <- lapply(WFU_Olivier_co[5:9], function(x) x[-1])

colnames(WFU_Olivier_co_sexQC[5:9]) <- lapply(WFU_Olivier_co[5:9], function(x) x[1])

## process Olivier data 
library(dplyr)
library(data.table)
library(tidyverse)
library(readxl)
setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/Olivier_George_U01DA043799 (Cocaine)/Olivier_George_U01/DATA Updated")
path <- getwd()

mylist.names <- rep(paste0("cohort", 1:8), each = 2) 
mylist <- sapply(mylist.names, function(x) NULL)

# u01.importxlsx.olivier <- function(path){
#   xlname <- list.files(path = path, pattern = "*.xlsx", full.names = F)
#   lapply(seq_along(mylist, function(i){
#     mylist[[i]] <- lapply(xlname, function(x){
#       read_excel(path = x)
#       path_sheetnames <- excel_sheets(x)
#     })
#   }))
# }
# Olivier_co_all <- u01.importxlsx.olivier(path)

u01.importxlsx <- function(xlname){
  df <- lapply(excel_sheets(path = xlname), read_excel, path = xlname, col_names = FALSE)
  names(df) <- excel_sheets(xlname)
  return(df)
}  # diff from u01 function because of col_names and need to transpose for observations each row


xlname <- list.files(path = path, pattern = "*.xlsx", full.names = F)
for(i in 1:length(mylist)){
  for(j in 1:length(xlname)){
    data <- u01.importxlsx(xlname[j]) %>% 
      as.data.table
  }
  mylist[[i]] <- data[[1]]
  mylist[[i + 1]] <- data[[2]]
  names(mylist) <- gsub(" ", "", paste0(names(mylist), sapply(xlname, excel_sheets)))
}

# transpose data tables
dcast(melt(mydata, id.vars = "col0"), variable ~ col0)

sapply(xlname, excel_sheets)

# testing the setup for data tables
empty_list <- list()
dt <- u01.importxlsx("C01_cocaine.xlsx")[[1]] %>% 
  drop_na() %>%
  as.data.table()
empty_list <- append(empty_list, list(dt))
dt2 <- u01.importxlsx("C01_cocaine.xlsx")[[2]] %>% as.data.table()
empty_list <- append(empty_list, list(dt2))
empty_list[2]

# separate approach
options(scipen = 100) # fixes scientific notation
test <- u01.importxlsx("C01_cocaine.xlsx")[[1]] %>%
  as.data.table %>% 
  na.omit(cols = seq_along('..1')) 

# test <- u01.importxlsx("C01_cocaine.xlsx")[[1]] %>%
#   as.data.table %>% 
#   na.omit(cols = seq_along(cols)) 
# cols <- grep("^C0|Rat|Description", names(test), value = TRUE) # doesn't apply once I remove the columns
# test_subset <- test[, ..cols]

test2 <- dcast(melt(test, id.vars = "..1"), variable ~ "..1")


# todo: 
# change the chr to num


# why are there 1,032,228 observations 


test <- u01.importxlsx("testtodt.xlsx") %>%
  as.data.table %>% 
  na.omit(cols = seq_along('Sheet1....1')) 
test2 <- test %>% # take the data.table
  summarise_all(list( ~ sum)) %>% # get the sum of each column
  gather(variable, sum) 
