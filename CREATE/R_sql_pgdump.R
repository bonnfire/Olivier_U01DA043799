## upload data onto wfu

library(RPostgreSQL)

# set up the connection
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, user='postgres', password='postgres', dbname='U01')

## PRACTICE FOR CODE SESSION 06/10/2020
upload_practice <- Olivier_Cocaine_df_sql %>% distinct(rfid, exp, .keep_all = T) ## temporary
dbWriteTable(con, c("u01_olivier_george_cocaine","olivier_rewards"), value = upload_practice, row.names = FALSE)
dbExecute(con,"ALTER TABLE u01_olivier_george_cocaine.olivier_rewards ADD PRIMARY KEY(rfid,exp)")

# troubleshoot
dbExecute(con,"drop table if exists u01_olivier_george_cocaine.olivier_rewards")



# write data frame to database
## to add dataframes (by cohort)
## write to a temporary table in the database with dbWriteTable
# then issue an SQL statement to insert into your table as you would within the db environment
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, user='postgres', password='postgres', dbname='U01')
# update_sql_db <- function(df, df2){
#   df_to_add <- anti_join(df2, df)
#   
# }



test1 <- Olivier_Cocaine_df_sql %>% distinct() %>% subset(cohort != "C11") 
test2 <- Olivier_Cocaine_df_sql %>% distinct() %>% subset(cohort == "C11")
dbWriteTable(con, c("u01_olivier_george_cocaine","olivier_rewards"), value = test1, row.names = FALSE)
dbWriteTable(con, c("u01_olivier_george_cocaine","olivier_rewards_temp"), value = test2, row.names = FALSE)
dbExecute(con,"ALTER TABLE u01_olivier_george_cocaine.olivier_rewards ADD PRIMARY KEY(rfid,exp)")
dbExecute(con,"insert into u01_olivier_george_cocaine.olivier_rewards select * from u01_olivier_george_cocaine.olivier_rewards_temp")
dbExecute(con,"drop table if exists u01_olivier_george_cocaine.olivier_rewards_temp")

# for troubleshooting 
dbExecute(con,"drop table if exists u01_olivier_george_cocaine.olivier_rewards")

# disconnect
dbDisconnect(con)




## in terminal
cd /tmp
sudo su postgres
pg_dump -d U01 -t u01_olivier_george_cocaine.olivier_rewards > olivier_rewards.sql
exit
cp olivier_rewards.sql /home/bonnie/Dropbox\ \(Palmer\ Lab\)/PalmerLab_Datasets/u01_george_oliviercocaine/database
## operation not permitted 
# To move a file or directory, you need to have write permissions on both SOURCE and DESTINATION. Otherwise, you will receive a permission denied error.


## 
mkdir {C01,C02,C03,C04,C05,C06,C07,C08,C09,C10,C11}
