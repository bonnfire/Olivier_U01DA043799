## upload data onto wfu

library(RPostgreSQL)

# set up the connection
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, user='postgres', password='postgres', dbname='U01')

## PICK UP AND REDO BC ALL OF THE REWARDS ARE NA


# write data frame to database
## to add dataframes (by cohort)
## write to a temporary table in the database with dbWriteTable
# then issue an SQL statement to insert into your table as you would within the db environment
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, user='postgres', password='postgres', dbname='U01')
update_sql_db <- function(df, df2){
  df_to_add <- anti_join(df2, df)
  
}



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
