## upload data onto wfu

library(RPostgreSQL)

# set up the connection
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, user='postgres', password='postgres', dbname='U01')


# drop table if it already exists
# if (dbExistsTable(con, c("u01_tom_jhou","jhou_wfu_master")))
#   dbRemoveTable(con, c("u01_tom_jhou","jhou_wfu_master"))

# write data frame to database

Olivier_Cocaine_df_sql$cohort %>% table()
test1 <- Olivier_Cocaine_df_sql %>% subset(cohort != "C11") 
test2 <- Olivier_Cocaine_df_sql %>% subset(cohort == "C11")
dbWriteTable(con, c("u01_olivier_george_cocaine","olivier_rewards"), value = test1, row.names = FALSE)
append_sql <- sqlAppendTable(con, c("u01_olivier_george_cocaine","olivier_rewards"), test2, row.names = F) 
dbExecute(con, append_sql)
dbExecute(con, "INSERT INTO u01_olivier_george_cocaine.olivier_rewards (cohort, rfid, labanimalid, exp, rewards) VALUES ('C11', '933000320187309', NULL, 'SHA08', NULL)")

upd_sql <- "UPDATE olivier_rewards
            SET t.Val = s.Val
            FROM dbo.tbl t 
            INNER JOIN dbo.tbl_staging s 
               ON t.[Key] = s.[Key]"

dbSendQuery(con, upd_sql)

sql <- "
    select rewards
    from u01_olivier_george_cocaine.olivier_rewards
    where cohort='C01'
"
results <- dbGetQuery(con, sql)

