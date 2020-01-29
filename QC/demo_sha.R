### demo new
SUBJECTS <- process_subjects_new(sha_new_files[10:15]) ## sha_subjects_new %>% dplyr::filter(!grepl("^[MF]", labanimalid)) %>% head # will need to use the comp df to name sessions
REWARDS <- lapply(sha_new_files[10:15], read_fread_new, "rewards") %>% unlist(recursive = F)
nrow(SUBJECTS) == length(REWARDS)
names(REWARDS) <- SUBJECTS$labanimalid
REWARDS_DF <- REWARDS %>% rbindlist(idcol = "labanimalid") %>% dplyr::filter(bin == "total") %>% separate(labanimalid, into = c("labanimalid", "cohort", "exp", "filename", "date", "time"), sep = "_") %>% mutate(date = lubridate::mdy(date), time = chron::chron(times = time))
  # mutate_at(vars("date"), lubridate::mdy) %>% mutate_at(vars("time"), chron::chron)
VALID <- date_time_subject_df_comp %>% select(labanimalid, cohort, exp, filename, valid, start_date, start_time) %>% rename("date" = "start_date", "time" = "start_time")
REWARDS_MERGE <- REWARDS_DF %>% left_join(., VALID, by = c("cohort", "exp", "filename", "date", "time"))

setDT(REWARDS_DF)
setDT(date_time_subject_df)
date_time_subject_df[, c("box","start_time", "end_time", "directory", "experiment_duration"):=NULL]
date_time_subject_df[REWARDS_DF, on = c('cohort','exp', 'start_date')]



### demo for lga new
SUBJECTS <- process_subjects_new(lga_new_files[190]) # grep("BSB273BC08HSLGA01", lga_new_files) for F801
read_rewards_new <- function(x){
  rewards <- fread(paste0("awk '/W:/{flag=1;next}/5:/{flag=0}flag' ", "'", x, "' | awk '/0:/{print $2}'"))
  return(rewards)
}
REWARDS_DF <- lapply(lga_new_files[190], read_rewards_new) %>% rbindlist() %>% bind_cols(SUBJECTS) %>% 
  separate(labanimalid, into = c("labanimalid", "cohort", "exp", "filename", "date", "time"), sep = "_") %>% 
  mutate(date = lubridate::mdy(date), time = chron::chron(times = time)) %>%  
  left_join(., date_time_subject_df_comp %>% 
              select(cohort, exp, filename, valid, start_date, start_time) %>% 
              rename("date" = "start_date", "time" = "start_time"), 
            by = c("cohort", "exp", "filename", "date", "time")) %>% 
  dplyr::filter(valid == "yes") %>% 
  rename("rewards" = "V1") %>% 
  mutate(time = as.character(time))
