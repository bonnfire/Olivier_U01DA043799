### demo 
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
