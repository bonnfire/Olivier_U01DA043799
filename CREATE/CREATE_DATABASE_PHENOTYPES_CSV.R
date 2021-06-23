## database objects

# extracted from Olivier_RAW_AND_EXCEL_QC.R
# CREATE A SELF ADMIN PHENOTYPES TABLES

## SHA
c01_sha <- cocaine_intermediates_xl_sha_corrected %>%
  subset(cohort == "C01") %>% 
  mutate(exp = tolower(exp)) %>%
  spread(exp, rewards) %>%
  mutate(mean_sha_last3 = rowMeans(select(., starts_with("sha")), na.rm = T))


## PR 


## for 06/29 meeting 

# use Excel values and then fill in blanks with raw

# sha (3 traits)

sha_rawgwas <- read.csv("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/U01/Olivier_George_U01DA043799 (Cocaine)/excel_and_csv_files/cocaine_sha_raw_c01_11_oldnewdirs.csv", stringsAsFactors = F) %>%  # join to raw file 
  mutate_at(vars(one_of("box", "session_duration")), as.character) %>% pivot_longer(cols = where(is.numeric), values_to = "raw") %>% 
  mutate(subject = ifelse(subject == "M7678", "M768", subject)) %>% 
  select(subject, box, exp, name, raw, cohort) %>% 
  left_join(cocaine_metadata_df %>% select(labanimalid, rfid), by = c("subject" = "labanimalid")) %>% # validate the rfid 
  rowwise() %>% 
  mutate(subject = ifelse(is.na(rfid), gsub("F", "M", subject), subject)) %>% # fix the sex in the ID and then rejoin to metadata
  ungroup() %>% 
  select(-rfid) %>% 
  left_join(cocaine_metadata_df %>% select(labanimalid, rfid), by = c("subject" = "labanimalid")) # validate the rfid 
  
## potentially needs to be corrected by the compromised table, bc they're naives?
# use to QC sha_rawgwas %>% subset(is.na(rfid)) %>% distinct(subject, cohort) %>% mutate(subject = parse_number(subject)) %>% select(subject) %>% unlist() %>% paste0(collapse = "|") %>% cat 

sha_gwasprelim2 <- oliviercocaine_excel_all %>% 
  select(cohort, measurement, labanimalid, rfid, matches("sha")) %>% 
  pivot_longer(cols = matches("sha")) %>% 
  mutate(measurement = ifelse(grepl("date", name), "date", measurement),
         name = ifelse(grepl("date", name), gsub("date_", "", name), name)) %>% 
  distinct() %>% 
  full_join(sha_rawgwas, by = c("labanimalid" = "subject", "measurement" = "name", "name" = "exp")) %>% # join to raw file
  mutate(cohort = coalesce(cohort.x, cohort.y)) %>% 
  select(-cohort.x, -cohort.y) %>% 
  mutate(value = coalesce(value, as.character(raw))) %>% 
  full_join(openxlsx::read.xlsx("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/U01/Olivier_George_U01DA043799 (Cocaine)/excel_and_csv_files/cocaine_sha_to_presses_final.xlsx") %>%  # join to TO file 
              clean_names %>% 
              mutate_at(vars(starts_with("sha")), as.numeric) %>% 
              pivot_longer(cols = where(is.numeric), values_to = "timeout") %>% 
              # mutate(subject = ifelse(subject == "M7678", "M768", subject)) %>% 
              select(labanimalid, sex, name, timeout), by = c("labanimalid", "name")) %>% 
  group_by(labanimalid) %>%
  # fill(cohort, .direction = "downup") %>% 
  fill(rfid, .direction = "downup") %>%
  ungroup() %>%
  subset(!is.na(value)) %>% 
  group_by(labanimalid, measurement) %>% do(tail(., n=3)) %>% # extract the last three sessions
  ungroup() %>% 
  left_join(cocaine_metadata_df %>% select(labanimalid, rfid), by = "labanimalid") %>%  # fill in rfid information
  mutate(rfid = coalesce(rfid.x, rfid.y))


# pr 

pr_rawgwas <- read.csv("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/U01/Olivier_George_U01DA043799 (Cocaine)/excel_and_csv_files/cocaine_pr_raw_c01_11_oldnewdirs.csv", stringsAsFactors = F) %>%  # join to raw file 
  mutate_at(vars(one_of("box")), as.character) %>% pivot_longer(cols = where(is.numeric), values_to = "raw") %>% 
  left_join(cocaine_metadata_df %>% select(labanimalid, rfid), by = c("subject" = "labanimalid")) # validate the rfid 

# fix wrong subjects
pr_rawgwas %>% subset(is.na(rfid))

# use this for all combos to find the animals pr_rawgwas %>% subset(cohort == "C09"&box == "16") %>% distinct(subject)
pr_rawgwas <- pr_rawgwas %>% 
  mutate(subject = ifelse(subject == "F951", "M951", subject),
         subject = ifelse(subject == "F952", "M952", subject),
         subject = ifelse(subject == "M3", "F3", subject)) %>% 
  select(-rfid) %>% # remove and rejoin to the metadata
  left_join(cocaine_metadata_df %>% select(labanimalid, rfid), by = c("subject" = "labanimalid")) # validate the rfid 

pr_rawgwas %>% subset(is.na(rfid))

pr_rawgwas <- oliviercocaine_excel_all %>%  # fill in rfid information
  select(cohort, measurement, labanimalid, rfid, matches("pr")) %>% 
  pivot_longer(cols = matches("pr")) %>% 
  mutate(measurement = ifelse(grepl("date", name), "date", measurement),
         name = ifelse(grepl("date", name), gsub("date_", "", name), name)) %>% 
  distinct() %>% 
  full_join(pr_rawgwas, by = c("labanimalid" = "subject", "measurement" = "name", "name" = "exp")) %>% 
  mutate(cohort = coalesce(cohort.x, cohort.y),
         rfid = coalesce(rfid.x, rfid.y),
         value = coalesce(value, as.character(raw))) %>% 
    select(-matches("[.]")) %>% 
  distinct() 

# calculate the pr related values 
pr_rawgwas_traits <- pr_rawgwas %>% 
  subset(measurement == "rewards") %>% distinct(cohort, rfid, labanimalid, box, room, name, value) %>% group_by(labanimalid) %>% 
  fill(c("box", "room"), .direction = "downup") %>% pivot_wider(names_from = name, values_from = value) %>% ungroup() %>% 
  select(-preshock,
         -pr04) %>% # only 56 animals in C07  
  left_join(pr_rawgwas %>% 
              subset(measurement == "rewards"&name %in% c("pr02", "pr03")) %>% 
              group_by(labanimalid) %>% 
              mutate(pr_max_02_03 = max(value)) %>% distinct(labanimalid, pr_max_02_03),
            by = "labanimalid")

# shock

shock_rawgwas <- read.csv("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/U01/Olivier_George_U01DA043799 (Cocaine)/excel_and_csv_files/cocaine_shock_raw_c01_11_oldnewdirs.csv", stringsAsFactors = F) %>%  # join to raw file 
  mutate_at(vars(one_of("box")), as.character) %>% pivot_longer(cols = where(is.numeric), values_to = "raw") %>% 
  left_join(cocaine_metadata_df %>% select(labanimalid, rfid), by = c("subject" = "labanimalid")) # validate the rfid 

# fix wrong subjects
shock_rawgwas %>% subset(is.na(rfid))

# using this to clean up 
shock_rawgwas <- shock_rawgwas %>% 
  rowwise() %>% 
  mutate(subject = ifelse(is.na(rfid), gsub("F", "M", subject), subject)) %>% # checked the Excel to verify
  ungroup() %>% 
  select(-rfid) %>% # remove and rejoin to the metadata
  left_join(cocaine_metadata_df %>% select(labanimalid, rfid), by = c("subject" = "labanimalid")) # validate the rfid 

# calculate the shock related values 
shock_rawgwas_traits <- shock_rawgwas %>% 
  subset(exp %in% c("preshock", "shock03")&name == "rewards_w_postshock1") %>%
  distinct(cohort, rfid, subject, box, exp, room, raw) %>% 
  rename("labanimalid" = "subject") %>% 
  group_by(labanimalid) %>% 
  fill(box, room) %>% 
  ungroup() %>% 
  distinct() %>% 
  pivot_wider(names_from = exp, values_from = raw) %>% 
  mutate(shock_03_pre = (shock03-preshock)/preshock)







# bind all objects into one object and create csv file 
gwas_prelim2 <- full_join()
  