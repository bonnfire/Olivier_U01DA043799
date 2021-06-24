## database objects

# extracted from Olivier_RAW_AND_EXCEL_QC.R
# CREATE A SELF ADMIN PHENOTYPES TABLES

# ## SHA
# c01_sha <- cocaine_intermediates_xl_sha_corrected %>%
#   subset(cohort == "C01") %>% 
#   mutate(exp = tolower(exp)) %>%
#   spread(exp, rewards) %>%
#   mutate(mean_sha_last3 = rowMeans(select(., starts_with("sha")), na.rm = T))
# 
# 


### TIER ONE TRAITS
## for 06/29 meeting 

sha_rawgwas <- read.csv("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/U01/Olivier_George_U01DA043799 (Cocaine)/excel_and_csv_files/cocaine_sha_raw_c01_11_oldnewdirs.csv", stringsAsFactors = F) %>%  # join to raw file 
  mutate_at(vars(one_of("box", "session_duration")), as.character) %>% pivot_longer(cols = where(is.numeric), values_to = "raw") %>% 
  select(subject, box, exp, name, raw, cohort) %>% 
  left_join(cocaine_metadata_df %>% select(labanimalid, rfid), by = c("subject" = "labanimalid")) 

sha_rawgwas %>% subset(is.na(rfid)) # validate the rfid ## there are ~30 cases that need to be fixed, but plan to come back to this after deadline 

## focus on the ones that are used for calculations
# sessions 1, 08-10
sha_rawgwas_1_08_10 <- sha_rawgwas %>% 
  subset(exp %in% c("sha01", "sha08", "sha09", "sha10"))

sha_rawgwas_1_08_10 %>% subset(is.na(rfid))

sha_rawgwas_1_08_10 <- sha_rawgwas_1_08_10 %>% 
  mutate(subject = ifelse(subject == "M7678", "M768", subject),
         rfid = ifelse(subject == "M768", "933000320046611", rfid)) %>% 
  rowwise() %>%
  mutate(subject = ifelse(is.na(rfid), gsub("F", "M", subject), subject)) %>% # fix the sex in the ID and then rejoin to metadata
  ungroup() %>%
  group_by(subject) %>% 
  fill(rfid, .direction = 'downup') %>% 
  ungroup() %>% 
  mutate(subject = ifelse(is.na(rfid), gsub("M", "F", subject), subject))# if it didn't fix, then switch back
  
# now use boxes to clean up
sha_rawgwas_1_08_10 %>% subset(is.na(rfid))

# sha_rawgwas_1_08_10 %>% subset(box == "14"&cohort == "C04") # check for non dupes in sessions

sha_rawgwas_1_08_10 <- sha_rawgwas_1_08_10 %>% 
  mutate(subject = ifelse(subject == "F414", "F424", subject)
         # subject = ifelse(subject == "F124", "F524", subject), # technically an error but there is already an existing session for this animal
         # subject = ifelse(subject == "F123", "F523", subject), # technically an error but there is already an existing session for this animal
         # subject = ifelse(subject == "F125", "F525", subject), # technically an error but there is already an existing session for this animal
         ) %>% 
  subset(!(subject %in% c("F124", "F123", "F125")&cohort == "C05")) %>% 
  group_by(subject) %>% 
  fill(rfid, .direction = "downup") %>% 
  ungroup()


sha_gwasprelim2 <- oliviercocaine_excel_all %>% 
  select(cohort, measurement, labanimalid, rfid, matches("sha(0[189]|10)")) %>% 
  pivot_longer(cols = matches("sha")) %>% 
  mutate(measurement = ifelse(grepl("date", name), "date", measurement),
         name = ifelse(grepl("date", name), gsub("date_", "", name), name)) %>% 
  distinct() %>% 
  full_join(sha_rawgwas_1_08_10, by = c("labanimalid" = "subject", "measurement" = "name", "name" = "exp")) %>% # join to raw file and fill NA's
  mutate(cohort = coalesce(cohort.x, cohort.y),
         rfid = coalesce(rfid.x, rfid.y)) %>% 
  select(-matches("[.][xy]$")) %>% 
  mutate(value = coalesce(value, as.character(raw))) %>%
  group_by(rfid) %>%
  # fill(cohort, .direction = "downup") %>% 
  # fill(sex, .direction = "downup") %>% 
  fill(labanimalid, .direction = "downup") %>%
  fill(box, .direction = "downup") %>% 
  ungroup() %>% 
  select(-raw) %>%
  distinct
  
## potentially needs to be corrected by the compromised table, bc they're naives?
# use to QC sha_rawgwas %>% subset(is.na(rfid)) %>% distinct(subject, cohort) %>% mutate(subject = parse_number(subject)) %>% select(subject) %>% unlist() %>% paste0(collapse = "|") %>% cat 

sha_gwasprelim2 <- sha_gwasprelim2 %>% 
  subset(measurement %in% c("inactive", "rewards", "date"))

sha_gwastraits <- sha_gwasprelim2 %>% 
  pivot_wider(names_from = c("measurement", "name"), values_from = "value") %>% 
  left_join(openxlsx::read.xlsx("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/U01/Olivier_George_U01DA043799 (Cocaine)/excel_and_csv_files/cocaine_sha_to_presses_final.xlsx") %>%  # join to TO file 
              clean_names %>% 
              mutate_at(vars(starts_with("sha")), as.numeric) %>% 
              select(labanimalid, sex, matches("sha(0[189]|10)")) %>% 
              select_all(~gsub("(sha.*)", "timeout_\\1", tolower(.))) , by = c("labanimalid"))

# calculate esc
sha_gwastraits <- sha_gwastraits %>% 
  mutate_at(vars(starts_with("rewards"), starts_with("inactive")), as.numeric) %>% 
  mutate_at(vars(matches("rewards_sha(0[89]|10)")), list(esc = ~.-rewards_sha01)) %>%
  mutate(mean_sha_delta_esc_08_10 = rowMeans(select(., ends_with("_esc")), na.rm = T)) %>% 
  mutate(mean_sha_08_10 = rowMeans(select(., matches("rewards_sha(0[89]|10)")), na.rm = T)) %>% 
  mutate(mean_inactive_sha_08_10 = rowMeans(select(., matches("inactive_sha(0[89]|10)")), na.rm = T)) %>% 
  mutate(mean_to_sha_08_10 = rowMeans(select(., matches("timeout_sha(0[89]|10)")), na.rm = T))


######
## LGA
######

lga_rawgwas <- read.csv("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/U01/Olivier_George_U01DA043799 (Cocaine)/excel_and_csv_files/cocaine_lga_raw_c01_11_oldnewdirs.csv", stringsAsFactors = F) %>%  # join to raw file 
  mutate_at(vars(one_of("box", "session_duration")), as.character) %>% pivot_longer(cols = where(is.numeric), values_to = "raw") %>% 
  select(subject, box, exp, name, raw, cohort) %>% 
  left_join(cocaine_metadata_df %>% select(labanimalid, rfid), by = c("subject" = "labanimalid")) 

lga_rawgwas %>% subset(is.na(rfid)) # validate the rfid ## there are ~30 cases that need to be fixed, but plan to come back to this after deadline 

## focus on the ones that are used for calculations
# sessions 1, 12-15
lga_rawgwas_1_12_14 <- lga_rawgwas %>% 
  subset(exp %in% c("lga01", "lga12", "lga13", "lga14"))

lga_rawgwas_1_12_14 %>% subset(is.na(rfid))

lga_rawgwas_1_12_14 <- lga_rawgwas_1_12_14 %>% 
  rowwise() %>%
  mutate(subject = ifelse(is.na(rfid), gsub("F", "M", subject), subject)) %>% # fix the sex in the ID and then rejoin to metadata
  ungroup() %>%
  group_by(subject) %>% 
  fill(rfid, .direction = 'downup') %>% 
  ungroup()

## potentially needs to be corrected by the compromised table, bc they're naives?
# use to QC lga_rawgwas %>% subset(is.na(rfid)) %>% distinct(subject, cohort) %>% mutate(subject = parse_number(subject)) %>% select(subject) %>% unlist() %>% paste0(collapse = "|") %>% cat 

lga_gwasprelim2 <- oliviercocaine_excel_all %>% 
    select(cohort, measurement, labanimalid, rfid, matches("lga(01|1[2-4])")) %>% 
  pivot_longer(cols = matches("lga")) %>% 
  mutate(measurement = ifelse(grepl("date", name), "date", measurement),
         name = ifelse(grepl("date", name), gsub("date_", "", name), name)) %>% 
  distinct() %>% 
  full_join(lga_rawgwas_1_12_14, by = c("labanimalid" = "subject", "measurement" = "name", "name" = "exp")) %>% # join to raw file and fill NA's
  mutate(cohort = coalesce(cohort.x, cohort.y),
         rfid = coalesce(rfid.x, rfid.y)) %>% 
  select(-matches("[.][xy]$")) %>% 
  mutate(value = coalesce(value, as.character(raw))) %>%
  group_by(rfid) %>%
  # fill(cohort, .direction = "downup") %>% 
  # fill(sex, .direction = "downup") %>% 
  fill(labanimalid, .direction = "downup") %>%
  fill(box, .direction = "downup") %>% 
  ungroup() %>% 
  select(-raw) %>%
  distinct

lga_gwasprelim2 <- lga_gwasprelim2 %>% 
  subset(measurement %in% c("inactive", "rewards", "date"))

lga_gwastraits <- lga_gwasprelim2 %>% 
  pivot_wider(names_from = c("measurement", "name"), values_from = "value") %>% 
  left_join(openxlsx::read.xlsx("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/U01/Olivier_George_U01DA043799 (Cocaine)/excel_and_csv_files/cocaine_lga_to_presses_final.xlsx") %>%  # join to TO file 
              clean_names %>% 
              mutate_at(vars(starts_with("lga")), as.numeric) %>% 
              select(labanimalid, sex, matches("lga(01|1[2-4])")) %>% 
              select_all(~gsub("(lga.*)", "timeout_\\1", tolower(.))) , by = c("labanimalid"))

# calculate esc
lga_gwastraits <- lga_gwastraits %>% 
  mutate_at(vars(starts_with("rewards"), starts_with("inactive")), as.numeric) %>% 
  mutate_at(vars(matches("rewards_lga1[234]")), list(esc = ~.-rewards_lga01)) %>%
  mutate(mean_lga_delta_esc_12_14 = rowMeans(select(., ends_with("_esc")), na.rm = T)) %>% 
  mutate(mean_lga_12_14 = rowMeans(select(., matches("rewards_lga(1[2-4])")), na.rm = T)) %>% 
  mutate(mean_inactive_lga_12_14 = rowMeans(select(., matches("inactive_lga(1[2-4])")), na.rm = T)) %>% 
  mutate(mean_to_lga_12_14 = rowMeans(select(., matches("timeout_lga(1[2-4])")), na.rm = T))
  



######
## PR 
######


pr_rawgwas <- read.csv("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/U01/Olivier_George_U01DA043799 (Cocaine)/excel_and_csv_files/cocaine_pr_raw_c01_11_oldnewdirs.csv", stringsAsFactors = F) %>%  # join to raw file 
  mutate_at(vars(one_of("box")), as.character) %>% pivot_longer(cols = where(is.numeric), values_to = "raw") %>% 
  select(subject, box, exp, name, raw, cohort) %>% 
  left_join(cocaine_metadata_df %>% select(labanimalid, rfid), by = c("subject" = "labanimalid")) 

pr_rawgwas %>% subset(is.na(rfid)) # validate the rfid ## there are ~30 cases that need to be fixed, but plan to come back to this after deadline 

# all sessions are important, so no subsetting

# use this for all combos to find the animals pr_rawgwas %>% subset(cohort == "C09"&box == "16") %>% distinct(subject)
pr_rawgwas <- pr_rawgwas %>% 
  mutate(subject = ifelse(subject == "F951", "M951", subject),
         subject = ifelse(subject == "F952", "M952", subject),
         subject = ifelse(subject == "M3", "F3", subject)) %>% 
  select(-rfid) %>% # remove and rejoin to the metadata
  left_join(cocaine_metadata_df %>% select(labanimalid, rfid), by = c("subject" = "labanimalid")) # validate the rfid 

pr_rawgwas %>% subset(is.na(rfid))

## potentially needs to be corrected by the compromised table, bc they're naives?
# use to QC pr_rawgwas %>% subset(is.na(rfid)) %>% distinct(subject, cohort) %>% mutate(subject = parse_number(subject)) %>% select(subject) %>% unlist() %>% paste0(collapse = "|") %>% cat 

pr_gwasprelim2 <- oliviercocaine_excel_all %>% 
  select(cohort, measurement, labanimalid, rfid, matches("pr(0[1-3])")) %>% 
  pivot_longer(cols = matches("pr")) %>% 
  mutate(measurement = ifelse(grepl("date", name), "date", measurement),
         name = ifelse(grepl("date", name), gsub("date_", "", name), name)) %>% 
  distinct() %>% 
  full_join(pr_rawgwas, by = c("labanimalid" = "subject", "measurement" = "name", "name" = "exp")) %>% # join to raw file and fill NA's
  mutate(cohort = coalesce(cohort.x, cohort.y),
         rfid = coalesce(rfid.x, rfid.y)) %>% 
  select(-matches("[.][xy]$")) %>% 
  mutate(value = coalesce(value, as.character(raw))) %>%
  group_by(rfid) %>%
  # fill(cohort, .direction = "downup") %>% 
  # fill(sex, .direction = "downup") %>% 
  fill(labanimalid, .direction = "downup") %>%
  fill(box, .direction = "downup") %>% 
  ungroup() %>% 
  select(-raw) %>%
  distinct

# check for any missing labanimalid's or rfid's
pr_gwasprelim2 %>% naniar::vis_miss()

pr_gwasprelim2 <- pr_gwasprelim2 %>% 
  subset(measurement %in% c("pr_breakpoint", "rewards", "date"))

# don't need to join to timeout, no timeout value
pr_gwastraits <- pr_gwasprelim2 %>% 
  pivot_wider(names_from = c("measurement", "name"), values_from = "value") 


# calculate max 
pr_gwastraits <- pr_gwastraits %>% 
  mutate(pr_max_02_03 = pmax(rewards_pr02,rewards_pr03),
         pr_max_02_03_breakpoint = pmax(pr_breakpoint_pr02, pr_breakpoint_pr03)) %>% 
  rename("pr_01_sha" = "rewards_pr01",
         "pr_02_lga" = "rewards_pr02" ,
         "pr_03_postshock" = "rewards_pr03",
         "pr_01_sha_breakpoint" = "pr_breakpoint_pr01",
         "pr_02_lga_breakpoint" = "pr_breakpoint_pr02" ,
         "pr_03_postshock_breakpoint" = "pr_breakpoint_pr03")

######
## SHOCK 
######

shock_rawgwas <- read.csv("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/U01/Olivier_George_U01DA043799 (Cocaine)/excel_and_csv_files/cocaine_shock_raw_c01_11_oldnewdirs.csv", stringsAsFactors = F) %>%  # join to raw file 
  mutate_at(vars(one_of("box", "session_duration")), as.character) %>% pivot_longer(cols = where(is.numeric), values_to = "raw") %>% 
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

shock_rawgwas %>% subset(is.na(rfid))

# since there are no excel values using the new formula, will be using the raw data as sole source 

# calculate the shock related values 
shock_gwastraits <- shock_rawgwas %>% 
  subset(exp %in% c("preshock", "shock03")&name == "rewards_w_postshock1") %>%
  group_by(rfid, subject, exp) %>% 
  fill(raw, .direction = "downup") %>% 
  distinct(cohort, rfid, subject, box, exp, room, raw) %>% 
  rename("labanimalid" = "subject") %>% 
  group_by(labanimalid) %>% 
  fill(box, room) %>% 
  ungroup() %>% 
  distinct() %>% 
  pivot_wider(names_from = exp, values_from = raw) %>% 
  mutate(shock_03_pre = (shock03-preshock)/preshock)

## XX join to the lga object to calculate another trait  
shock_gwastraits <- shock_rawgwas_traits %>% 
  left_join(XX lga_gwastraits %>% , by = "labanimalid") %>% 
  mutate(shock_03_avg1h = (shock03-avghour1_lgalast3)/avghour1_lgalast3)









# bind all objects into one object and create csv file 
gwas_prelim2 <- full_join()
  