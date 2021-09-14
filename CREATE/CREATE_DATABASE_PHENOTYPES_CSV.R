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
  mutate(mean_sha_08_10 = rowMeans(select(., matches("rewards_sha(0[89]|10)$")), na.rm = T)) %>% 
  mutate(mean_inactive_sha_08_10 = rowMeans(select(., matches("inactive_sha(0[89]|10)$")), na.rm = T)) %>% 
  mutate(mean_to_sha_08_10 = rowMeans(select(., matches("timeout_sha(0[89]|10)$")), na.rm = T))

# check for dupes
sha_gwastraits %>% get_dupes(labanimalid)



## create the same for SHA01-03
sha_rawgwas_01_03 <- sha_rawgwas %>% 
  subset(exp %in% c("sha01", "sha02", "sha03"))

sha_rawgwas_01_03 %>% subset(is.na(rfid))

sha_rawgwas_01_03 <- sha_rawgwas_01_03 %>% 
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
sha_rawgwas_01_03 %>% subset(is.na(rfid))

# sha_rawgwas_01_03 %>% subset(box == "14"&cohort == "C04") # check for non dupes in sessions

sha_rawgwas_01_03 <- sha_rawgwas_01_03 %>% 
  mutate(subject = ifelse(subject == "F414", "F424", subject)
         # subject = ifelse(subject == "F124", "F524", subject), # technically an error but there is already an existing session for this animal
         # subject = ifelse(subject == "F123", "F523", subject), # technically an error but there is already an existing session for this animal
         # subject = ifelse(subject == "F125", "F525", subject), # technically an error but there is already an existing session for this animal
  ) %>% 
  subset(!(subject %in% c("F124", "F123", "F125")&cohort == "C05")) %>% 
  group_by(subject) %>% 
  fill(rfid, .direction = "downup") %>% 
  ungroup()


sha_gwasprelim3 <- oliviercocaine_excel_all %>% 
  select(cohort, measurement, labanimalid, rfid, matches("sha(0[123])")) %>% 
  pivot_longer(cols = matches("sha")) %>% 
  mutate(measurement = ifelse(grepl("date", name), "date", measurement),
         name = ifelse(grepl("date", name), gsub("date_", "", name), name)) %>% 
  distinct() %>% 
  full_join(sha_rawgwas_01_03, by = c("labanimalid" = "subject", "measurement" = "name", "name" = "exp")) %>% # join to raw file and fill NA's
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

sha_gwasprelim3 <- sha_gwasprelim3 %>% 
  subset(measurement %in% c("inactive", "rewards", "date"))

sha_gwastraits_2 <- sha_gwasprelim3 %>% 
  pivot_wider(names_from = c("measurement", "name"), values_from = "value") %>% 
  left_join(openxlsx::read.xlsx("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/U01/Olivier_George_U01DA043799 (Cocaine)/excel_and_csv_files/cocaine_sha_to_presses_final.xlsx") %>%  # join to TO file 
              clean_names %>% 
              mutate_at(vars(starts_with("sha")), as.numeric) %>% 
              select(labanimalid, sex, matches("sha(0[123])")) %>% 
              select_all(~gsub("(sha.*)", "timeout_\\1", tolower(.))) , by = c("labanimalid"))

# calculate traits for 01-03 sessions
sha_gwastraits_2 <- sha_gwastraits_2 %>% 
  mutate_at(vars(starts_with("rewards"), starts_with("inactive")), as.numeric) %>% 
  mutate(mean_sha_01_03 = rowMeans(select(., matches("rewards_sha(0[123])$")), na.rm = T)) %>% 
  mutate(mean_inactive_sha_01_03 = rowMeans(select(., matches("inactive_sha(0[123])$")), na.rm = T)) %>% 
  mutate(mean_to_sha_01_03 = rowMeans(select(., matches("timeout_sha(0[123])$")), na.rm = T)) %>% 
  select(labanimalid, mean_sha_01_03:mean_to_sha_01_03)

# check for dupes
sha_gwastraits_2 %>% get_dupes(labanimalid)







# join columns for iti data
sha_gwastraits_all <- left_join(sha_gwastraits, sha_gwasiti_df_agg, by = c("labanimalid" = "subject")) %>% 
  left_join(sha_gwastraits_2, by = "labanimalid") %>% 
  left_join(sha_gwasiti_df_agg_2, by = c("labanimalid" = "subject"))


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
              select_all(~gsub("(lga.*)", "timeout_\\1", tolower(.))) , by = c("labanimalid")) %>% 
  select(-box) %>% 
  distinct() 


#check for dupes 
lga_gwastraits %>% get_dupes(labanimalid)

lga_gwastraits <- lga_gwastraits %>%
  subset(!is.na(sex)) %>% 
  fill(-labanimalid, -rfid, -cohort, -sex) %>% 
  distinct() 

## leave as NA's for now 
lga_gwastraits %>% get_dupes(labanimalid)


# calculate esc
lga_gwastraits_2 <- lga_gwastraits %>% 
  add_count(labanimalid) %>% 
  subset(n==1) %>% 
  select(-n) %>% 
  mutate_at(vars(starts_with("rewards"), starts_with("inactive")), as.numeric) %>% 
  mutate_at(vars(matches("rewards_lga1[234]")), list(esc = ~.-rewards_lga01)) %>%
  mutate(mean_lga_delta_esc_12_14 = rowMeans(select(., ends_with("_esc")), na.rm = T)) %>% 
  mutate(mean_lga_12_14 = rowMeans(select(., matches("rewards_lga(1[2-4])$")), na.rm = T)) %>% 
  mutate(mean_inactive_lga_12_14 = rowMeans(select(., matches("inactive_lga(1[2-4])$")), na.rm = T)) %>% 
  mutate(mean_to_lga_12_14 = rowMeans(select(., matches("timeout_lga(1[2-4])$")), na.rm = T))
  
lga_gwastraits_2 %>% get_dupes(labanimalid)

# join columns for iti data
lga_gwastraits_all <- left_join(lga_gwastraits_2, lga_gwasiti_df_agg, by = c("labanimalid" = "subject")) 

lga_gwastraits_all %>% get_dupes(labanimalid)

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

# check dupes
pr_gwastraits %>% get_dupes(labanimalid)

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

# remove duplicate animals for now 
shock_gwastraits2 <- shock_gwastraits %>% 
  add_count(labanimalid) %>% 
  subset(n == 1) %>% 
  select(-n)

## join to the lga object to calculate another trait  
shock_gwastraits2 <- shock_gwastraits2 %>% 
  left_join(lga_gwastraits_all %>% 
              select(labanimalid, lga_rewards_first1hr), by = "labanimalid") %>% 
  mutate(shock_03_avg1h = (shock03-lga_rewards_first1hr)/lga_rewards_first1hr)


## get dupes
shock_gwastraits2 %>% distinct() %>% get_dupes(labanimalid)




## IRRITABILITY 

irritability_gwastraits <- irritability_df %>% 
  left_join(cocaine_metadata_df %>% select(rfid, d_o_b), by = "rfid") %>% # to get the age 
  mutate(irr_age = difftime(date_defensive01, d_o_b, units = c("days")) %>% as.numeric) %>% # covariate 
  mutate(irr_agg_change = aggressive02 - aggressive01,
         irr_def_change = defensive02 - defensive01,
         irr_total_change = (defensive02+aggressive02) - (defensive01+aggressive01)) %>% 
  select(rfid, irr_age:irr_total_change)






## covariate weight

#cocaine_metadata_df %>% select(cohort, rfid, matches("surgery")) %>% subset(!is.na(weight_2_surgery)) %>% View()
surgery_weights <- cocaine_metadata_df %>% 
  mutate(weight_2_surgery = coalesce(weight_2_sugery, weight_2_surgery)) %>% 
  mutate(surgery_weight = coalesce(weight_1_surgery, weight_2_surgery) %>% as.numeric) %>% 
  select(rfid, surgery_weight)
  
# to use first weight or surgery weight
  #mutate(weight_pre_shipment = coalesce(weight_1_pre_shipment, weight_1_11_26_2018) %>% 
  #          coalesce(weight_1_4_5_2018)) %>% 
  # select(cohort, rfid, matches("surgery|shipment")) %>%
  # mutate(surgery_weight = coalesce(weight_1_surgery, weight_2_surgery) %>% as.numeric) %>% 
  # select(cohort, rfid, surgery_weight, weight_pre_shipment) %>% 
  # pivot_longer(cols = surgery_weight:weight_pre_shipment, names_to = "weight_event", values_to = "weight") %>% 
  # mutate(weight_event_cohort = paste0(cohort, "_", weight_event)) %>% 
  # mutate(weight_event_cohort_2 = ifelse(grepl("surgery", weight_event), "surgery", "non_surgery")) %>% 
  # mutate(weight_event_cohort_3 = ifelse(grepl("surgery", weight_event), "surgery", weight_event_cohort))
# 
# 
# ggplot(surgery_weights, aes(x = weight, color = weight_event_cohort)) + geom_density() 
# ggplot(surgery_weights, aes(x = weight, group = weight_event_cohort, color = weight_event_cohort_2)) + geom_density() 
# ggplot(surgery_weights, aes(x = weight, group = weight_event_cohort, color = weight_event_cohort_3)) + geom_density() 

# check for no dupes 
surgery_weights %>% get_dupes(rfid)


# bind all objects into one object and create csv file 
gwas_prelim2 <- full_join(sha_gwastraits_all %>% 
                            select(cohort, labanimalid, rfid, sex, box, date_sha01, date_sha08, mean_sha_delta_esc_08_10:sha_sd_iti_na, mean_sha_01_03:mean_to_sha_01_03, sha_rewards_first10min_2:sha_rewards_first1hr_2), 
                          lga_gwastraits_all %>% 
                            select(labanimalid, date_lga12, mean_lga_delta_esc_12_14:lga_sd_iti_na), by = "labanimalid") %>% 
  full_join(shock_gwastraits2 %>% 
              select(labanimalid, shock03, shock_03_pre, shock_03_avg1h), by = "labanimalid") %>% 
  full_join(irritability_gwastraits, by = "rfid") %>% 
  left_join(oliviercocaine_excel_all %>% select(labanimalid, matches("date_shock03")) %>% distinct(), by = "labanimalid") %>% 
  full_join(pr_gwastraits %>% 
              select(labanimalid, pr_01_sha:pr_max_02_03_breakpoint), by = "labanimalid") %>% 
  left_join(cocaine_metadata_df %>% 
              distinct(rfid, d_o_b) %>% 
              mutate(d_o_b = as.Date(d_o_b)), by = "rfid") %>% 
  left_join(surgery_weights, by = "rfid") %>% 
  mutate_at(vars(matches("date")), list(age = ~ difftime(., d_o_b, units = c("days")) %>% as.numeric)) %>% 
  select_all(~gsub("^date_(.*age$)", "\\1", tolower(.))) %>% 
  select(-matches("date|d_o_b")) %>% 
  mutate(sha_coefficient_of_variation = sha_sd_iti/sha_mean_iti, 
         sha_coefficient_of_variation_2 = sha_sd_iti_2/sha_mean_iti_2, 
         esc_coefficient_of_variation = lga_sd_iti/lga_mean_iti) %>% 
  rename("iti_median_sha_08_10" = "sha_median_iti",
         "behavior_stability_sha_08_10" = "sha_sd_iti",
         "behavior_stability_sha_08_10_na" = "sha_sd_iti_na",
         
         "iti_median_sha_01_03" = "sha_median_iti_2",
         "behavior_stability_sha_01_03" = "sha_sd_iti_2",

         "iti_median_esc_12_14" = "lga_median_iti",
         "behavior_stability_esc_12_14" = "lga_sd_iti",
         "behavior_stability_esc_12_14_na" = "lga_sd_iti_na",
         
         "loading_phase_intake_sha_08_10" = "sha_rewards_first10min",
         "titration_phase_sha_08_10" = "sha_rewards_last60min",
         
         "loading_phase_intake_sha_01_03" = "sha_rewards_first10min_2",
         "titration_phase_sha_01_03" = "sha_rewards_last60min_2",
         
         "loading_phase_intake_esc_12_14" = "lga_rewards_first10min",
         "titration_phase_esc_12_14" = "lga_rewards_last60min"
         ) 




# take care of later 
# get_dupes(rfid) results 
gwas_prelim2 <- gwas_prelim2 %>% 
  subset(labanimalid != "F91") 


# addiction indices
gwas_prelim2_indices <- gwas_prelim2 %>%
  mutate(pr_max_02_03 = as.numeric(pr_max_02_03)) %>%
  
  # without sex
  group_by(cohort) %>% 
  mutate(mean_lga_delta_esc_12_14_zscore = scale(mean_lga_delta_esc_12_14),
         pr_max_02_03_zscore = scale(pr_max_02_03),
         shock03_zscore = scale(shock03)) %>% 
  ungroup() %>% 
  mutate(addiction_index_noSexZ = rowMeans(select(., ends_with("zscore")), na.rm = TRUE)) %>% 
  select(-matches("zscore")) %>% 
  
  # with sex
  group_by(cohort, sex) %>% 
  mutate(mean_lga_delta_esc_12_14_zscore = scale(mean_lga_delta_esc_12_14),
         pr_max_02_03_zscore = scale(pr_max_02_03),
         shock03_zscore = scale(shock03)) %>% 
  ungroup() %>% 
  mutate(add_index_sexcohortZ = rowMeans(select(., ends_with("zscore")), na.rm = TRUE)) %>% 
  select(-matches("zscore")) %>% 
  
  # without sex or cohort
  # mutate(add_ind_calc_noZ_shock_03 = rowMeans(select(., matches("(mean_lga_delta_esc_12_14|pr_max_02_03|shock03)$")), na.rm = TRUE)) %>% 
  
  # without sex or cohort
  mutate(mean_lga_delta_esc_12_14_zscore = scale(mean_lga_delta_esc_12_14),
         pr_max_02_03_zscore = scale(pr_max_02_03),
         shock03_zscore = scale(shock03)) %>% 
  mutate(add_ind_calc_noZ_shock_03 = rowMeans(select(., ends_with("zscore")), na.rm = TRUE)) %>% 
  select(-matches("zscore")) %>% 
  
  # repeat with a diff variable for shock
  group_by(cohort) %>% 
  mutate(mean_lga_delta_esc_12_14_zscore = scale(mean_lga_delta_esc_12_14),
         pr_max_02_03_zscore = scale(pr_max_02_03),
         shock03_zscore = scale(shock_03_avg1h)) %>% 
  ungroup() %>% 
  mutate(addiction_index_noSexZ_shock_03_avg1h = rowMeans(select(., ends_with("zscore")), na.rm = TRUE)) %>% 
  select(-matches("zscore")) %>% 
  
  group_by(cohort, sex) %>% 
  mutate(mean_lga_delta_esc_12_14_zscore = scale(mean_lga_delta_esc_12_14),
         pr_max_02_03_zscore = scale(pr_max_02_03),
         shock03_zscore = scale(shock_03_avg1h)) %>% 
  ungroup() %>% 
  mutate(add_index_palmer_shock_03_avg1h = rowMeans(select(., ends_with("zscore")), na.rm = TRUE)) %>% 
  select(-matches("zscore")) %>% 
  
  # without sex or cohort
  # mutate(add_ind_calc_noZ_shock_03_avg1h = rowMeans(select(., matches("(mean_lga_delta_esc_12_14|pr_max_02_03|shock_03_avg1h)$")), na.rm = TRUE)) 
  
  # without sex or cohort
  mutate(mean_lga_delta_esc_12_14_zscore = scale(mean_lga_delta_esc_12_14),
         pr_max_02_03_zscore = scale(pr_max_02_03),
         shock03_zscore = scale(shock_03_avg1h)) %>% 
  mutate(add_ind_calc_noZ_shock_03_avg1h = rowMeans(select(., ends_with("zscore")), na.rm = TRUE)) %>% 
  select(-matches("zscore")) %>% 
  
  ## assign NA's if any of the three values in the index are missing
  mutate_at(vars(matches("^add")), ~ ifelse(is.na(mean_lga_delta_esc_12_14)|is.na(pr_max_02_03)|is.na(shock_03_avg1h), NA, .))



  
 # calculate the add_ind with and without normalizing for sex (add_ind_calc and add_ind_calc-noSexZ, respectively) and without normalizing for sex or cohort (add_ind_calc-noZscoreatALL). 
  



# join to Olivier's addiction indices
  
gwas_prelim2_indices <- gwas_prelim2_indices %>% 
  left_join(openxlsx::read.xlsx("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/U01/Olivier_George_U01DA043799 (Cocaine)/excel_and_csv_files/cocaine_gwas_addind_06242021.xlsx") %>% 
              clean_names %>% 
              mutate(rfid = as.character(rfid)) %>%
              rename("add_ind_olivier" = "add_ind") %>% 
              select(rfid, add_ind_olivier), by = "rfid") %>% 
  select(-sex) %>% 
  left_join(cocaine_metadata_df %>% distinct(rfid, sex), by = "rfid") %>% 
  mutate_at(vars(-one_of("cohort", "labanimalid", "rfid", "box", "sex")), as.numeric) %>% 
  select(cohort, labanimalid, rfid, sex, box, everything()) %>% 
  mutate_if(is.numeric, list(~na_if(., Inf)))

gwas_prelim2_indices %>% naniar::vis_miss()


## rename for database organization 
gwas_prelim2_indices <- gwas_prelim2_indices %>% 
  rename("sha_mean_delta_esc_08_10" = "mean_sha_delta_esc_08_10",
         "sha_mean_08_10" = "mean_sha_08_10", 
         "sha_mean_inactive_08_10" = "mean_inactive_sha_08_10",
         "sha_mean_to_08_10" = "mean_to_sha_08_10",                     
         "sha_loading_phase_intake_08_10" = "loading_phase_intake_sha_08_10",        
         "sha_titration_phase_08_10" = "titration_phase_sha_08_10",             
         "sha_iti_median_08_10" = "iti_median_sha_08_10",   
         "sha_mean_iti" = "sha_mean_iti",                          
         "sha_behavior_stability_08_10" = "behavior_stability_sha_08_10",          
         "sha_behavior_stability_08_10_na" = "behavior_stability_sha_08_10_na",       
         "sha_mean_01_03" = "mean_sha_01_03",
         "sha_mean_inactive_01_03" = "mean_inactive_sha_01_03",               
         "sha_mean_to_01_03" = "mean_to_sha_01_03",                     
         "sha_loading_phase_intake_01_03" = "loading_phase_intake_sha_01_03",        
         "sha_titration_phase_01_03" = "titration_phase_sha_01_03", 
         "sha_iti_median_01_03" = "iti_median_sha_01_03",                  
         "sha_mean_iti_2" = "sha_mean_iti_2",                        
         "sha_behavior_stability_01_03" = "behavior_stability_sha_01_03",          
         "sha_sd_iti_na_2" = "sha_sd_iti_na_2",
         "sha_rewards_first1hr_2" = "sha_rewards_first1hr_2",                
         "lga_mean_delta_esc_12_14" = "mean_lga_delta_esc_12_14",              
         "lga_mean_12_14" = "mean_lga_12_14",                        
         "lga_mean_inactive_12_14" = "mean_inactive_lga_12_14",
         "lga_mean_to_12_14" = "mean_to_lga_12_14",                     
         "lga_loading_phase_intake_12_14" = "loading_phase_intake_esc_12_14",       
         "lga_titration_phase_12_14" = "titration_phase_esc_12_14",            
         "lga_iti_median_12_14" = "iti_median_esc_12_14",  
         "lga_mean_iti" = "lga_mean_iti",                          
         "lga_behavior_stability_12_14" = "behavior_stability_esc_12_14",          
         "lga_behavior_stability_12_14_na" = "behavior_stability_esc_12_14_na",       
         "shock_03" = "shock03", 
         "shock_03_pre" = "shock_03_pre",                          
         "shock_03_avg1h" = "shock_03_avg1h",                        
         "irr_age" = "irr_age",                              
         "irr_agg_change" = "irr_agg_change", 
         "irr_def_change" = "irr_def_change",                        
         "irr_total_change" = "irr_total_change",                      
         "pr_01_sha" = "pr_01_sha",                             
         "pr_02_lga" = "pr_02_lga",
         "pr_03_postshock" = "pr_03_postshock",                       
         "pr_01_sha_breakpoint" = "pr_01_sha_breakpoint",                  
         "pr_02_lga_breakpoint" = "pr_02_lga_breakpoint",                  
         "pr_03_postshock_breakpoint" = "pr_03_postshock_breakpoint", 
         "pr_max_02_03" = "pr_max_02_03",                          
         "pr_max_02_03_breakpoint" = "pr_max_02_03_breakpoint",               
         "weight_surgery" = "surgery_weight",  
         "sha_01_age" = "sha01_age",  
         "sha_08_age" = "sha08_age",  
         "lga_12_age" = "lga12_age",                             
         "shock_03_age" = "shock03_age",                           
         "pr_01_age" = "pr01_age",                             
         "pr_02_age" = "pr02_age",
         "pr_03_age" = "pr03_age",                             
         "sha_coefficient_of_variation" = "sha_coefficient_of_variation",         
         "sha_coefficient_of_variation_2" = "sha_coefficient_of_variation_2",        
         "lga_coefficient_of_variation" = "esc_coefficient_of_variation", 
         "addiction_index_noSexZ" = "addiction_index_noSexZ",                
         "add_index_sexcohortZ" = "add_index_sexcohortZ",                  
         "add_ind_calc_noZ_shock_03" = "add_ind_calc_noZ_shock_03",            
         "addiction_index_noSexZ_shock_03_avg1h" = "addiction_index_noSexZ_shock_03_avg1h",
         "add_index_palmer_shock_03_avg1h" = "add_index_palmer_shock_03_avg1h",       
         "add_ind_calc_noZ_shock_03_avg1h" = "add_ind_calc_noZ_shock_03_avg1h",      
         "add_ind_olivier" = "add_ind_olivier"
         )






write.csv(gwas_prelim2_indices, "~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/U01/Olivier_George_U01DA043799 (Cocaine)/excel_and_csv_files/tier1_gwas_prelim_n545_v2.csv", row.names = F)
  







  