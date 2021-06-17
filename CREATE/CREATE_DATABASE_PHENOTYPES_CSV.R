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
