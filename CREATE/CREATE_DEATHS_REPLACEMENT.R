### EXTRACT PROBLEMATIC RAT (Excel provided by the Olivier lab)


### COCAINE ### 


setwd("~/Dropbox (Palmer Lab)/Olivier_George_U01/Rat Information")



ratinfo_list_deaths <- u01.importxlsx("Rat Information - All.xlsx")[[2]]
ratinfo_list_replacements <- u01.importxlsx("Rat Information - All.xlsx")[[1]]

ratinfo_list_deaths_processed <- ratinfo_list_deaths %>% 
  clean_names() %>% 
  mutate(cohort = paste0("C", str_pad(cohort, 2, "left", "0")),
         rfid = as.character(rfid), 
         tailmark = str_match(tail_mark, "(M|F)[[0-9]]+")[, 1],
         naive = ifelse(grepl("Naive", tail_mark, ignore.case = T), "yes", "no"),
         datedropped = openxlsx::convertToDateTime(day_excluded),
         datedropped = replace(datedropped, day_excluded == "9/12/207", "2017-09-12")) %>% 
  rename("reasoning" = "reasoning") %>% 
  select(cohort, rfid, tailmark, naive, datedropped, reasoning)


## XX COME BACK TO FIX THIS 
ratinfo_list_replacements_processed <- ratinfo_list_replacements %>% 
  clean_names() %>% ### XXX FIX THE NAMES AFTER THIS LINE
  mutate(cohort = paste0("C", str_pad(Cohort, 2, "left", "0")),
         originalrat = str_match(`Original Rat`, "(M|F)[[0-9]]+")[, 1],
         replacement = str_match(`Replaced with`, "(M|F)[[0-9]]+")[, 1],
         rfidreplacement = as.character(`RFID of Replaced`)) %>% 
  rename("comment" = `...5`) %>% 
  select(cohort, originalrat, replacement, rfidreplacement, comment)




# ratinfo_excel <- list.files(pattern = ".*Cocaine.*")
# ratinfo_list <- lapply(ratinfo_excel, u01.importxlsx)
# names(ratinfo_list) <- ratinfo_excel

# lapply(ratinfo_list, function(x){
#   x <- rbindlist()
# })
# 
# cocaine_deaths <- lapply(list, `[[`, 1) %>% bind_rows()
# cocaine_replacement <- lapply(list, `[[`, 1) %>% bind_rows()