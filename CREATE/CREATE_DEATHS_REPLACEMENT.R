### EXTRACT PROBLEMATIC RAT (Excel provided by the Olivier lab)


### COCAINE ### 


setwd("~/Dropbox (Palmer Lab)/Olivier_George_U01/Rat Information")



ratinfo_list_deaths <- u01.importxlsx("Rat Information - All.xlsx")[[1]]
ratinfo_list_replacements <- u01.importxlsx("Rat Information - All.xlsx")[[2]]

ratinfo_list_deaths_processed <- ratinfo_list_deaths %>% 
  mutate(cohort = paste0("C", str_pad(Cohort, 2, "left", "0")),
         rfid = as.character(RFID), 
         tailmark = str_match(`Tail.Mark`, "(M|F)[[0-9]]+")[, 1],
         naive = ifelse(grepl("Naive", `Tail.Mark`, ignore.case = T), "yes", "no"),
         datedropped = openxlsx::convertToDateTime(`Day.Excluded`),
         datedropped = replace(datedropped, `Day.Excluded` == "9/12/207", "2017-09-12")) %>% 
  rename("reasoning" = "Reasoning") %>% 
  select(cohort, rfid, tailmark, naive, datedropped, reasoning )


ratinfo_list_replacements_processed <- ratinfo_list_replacements %>% 
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