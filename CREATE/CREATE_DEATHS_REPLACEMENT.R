### EXTRACT PROBLEMATIC RAT (Excel provided by the Olivier lab)


### COCAINE ### 


setwd("~/Dropbox (Palmer Lab)/Olivier_George_U01/Rat Information")



ratinfo_list <- u01.importxlsx("Rat Information - All Cocaine and Oxycodone.xlsx")
names(ratinfo_list) <- excel_sheets("Rat Information - All Cocaine and Oxycodone.xlsx")

# ratinfo_excel <- list.files(pattern = ".*Cocaine.*")
# ratinfo_list <- lapply(ratinfo_excel, u01.importxlsx)
# names(ratinfo_list) <- ratinfo_excel

lapply(ratinfo_list, function(x){
  x <- rbindlist()
})

cocaine_deaths <- lapply(list, `[[`, 1) %>% bind_rows()
cocaine_replacement <- lapply(list, `[[`, 1) %>% bind_rows()