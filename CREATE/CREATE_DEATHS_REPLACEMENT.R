### EXTRACT PROBLEMATIC RAT (Excel provided by the Olivier lab)


### COCAINE ### 


setwd("~/Dropbox (Palmer Lab)/Olivier_George_U01/Rat Information/Cocaine")



ratinfo_list_deaths <- lapply(list.files(pattern = ".xlsx"), function(x){ 
                              x <- u01.importxlsx(x)[[2]] %>% 
                                clean_names() %>%
                                mutate(cohort = paste0("C", str_pad(cohort, 2, "left", "0")),
                                       rfid = as.character(rfid),
                                       tailmark = str_match(tail_mark, "(M|F)[[0-9]]+")[, 1],
                                       naive = ifelse(grepl("Naive", tail_mark, ignore.case = T), "yes", "no"))
                              if(!lubridate::is.POSIXct(x$day_excluded)){
                                x <- x %>% 
                                  mutate(datedropped = openxlsx::convertToDateTime(day_excluded)) %>% 
                                  select(-day_excluded)
                              }
                              else{
                                x <- x %>% 
                                  mutate(datedropped = day_excluded)
                              }
                              x <- x %>%
                                mutate_all(as.character) %>% 
                                mutate(datedropped = replace(datedropped, datedropped == "9/12/207", "2017-09-12")) %>%
                                rename("reasoning" = "reasoning") %>%
                                select(cohort, rfid, tailmark, naive, datedropped, reasoning)
                              return(x)
}) %>% rbindlist() %>% 
  mutate(tailmark = replace(tailmark, tailmark == "F959", "M959"))


ratinfo_list_replacements <- lapply(list.files(pattern = ".xlsx"), function(x){ 
  if(length(u01.importxlsx(x)) > 2 & !grepl("C10", x)){ # C10 No replacements, no deaths during surgery and other cohorts don't have this tab
  x <- u01.importxlsx(x)[[3]] %>% 
    clean_names() %>%
    mutate(cohort = paste0("C", str_pad(cohort, 2, "left", "0")),
           originalrat = str_match(original_rat, "(M|F)[[0-9]]+")[, 1],
           replacement = str_match(replaced_with, "(M|F)[[0-9]]+")[, 1],
           rfidreplacement = as.character(rfid_of_replaced)) 
  
  names(x) <- gsub("x5", "comment", names(x))

  return(x)
  }
}) %>% rbindlist(fill = T) %>%
    select(cohort, originalrat, replacement, rfidreplacement, comment)





# check to see that all rats in replacements list are accounted for in deaths table
if(ratinfo_list_replacements %>% subset(!originalrat %in% ratinfo_list_deaths$tailmark) %>% nrow() != 0){
  ratinfo_list_replacements %>% subset(!originalrat %in% ratinfo_list_deaths$tailmark)
}
# if 0, proceed to left the replacements to dead animals table 

setwd("~/Desktop/Database/csv files/u01_olivier_george_cocaine")
compromised_rats <- left_join(ratinfo_list_deaths, ratinfo_list_replacements %>% 
                                select(-cohort), by = c("tailmark" = "originalrat")) %>% 
  rename("death_comment" = "reasoning",
         "rfid_compromised" = "rfid") %>% 
  mutate(comment = str_to_title(comment)) %>% 
  mutate(labanimalid = ifelse(grepl("^Not Renumbered", comment), replacement, tailmark),
         rfid = ifelse(grepl("Renumbered", comment), rfidreplacement, rfid_compromised)) %>% 
  # mutate(labanimalid = coalesce(rfidreplacement, tailmark)) %>%  # labanimalid is the one ultimately used 
  mutate(death_comment = gsub("^ | $", "", death_comment))
# %>%  ## 08/07/2020 RUN THIS LINE OF CODE WHEN THE TWO COHORT 11 ANIMALS DEATH INFO IS RECORDED
  # write.csv("compromised_rats_n70_c01_10.csv")
  
  




# ratinfo_excel <- list.files(pattern = ".*Cocaine.*")
# ratinfo_list <- lapply(ratinfo_excel, u01.importxlsx)
# names(ratinfo_list) <- ratinfo_excel

# lapply(ratinfo_list, function(x){
#   x <- rbindlist()
# })
# 
# cocaine_deaths <- lapply(list, `[[`, 1) %>% bind_rows()
# cocaine_replacement <- lapply(list, `[[`, 1) %>% bind_rows()