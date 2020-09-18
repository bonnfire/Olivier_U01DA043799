# plotting olivier's data

##############################################
## plotting stripped phenotypes
##############################################
cocaine_phenotypes_merge_stripped %>% 






##############################################
## plotting excel self admin data
##############################################

setwd("~/Dropbox (Palmer Lab)/Palmer Lab/Bonnie Lin/github/Olivier_U01Cocaine/CREATE")
cocaine_indices_xl <- u01.importxlsx("Addiction indices for C01-C07 Cocaine.xlsx") %>% 
  lapply(function(x){
    if(grepl("F701", x[1,1])){ # if cohort 7 table
      names(x)[1] <- "rats"
      names(x)[44] <- "rats"
    }
    x <- x %>% clean_names
    i <- 1
    if(names(x)[1] == "x1"){
      names(x) <- x[i,]
      x <- x %>% clean_names
      x <- x[-i, ]
      i = i + 1 # go down the rows until the test expression is false
    }
    
    # extract the last section of data
    last_rats_section_index <- grep("rat", names(x), ignore.case = T) %>% tail(1) # find the last section
    x <- x[, last_rats_section_index:ncol(x)]

    # rename the variables for uniformity
    names(x) <- c("labanimalid", "escalation_index", "pr_index", "shock_index", "addiction_index")
    
    # clean up the data
    x <- x %>% 
      mutate(labanimalid = str_extract(toupper(labanimalid), "[MF]\\d+")) %>% 
      subset(grepl("[MF]\\d+", labanimalid)) %>% 
      mutate_at(vars(ends_with("index")), as.numeric) %>% 
      mutate(sex = str_extract(labanimalid, "[MF]"))
    
    return(x)
  }) %>% rbindlist(idcol = "cohort") 

## XX figure out how to join box information

pdf("olivier_cocaine_indices_xl.pdf",onefile = T)

plot_list_main = list()
plot_list = list()
plot_list2 = list()

phenotypes_xl_vars <- grep("index$", names(cocaine_indices_xl), value = T)
for (i in seq_along(phenotypes_xl_vars)){
  
  plot_list_main[[i]] <- cocaine_indices_xl %>% 
    ggplot() + 
    geom_density(aes_string(phenotypes_xl_vars[i])) +
    theme(axis.text=element_text(size=12))
  plot_list[[i]] <- cocaine_indices_xl %>% 
    ggplot(aes(x = cohort)) + 
    geom_boxplot(aes_string(y = phenotypes_xl_vars[i])) + 
    theme(axis.text=element_text(size=12), axis.text.x = element_text(angle = 45))
  plot_list2[[i]] <- cocaine_indices_xl %>% 
    ggplot() + 
    geom_density(aes_string(phenotypes_xl_vars[i])) + 
    facet_grid(rows = vars(cohort)) + 
    theme(axis.text=element_text(size=12), axis.text.x = element_text(angle = 45))
  
  print(plot_list_main[[i]])
  print(plot_list[[i]])
  print(plot_list2[[i]])
  
}
dev.off()


# pdftools::pdf_combine() # needed to open terminal and run , before install.packages("pdftools") 09/02/2020


## PLOT 2 
## Raw traits from Excel data, calculate mean and max locally, and plot

pdf("olivier_cocaine_rawtraits_xl.pdf",onefile = T)

plot_list_main = list()
plot_list = list()
plot_list2 = list()

phenotypes_xl_vars <- grep("cohort|labanimalid", names(cocaine_intermediates_xl_bind), value = T, invert = T)
for (i in seq_along(phenotypes_xl_vars)){
  
  plot_list_main[[i]] <- cocaine_intermediates_xl_bind %>% 
    ggplot() + 
    geom_density(aes_string(phenotypes_xl_vars[i])) +
    theme(axis.text=element_text(size=12))
  plot_list[[i]] <- cocaine_intermediates_xl_bind %>% 
    ggplot(aes(x = cohort)) + 
    geom_boxplot(aes_string(y = phenotypes_xl_vars[i])) + 
    theme(axis.text=element_text(size=12), axis.text.x = element_text(angle = 45))
  plot_list2[[i]] <- cocaine_intermediates_xl_bind %>% 
    ggplot() + 
    geom_density(aes_string(phenotypes_xl_vars[i])) + 
    facet_grid(rows = vars(cohort)) + 
    theme(axis.text=element_text(size=12), axis.text.x = element_text(angle = 45))
  
  print(plot_list_main[[i]])
  print(plot_list[[i]])
  print(plot_list2[[i]])
  
}
dev.off()

##############################################
## plotting excel self admin data
##############################################

pdf("cohort1_olivier_selfadmin.pdf",onefile = T)

plot_list = list()

values <- grep("^(def|agg)", names(cohort1$tirritability), value = T)
df <- cohort1$tirritability
for (i in seq_along(values)){

  plot_list[[i]] <- ggplot(df, aes(x=rat))+ geom_point(aes_string(y = values[i])) + labs(title = values[i])
  print(plot_list[[i]])

}

df2 <- cohort1$tselfadmin
values2 <- grep(pattern = "^(?!rf|date|comment|lab|cohort)", names(cohort1$tselfadmin), perl = T, value = T)
for (i in seq_along(values2)){
  
  plot_list[[i]] <- ggplot(df2, aes(x=labanimalid))+ geom_point(aes_string(y = values2[i])) + labs(title = values2[i]) 
  print(plot_list[[i]])
  
}
dev.off()

##### 

# redo pdf to create one table with all cohorts


list <- list("cohort1" = cohort1, 
             "cohort2" = cohort2,
             "cohort3" = cohort3,
             "cohort4" = cohort4,
             "cohort5" = cohort5,
             "cohort7" = cohort7,
             "cohort8" = cohort8)

selfadmin_list <- lapply(list, `[[`, 1)
selfadmin_df <- bind_rows(selfadmin_list) # check dplyr documentation for the column compatibility 

datadic_list <- lapply(list, `[[`, 2)
datadictionary_df <- bind_rows(datadic_list, .id = "cohort")
# Create subset for unique entries 
datadictionary_df_subset <- subset(datadictionary_df, !duplicated(Rat)) ## check in on [47] formatting cohort7	Shock	Real	30% of rewards associated with 0.3 mA footshook 

pdf("olivier_selfadmin2.pdf", onefile = T)
a = "Note from Giordano email: 
The numbers of LgA days varies because we want to keep running our rats before dissection and between treatments (PR or Shock). Since sometimes there are weekends or vacations the number of LgA days will be different on every cohort.
For the final analysis and to determine our addiction index for each rat we only use the first 14 days of long access.
Regarding the shock sessions you are right. We decided to do only one shock session (0.3mA)  since the animals were responding normally to the other 2 intensities.
The  preshock session in cohort 8 was a mistake by the technician. She ran a 6h session instead og a 1h that’s why the numbers look higher."
text(1,4,a, pos=4)

plot_list = list()
plot_list2 = list()

colnames(selfadmin_df) <- make.unique(names(selfadmin_df)) ## temporary fix, the colnames need fixing
wrapper <- function(x, ...) {
  paste(strwrap(x, ...), collapse = "\n")
}

selfadminmeasures <- grep(pattern = "^(?!rf|date|comment|lab|cohort)", names(selfadmin_df), perl = T, value = T)
for (i in seq_along(selfadminmeasures)){
  
  plot_list[[i]] <- ggplot(selfadmin_df, aes(x=cohort, color = cohort))+ geom_boxplot(aes_string(y = selfadminmeasures[i])) + labs(title = paste0(selfadminmeasures[i], "_SA_U01_Olivier", "\n", wrapper(datadictionary_df_subset$Description[i+1], width = 80)))
  plot_list2[[i]] <- ggplot(selfadmin_df, aes(color = cohort))+ geom_density(aes_string(selfadminmeasures[i])) + labs(title = paste0(selfadminmeasures[i], "_SA_U01_Olivier", "\n", wrapper(datadictionary_df_subset$Description[i+1], width = 80)))
  print(plot_list[[i]])
  print(plot_list2[[i]])
  
}


dev.off()



# ggplot(df, aes(x=rat, y = defensive))+ geom_point()

##############################################
## plotting raw self admin data
##############################################



# redo pdf to create one table with all cohorts

pdf("olivier_raw_selfadmin_bysex.pdf",onefile = T)

plot_list = list()
plot_list2 = list()

selfadminmeasures <- grep(pattern = "^(?!rfid|cohort|labanimalid|sex)", names(Olivier_Cocaine_C01_09), perl = T, value = T)
for (i in seq_along(selfadminmeasures)){
  
  plot_list[[i]] <- ggplot(Olivier_Cocaine_C01_09, aes(x = sex)) + 
    geom_boxplot(aes_string(y = selfadminmeasures[i]))
  plot_list2[[i]] <- ggplot(Olivier_Cocaine_C01_09, aes(color = sex)) + 
    geom_density(aes_string(selfadminmeasures[i])) 

  print(plot_list[[i]])
  print(plot_list2[[i]])
  
}

dev.off()



pdf("olivier_raw_selfadmin_bycohort.pdf",onefile = T)

plot_list = list()

selfadminmeasures <- grep(pattern = "^(?!rfid|cohort|labanimalid|sex)", names(Olivier_Cocaine_C01_09), perl = T, value = T)
for (i in seq_along(selfadminmeasures)){
  
  plot_list[[i]] <- Olivier_Cocaine_C01_09 %>% subset(!grepl("C0[69]|1[012]", cohort)) %>% # subset because the data are not loaded yet 
    ggplot(aes(x = cohort)) + 
    geom_boxplot(aes_string(y = selfadminmeasures[i]))
  plot_list2[[i]] <- Olivier_Cocaine_C01_09 %>% subset(!grepl("C0[69]|1[012]", cohort)) %>% 
    ggplot(aes(color = cohort)) + 
    geom_density(aes_string(selfadminmeasures[i])) 
  
  print(plot_list[[i]])
  print(plot_list2[[i]])
  
}


dev.off()

# cohort x sex by sha day
Olivier_Cocaine_C01_09 %>% subset(!grepl("C0[69]|1[012]", cohort)) %>% 
  select(cohort, sex, starts_with("SHA")) %>% 
  gather("exp", "rewards", -cohort, -sex) %>% 
  group_by(cohort, sex, exp) %>% 
  summarize(mean = mean(rewards, na.rm = T)) %>%  
  ggplot(aes(x = exp, y = mean)) +
  geom_point(aes(color = sex)) +
  geom_line(aes(color = sex, group = sex)) +
  facet_grid(rows = "cohort")




pdf("olivier_raw_selfadmin_bybox.pdf",onefile = T)

plot_list = list()

selfadminmeasures <- grep(pattern = "^(?!rfid|box|labanimalid|sex)", names(Olivier_Cocaine_C01_09), perl = T, value = T)
for (i in seq_along(selfadminmeasures)){
  
  plot_list[[i]] <- Olivier_Cocaine_C01_09 %>% subset(!grepl("C0[69]|1[012]", box)) %>% # subset because the data are not loaded yet 
    ggplot(aes(x = box)) + 
    geom_boxplot(aes_string(y = selfadminmeasures[i]))
  plot_list2[[i]] <- Olivier_Cocaine_C01_09 %>% subset(!grepl("C0[69]|1[012]", box)) %>% 
    ggplot(aes(color = box)) + 
    geom_density(aes_string(selfadminmeasures[i])) 
  
  print(plot_list[[i]])
  print(plot_list2[[i]])
  
}


dev.off()


### PLOTTING EXCEL INDICES


pdf("olivier_xl_selfadmin.pdf",onefile = T)

plot_list = list()

selfadminmeasures <- grep(pattern = "^(?!rfid|cohort|labanimalid|sex)", names(cocaine_gwas_xl_df), perl = T, value = T)
for (i in seq_along(selfadminmeasures)){
  
  plot_list[[i]] <- cocaine_gwas_xl_df %>% subset(!grepl("C0[69]|1[012]", cohort)) %>% # subset because the data are not loaded yet 
    ggplot(aes(x = cohort)) + 
    geom_boxplot(aes_string(y = selfadminmeasures[i]))
  plot_list2[[i]] <- cocaine_gwas_xl_df %>% subset(!grepl("C0[69]|1[012]", cohort)) %>% 
    ggplot(aes(color = cohort)) + 
    geom_density(aes_string(selfadminmeasures[i])) 
  
  print(plot_list[[i]])
  print(plot_list2[[i]])
  
}


dev.off()



pdf("olivier_xl_selfadmin.pdf",onefile = T)

plot_list = list()

selfadminmeasures <- grep(pattern = "^(?!rfid|cohort|labanimalid|sex)", names(cocaine_gwas_xl_df), perl = T, value = T)
for (i in seq_along(selfadminmeasures)){
  
  plot_list[[i]] <- cocaine_gwas_xl_df %>% subset(!grepl("C0[69]|1[012]", cohort)) %>% # subset because the data are not loaded yet 
    ggplot(aes(x = cohort)) + 
    geom_boxplot(aes_string(y = selfadminmeasures[i]))
  plot_list2[[i]] <- cocaine_gwas_xl_df %>% subset(!grepl("C0[69]|1[012]", cohort)) %>% 
    ggplot(aes(color = cohort)) + 
    geom_density(aes_string(selfadminmeasures[i])) 
  
  print(plot_list[[i]])
  print(plot_list2[[i]])
  
}


dev.off()