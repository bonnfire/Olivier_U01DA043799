# plotting olivier's data

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
selfadmin_df <- bind_rows(selfadmin_list, .id = "cohort") # check dplyr documentation for the column compatibility

datadic_list <- lapply(list, `[[`, 2)
datadictionary_df <- bind_rows(datadic_list, .id = "cohort")
# Create subset for unique entries 
datadictionary_df_subset <- subset(datadictionary_df, !duplicated(Rat)) ## check in on [47] formatting cohort7	Shock	Real	30% of rewards associated with 0.3 mA footshook 

pdf("olivier_selfadmin.pdf", onefile = T)
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
  plot_list2[[i]] <- ggplot(selfadmin_df, aes(color = cohort, alpha = 0.8))+ geom_density(aes_string(selfadminmeasures[i])) + labs(title = paste0(selfadminmeasures[i], "_SA_U01_Olivier", "\n", wrapper(datadictionary_df_subset$Description[i+1], width = 80)))
  print(plot_list[[i]])
  print(plot_list2[[i]])
  
}


dev.off()



# ggplot(df, aes(x=rat, y = defensive))+ geom_point()
