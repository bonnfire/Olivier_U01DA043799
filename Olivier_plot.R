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


list <- list("cohort1" = cohort1, 
             "cohort2" = cohort2,
             "cohort3" = cohort3,
             "cohort4" = cohort4,
             "cohort5" = cohort5,
             "cohort7" = cohort7,
             "cohort8" = cohort8)

selfadmin_list <- lapply(list, `[[`, 1)
selfadmin_df <- bind_rows(selfadmin_list, .id = "cohort") # check dplyr documentation for the column compatibility
names(selfadmin_df)

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


# ggplot(df, aes(x=rat, y = defensive))+ geom_point()
