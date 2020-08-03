## CREATE TRAIT DESCRIPTION TABLE

rewards %>% 
  distinct(exp) %>% 
  mutate(variable = paste0(exp, "_rewards"),
         variable_description = "NA", 
         cov_or_trait = "trait") %>% 
  arrange(variable) %>% 
  bind_rows(rewards %>% 
              select(-c("rewards", "exp", "valid", "filename")) %>% 
              names %>% as.data.frame() %>% 
              mutate_all(as.character) %>% 
              rename("variable" = ".") %>% 
              mutate(cov_or_trait = "cov")) %>% 
  mutate(project_name = "u01_olivier_george_cocaine", 
         exp = "self_administration")
