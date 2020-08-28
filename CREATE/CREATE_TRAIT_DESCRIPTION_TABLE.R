## CREATE TRAIT DESCRIPTION TABLE

trait_description_table <- data.frame(
  variable = c("bodyweight"),
  variable_description = c("abdominal fat weight measured in grams"),
  cov_or_trait = "trait",
  general_category = c("Whole Body"),
  trait_name = c("Body Weight"),
  subtrait_name = c("abdominal fat weight")
) %>%
  mutate_all(as.character) %>% 
  bind_rows(rewards %>% 
              distinct(exp) %>%
              subset(!grepl("LGA([2]\\d+|[1][5-9])", exp)) %>% 
              mutate(exp = gsub("(\\D+)(\\d+)", "\\1_\\2", exp) %>% tolower) %>% 
              mutate(variable = paste0(exp, "_rewards"),
                     variable_description = "number of cocaine infusions received during", 
                     cov_or_trait = "trait",
                     general_category = "Behavior", 
                     trait_name = "Consumption level",
                     subtrait_name = "cocaine preference") %>% 
              arrange(variable) %>% 
              mutate(variable_description = case_when(
                grepl("lga.*_rewards", variable) ~ paste(variable_description, "day",
                                                parse_number(variable), "of long access self administration (6 hours)"),
                grepl("sha.*_rewards", variable) ~ paste(variable_description, "day",
                                               parse_number(variable), "of short access self administration (2 hours)"),
                grepl("^shock.*_rewards", variable) ~ paste(variable_description, "foot shocks of",
                                                  paste0("0.", parse_number(variable), " mA")),
                grepl("^pr_.*_rewards", variable) ~ paste(variable_description,  "day",
                                                  parse_number(variable), "of progressive ratio at breakpoint"),
                grepl("^preshock.*_rewards", variable) ~ "NA"))
            ) %>% 
  bind_rows(rewards %>% 
              distinct(exp) %>% 
              subset(!grepl("LGA([2]\\d+|[1][5-9])", exp)) %>%
              mutate(exp = gsub("(\\D+)(\\d+)", "\\1_\\2", exp) %>% tolower) %>% 
              mutate(variable = paste0(exp, "_left"),
                     variable_description = "number of active lever presses during", 
                     cov_or_trait = "trait",
                     general_category = "Behavior", 
                     trait_name = "Consumption level",
                     subtrait_name = "cocaine preference") %>% 
              arrange(variable) %>% 
              mutate(variable_description = case_when(
                grepl("lga.*_left", variable) ~ paste(variable_description, "day",
                                                         parse_number(variable), "of long access self administration (6 hours)"),
                grepl("sha.*_left", variable) ~ paste(variable_description, "day",
                                                         parse_number(variable), "of short access self administration (2 hours)"),
                grepl("^shock.*_left", variable) ~ paste(variable_description, "foot shocks of",
                                                            paste0("0.", parse_number(variable), " mA")),
                grepl("^pr_.*_left", variable) ~ paste(variable_description,  "day",
                                                          parse_number(variable), "of progressive ratio at breakpoint"),
                grepl("^preshock.*_left", variable) ~ "NA"))
  ) %>% 
  bind_rows(rewards %>% 
              distinct(exp) %>%
              subset(!grepl("LGA([2]\\d+|[1][5-9])", exp)) %>% 
              mutate(exp = gsub("(\\D+)(\\d+)", "\\1_\\2", exp) %>% tolower) %>% 
              mutate(variable = paste0(exp, "_right"),
                     variable_description = "number of inactive lever presses during", 
                     cov_or_trait = "trait",
                     general_category = "Behavior", 
                     trait_name = "Consumption level",
                     subtrait_name = "cocaine preference") %>% 
              arrange(variable) %>% 
              mutate(variable_description = case_when(
                grepl("lga.*_right", variable) ~ paste(variable_description, "day",
                                                      parse_number(variable), "of long access self administration (6 hours)"),
                grepl("sha.*_right", variable) ~ paste(variable_description, "day",
                                                      parse_number(variable), "of short access self administration (2 hours)"),
                grepl("^shock.*_right", variable) ~ paste(variable_description, "foot shocks of",
                                                         paste0("0.", parse_number(variable), " mA")),
                grepl("^pr_.*_right", variable) ~ paste(variable_description,  "day",
                                                       parse_number(variable), "of progressive ratio at breakpoint"),
                grepl("^preshock.*_right", variable) ~ "NA"))
  ) %>% 
  bind_rows(rewards %>% 
              distinct(exp) %>% 
              subset(!grepl("LGA([2]\\d+|[1][5-9])", exp)) %>% 
              mutate(exp = gsub("(\\D+)(\\d+)", "\\1_\\2", exp) %>% tolower) %>% 
              mutate(variable = paste0(exp, "_age"),
                     variable_description = "age of animal on", 
                     cov_or_trait = "cov",
                     general_category = "NA", 
                     trait_name = "NA",
                     subtrait_name = "NA") %>% 
              arrange(variable) %>% 
              mutate(variable_description = case_when(
                grepl("lga.*_age", variable) ~ paste(variable_description, "day",
                                                       parse_number(variable), "of long access self administration (6 hours)"),
                grepl("sha.*_age", variable) ~ paste(variable_description, "day",
                                                       parse_number(variable), "of short access self administration (2 hours)"),
                grepl("^shock.*_age", variable) ~ paste(variable_description, "foot shocks of",
                                                          paste0("0.", parse_number(variable), " mA")),
                grepl("^pr_.*_age", variable) ~ paste(variable_description,  "day",
                                                        parse_number(variable), "of progressive ratio at breakpoint"),
                grepl("^preshock.*_age", variable) ~ "NA"))
  ) %>% 
  bind_rows(rewards %>% 
              distinct(exp) %>% 
              subset(grepl("LGA([0]\\d+|[1][0-3])", exp)) %>% 
              mutate(exp = gsub("(\\D+)(\\d+)", "\\1_\\2", exp) %>% tolower) %>% 
              mutate(variable = paste0(exp, "_esc"),
                     variable_description = paste("difference between the number of cocaine infusions received on day 1 and day", 
                                                   parse_number(variable) + 1, "of long access self administration (6 hours)"), 
                     cov_or_trait = "trait",
                     general_category = "Behavior", 
                     trait_name = "Consumption level",
                     subtrait_name = "cocaine preference") %>% 
              arrange(variable)) %>% 
  bind_rows(data.frame(variable = c("sha_days_until_5", 
                 "sha_mean_last_3",
                 "lga_esc_mean_last_3",
                 "pr_mean_last_2"),
    variable_description = c("number of days of short access self administration (2 hours) to reach 5 or more cocaine infusions",
                             "average of the number of cocaine infusions during the last three sessions of short access self administration (2 hours)",
                             "average of the escalation values during the last three sessions of long access self administration (6 hours)",
                             "average of the number of cocaine infusions at breakpoint for the last two sessions of progressive ratio")) %>% 
      mutate_all(as.character) %>% 
              mutate(cov_or_trait = "trait",
                     general_category = "Behavior", 
                     trait_name = "Consumption level",
                     subtrait_name = "cocaine preference") %>% 
              arrange(variable)) %>% 
  bind_rows(rewards %>% 
              select(-c("rewards", "exp", "valid", "filename")) %>% 
              names %>% as.data.frame() %>% 
              mutate_all(as.character) %>% 
              rename("variable" = ".") %>% 
              mutate(cov_or_trait = "cov")) %>% 
  mutate(project_name = "u01_olivier_george_cocaine", 
         exp = "self_administration") 

trait_description_table %>% 
  openxlsx::write.xlsx(file = "data_dictionary_olivier_cocaine.xlsx")
