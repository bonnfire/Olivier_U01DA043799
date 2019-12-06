# extract raw


### COCAINE
install.packages("remotes")
remotes::install_github("gkane26/rmedpc")
example_file = system.file("extdata", "example_data", package = "rmedpc")
import_medpc(example_file)

## OXYCODEINE