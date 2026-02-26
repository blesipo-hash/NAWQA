
# remotes::install_github('andschar/standartox')

library(standartox)
library(dplyr)

taxa <- stx_taxa() # taxonomic info

chem <- stx_chem() # chemical info

dat <- stx_data() # toxicity data

# join chem and tax info to tox data

full_dat <- dat %>% 
  left_join(chem, by = "cl_id") %>%
  left_join(taxa, by = "tl_id")
