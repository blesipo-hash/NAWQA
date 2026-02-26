library(data.table)
library(tidyverse)

##### prepare chemical and biological data for risk quotient calculations: #####

# read in fish density and chemical data
fish <- fread("../data/fishData_filtered_fishToxBasis_012226.csv") %>%
  mutate(CollectionDate = as_date(CollectionDate))
chem <- fread("../data/chemData_filtered_fishToxBasis_012226.csv") %>%
  mutate(CollectionDate = as_date(CollectionDate))


# aggregate data by fish sampling event. for each sampling event,
# average all chemical measurements collected within X days of fishing
buff <- 90

chem_list_fish <- list()

for (i in 1:nrow(fish)) {
  
  # get site of sampling event
  site <- fish[i,]$SiteNumber
  
  # get sampling date, calculate min/max sampling date using buffer
  sampDate <- fish[i,]$CollectionDate
  minDate <- sampDate - buff
  maxDate <- sampDate + buff
  
  # get all chem data between min and max dates
  chemDat <- chem %>% 
    filter(CollectionDate >= minDate,
           CollectionDate <= maxDate,
           SiteNumber == site)
  
  activities <- unique(chemDat$ActivityIdentifier)
  
  chemDat_avg <- NA
  
  if (nrow(chemDat) > 0) {
    
    chemDat_avg <- as.data.frame(matrix(colMeans(chemDat[,6:ncol(chemDat)]), ncol = ncol(chemDat) - 5))
    colnames(chemDat_avg) <- colnames(chemDat)[6:ncol(chemDat)]
    
    chemDat_avg$SiteNumber <- site
    chemDat_avg$CollectionDate <- sampDate
    chemDat_avg$ActivityIdentifiers <- str_flatten(activities, collapse = " | ")
    
    chemDat_avg <- chemDat_avg %>%
      relocate(any_of(c("SiteNumber", "CollectionDate", "ActivityIdentifiers")), 
               .before = starts_with("chem_"))
    
  }
  
  chem_list_fish[[i]] <- list(
    site = site,
    sampDate = sampDate,
    activities = activities,
    samplingEvents_chem = nrow(chemDat),
    chemDat_avg = chemDat_avg
  )
  
}

# number of fish sampling events with chemical data collected within `buff` days
length(which(sapply(chem_list_fish, \(x) all(!is.na(x$chemDat_avg)))))

# remove sampling events with no chem data within `buff` days
chem_list_fish_abbrev <- chem_list_fish[sapply(chem_list_fish, \(x) all(!is.na(x$chemDat_avg)))]

# average number of chemical sampling events within `buff` days
mean(na.omit(sapply(chem_list_fish_abbrev, \(x) x$samplingEvents_chem)))

# filter chem and fish data to only have sampling events where chemical and fish
# data were collected within `buff` days
chem_final <- do.call(rbind, lapply(chem_list_fish_abbrev, \(x) x$chemDat_avg))
fish_final <- fish %>%
  filter(paste(CollectionDate, SiteNumber) %in% 
           paste(chem_final$CollectionDate, chem_final$SiteNumber))

# some cells of fish data have `Inf` for density, replace these values with the average density for that species
which(apply(fish_final, 1, \(x) any(is.infinite(as.numeric(x[23:length(x)])))))

fish_final[,23:(ncol(fish_final))] <- apply(fish_final[,23:(ncol(fish_final))], 2, \(x) {
  
  if(any(is.infinite(x))) {
    
    meanConc <- mean(x[!is.infinite(x)])
    
    x[is.infinite(x)] <- meanConc
    
    return(x)
    
  } else return(x)
  
}) %>% as.data.frame()



#### now that bio and chem data are filtered and matched up,
#### we can start to join any relevant toxicity information
fish_names <- str_replace_all(colnames(fish_final)[23:ncol(fish_final)], "\\.", " ")
fish_tax <- fread("../data/FishTaxonomy_42024.csv")

fish_genera <- unique(fish_tax$Genus[fish_tax$Species %in% fish_names])
fish_families <-  unique(fish_tax$Family[fish_tax$Species %in% fish_names])

tox_data <- fread("../data/standartox_allData_090825.csv")

fish_genera[fish_genera %in% tox_data$genus]
fish_families[fish_families %in% tox_data$family]

length(fish_genera[fish_genera %in% tox_data$genus])/length(fish_genera)
length(fish_families[fish_families %in% tox_data$family])/length(fish_families)

# look at most frequently reported tox metrics (NOT LOEC,NOEC, etc...) to decide what standard to use
# 96-hr LC50 is most common, although might be possible to estimate additional 96-hr LC50s using Haber's rule
# to estimate 96-hr LC50s from non-96-hr LC50s (e.g. 24 or 48 hr LC50s), but this is a rabbit hole: 
# https://pubmed.ncbi.nlm.nih.gov/10963857/
# https://onlinelibrary.wiley.com/doi/pdfdirect/10.1111/j.1539-6924.2011.01663.x


# get all cas numbers for our 66 chemicals
cas_nums <- fread("../data/chemicalInfo.csv")[,c("Chemical name", "CAS Number")] %>% filter(`CAS Number` != "")

# select relevant LC50s
lc50_fish_family <- tox_data %>%
  filter(family %in% fish_families,
         endpoint %in% c("LC50", "LC50*"),
         duration_unit == "h",
         duration == "96",
         qualifier == "=",
         concentration_unit %in% c("g/l", "ppb", "g/m3"))


# now let's do some unit standardization 
# (some of these might not matter depending on how you filtered data in the previous step)
table(lc50_fish_family$concentration_unit)

# initialize new column for standardized (ug/L) concentration
lc50_fish_family$concentration_std <- NA

# g/l -> ug/L
lc50_fish_family$concentration_std[lc50_fish_family$concentration_unit == "g/l"] <- lc50_fish_family$concentration[lc50_fish_family$concentration_unit == "g/l"] * 1e6
  
# ppb -> ug/L (assume these are equal in water)
lc50_fish_family$concentration_std[lc50_fish_family$concentration_unit == "ppb"] <- lc50_fish_family$concentration[lc50_fish_family$concentration_unit == "ppb"]

# g/m3 -> ug/L
lc50_fish_family$concentration_std[lc50_fish_family$concentration_unit == "g/m3"] <- lc50_fish_family$concentration[lc50_fish_family$concentration_unit == "g/m3"] * 1e3


# aggregate (average) lc50s by fish family
cas_chemname_key <- fread("../data/cas_chemname_key.csv") # key for linking CAS nums and NAWQA characteristic names
fishChems_full_cas <- qs::qread("../data/fishChems_casNums.qs") # CAS nums of 10 chemicals for fish tox analyses

lc50_fishFam_agg <- lc50_fish_family %>%
  filter(cas %in% fishChems_full_cas) %>%
  group_by(cas, family) %>%
  reframe(lc50_fam = mean(concentration_std)) %>%
  left_join(chem_long2 %>% 
              dplyr::select(cas, CharacteristicName) %>% 
              group_by(cas) %>% slice(1) %>% ungroup(), 
            by = "cas")

# determine which fish families have data for all chemicals
fish_fam_table <- table(lc50_fishFam_agg$CharacteristicName, lc50_fishFam_agg$family)

# 6 families with 96-hr lc50 for all 10 chemicals
fishFams_final <- colnames(fish_fam_table)[apply(fish_fam_table, 2, \(x) all(x == 1))]

# which means we have toxicity data for 52% of all species in our dataset!
length(which(fish_names %in% fish_tax$Species[fish_tax$Family %in% fishFams_final]))/length(fish_names)

fish_species_final <- fish_names[fish_names %in% 
                                   fish_tax$Species[fish_tax$Family %in% 
                                                      fishFams_final]]

# filter fish data to only include species within those 6 families
fish_dat_tox <- fish_final %>%
  dplyr::select(colnames(fish_final)[1:23], 
                str_replace(fish_species_final, " ", "\\."))

##### calculate risk quotients for each species based on family-level toxicity
# create 3-dimensional array to hold risk quotients for each fish species at each sampling event
fish_rq_dat <- array(data = NA, 
                     dim = c(length(fish_species_final), 
                             length(fishChems_full), 
                             nrow(fish_dat_tox)),
                     dimnames = list(Species = str_replace(fish_species_final, " ", "\\."),
                                     Chemical = colnames(chem_final)[4:ncol(chem_final)],
                                     SamplingEvent = paste(fish_dat_tox$SiteNumber, fish_dat_tox$CollectionDate, sep = "|"))
)

# loop through chemicals and calculate RQ for each species, dividing concentration by the family-level lc50
for (i in 1:length(fish_rq_dat[1,,1])) {
  
  # chem name
  ch <- dimnames(fish_rq_dat)$Chemical[i]
  
  for (j in 1:length(fish_rq_dat[,1,1])) {
    
    # species and family name
    sp <- dimnames(fish_rq_dat)$Species[j]
    fam <- fish_tax$Family[fish_tax$Species == str_replace(sp, "\\.", " ")]
    
    # lc50 corresponing to chemical and family
    lc50 <- lc50_fishFam_agg$lc50_fam[lc50_fishFam_agg$family == fam & 
                                        paste("chem_", lc50_fishFam_agg$CharacteristicName, 
                                              sep = "") == ch]
    
    # divide environmental concentrations at each sampling event by lc50
    ch_RQs <- chem_final[,ch]/lc50
    
    fish_rq_dat[j,i,] <- ch_RQs
    
  }
  
}
                     

# sum RQs for each chemical to get TUs for each fish species at each sampling event
fish_tu_dat <- t(apply(fish_rq_dat, c(1,3), sum)) %>% as.data.frame

# re-join sampling event info
fish_tu_dat$SiteNumber <- sapply(rownames(fish_tu_dat), \(x) str_split_1(x,"\\|")[1])
fish_tu_dat$CollectionDate <- sapply(rownames(fish_tu_dat), \(x) str_split_1(x,"\\|")[2])

fish_tu_dat2 <- fish_tu_dat %>%
  relocate(any_of(c("SiteNumber","CollectionDate")))


qs::qsave(fish_rq_dat, "fish_riskQuotientArray_012326.qs")
fwrite(fish_tu_dat2, "../data/fish_toxicUnits_012326.csv")
fwrite(fish_dat_tox, "../data/fish_speciesWithToxicUnits_012326.csv")









