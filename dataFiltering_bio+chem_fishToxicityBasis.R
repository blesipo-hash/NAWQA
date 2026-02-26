

library(finsyncR)
library(tidyverse)
library(data.table)
library(dataRetrieval)


# get fish data and get vector of families and genera for which to look for tox data
fish <- getFishData(dataType = "abun", agency = "USGS", standardize = "MGMS")

fish_names <- str_replace_all(colnames(fish)[23:ncol(fish)], "\\.", " ")
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


##### explore and subset/filter tox data: #####

# all toxicological endpoints in dataset (e.g. LC50, EC50, NOEC, etc..)
table(tox_data$endpoint)

# more general grouping for endpoints (e.g. XX50 which includes EC LC LD ED IC ID...)
table(tox_data$endpoint_group)
table(tox_data$endpoint[tox_data$endpoint_group == "XX50"])

# all effects (e.g. mortality, reproduction, growth, behavior...)
table(tox_data$effect)

# exposure type for toxicity test, e.g. static, renewal, flow-through
table(tox_data$exposure)

# duration of toxicity test (histogram only shows distribution of durations for LC50 studies)
hist(tox_data$duration[tox_data$duration_unit == "h" & 
                         tox_data$duration <= 168 &  
                         tox_data$duration >= 0 & 
                         tox_data$endpoint == "LC50"])

# units of duration (mostly in hours "h")
table(tox_data$duration_unit)

# taxonomic rank of organism
# you will want to aggregate (average) these values if you use tox data at the genus or family level
table(tox_data$rank)

# broad taxonomic grouping (e.g. "fish", "insect", "invertebrate")
table(tox_data$group)

# concentration type (e.g. "total", "active ingredient", "dissolved",...)
table(tox_data$concentration_type)

# concentration qualifier (i.e. is the "actual" value equal to, less than, or greater than the listed conc)
table(tox_data$qualifier)

### now that we see what we're working with, let's filter the dataset 
### to include only 96-hour LC50s for fish genera and families in our dataset

# first get LC50 data for fish, selecting only the LC50s reported in "hours"
# and any concentration unit that is possible to convert to ug/L (same unit as NAWQA data)
# only select data for fish families and chemicals (cas numbers) in our dataset

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



# chem data in long form to get chemical info (e.g. CAS num)
chem_long <- data.table::fread("../data/chemData_long_071025.csv") %>% mutate(USGSPCode = as.character(USGSPCode))
cas_dat <- data.frame(cas = parameterCdFile$casrn[parameterCdFile$parameter_cd %in% 
                                                    unique(chem_long$USGSPCode)],
                      USGSPCode = parameterCdFile$parameter_cd[parameterCdFile$parameter_cd %in% 
                                                                               unique(chem_long$USGSPCode)])

# out of the top 50 chemicals for which fish toxicity data were collected, there are NAWQA data for 29
fish_top50 <- names(sort(table(lc50_fish_family$cas), decreasing = T)[1:50])

length(which(cas_dat$cas %in% fish_top50))

chem_long2 <- chem_long %>%
  left_join(cas_dat, by = "USGSPCode")

# what is the NAWQA sampling effort for the chemicals with the most fish tox data? 
sort(table(chem_long2$CharacteristicName[chem_long2$cas %in% fish_top50]), decreasing = T)

# how many times were each of these chemicals actually detected?
fishChems_detection <- table(chem_long2$CharacteristicName[chem_long2$cas %in% fish_top50],chem_long2$ResultDetectionConditionText[chem_long2$cas %in% fish_top50])
sort(fishChems_detection[,1], decreasing = T) # 11 chems detected > 100 times

# looks like there are 13 chemicals that were measured many times (>2000) with 11 of them actually being detected > 100 times
fishChems_final <- paste("chem_", names(fishChems_detection[,1][fishChems_detection[,1] > 100]), sep = "")

# we want to find optimal set of chems to include, so we loop through all combinations of 5 or more of the chems
# and determine the number of sampling events where all chems in each subset were measured
fishChems_sets <- rje::powerSet(fishChems_final)[sapply(rje::powerSet(fishChems_final), length) >=5]

# get chemical data in wide form
chem <- data.table::fread("../data/chemData_wide_071025.csv") %>%
  rename("SiteNumber" = "MonitoringLocationIdentifier",
         "CollectionDate" = "ActivityStartDate",
         "Latitude_dd" = "LatitudeMeasure",
         "Longitude_dd" = "LongitudeMeasure")


# get all sampling activities where all 11 fish chemicals were measured
# do this also for each set 
unique_activities <- unique(chem$ActivityIdentifier)


# (NOTE!!! this takes a while to run so we will use parallel processing to speed things up
# IMPORTANT!!! Disable OpenBLAS multithreading before starting parallel engine
library(parallel)
RhpcBLASctl::blas_set_num_threads(1)

if (RhpcBLASctl::blas_get_num_procs() == 1) {
  # specify number of threads to leave free
  free_threads <- 6
  
  # get number of threads to use in mclapply
  n.cores <- parallel::detectCores() - free_threads
  
  # use an abbreviated chem df to reduce computing
  # only use chems in list
  fishCols <- c(colnames(chem)[1:5], fishChems_final)
  chem_fishDf <- chem %>% dplyr::select(fishCols)

  activitiesBySet <- mclapply(X = 1:length(fishChems_sets), FUN = \(i) {
    
    fishSet <- fishChems_sets[[i]]
    
    activs <- unique_activities[sapply(unique_activities, \(x) {
      activ <- as.data.frame(chem_fishDf[chem_fishDf$ActivityIdentifier == x,])
      check_shared <- all(!is.na(activ[, fishSet]))
      return(check_shared)})]
    
    return(activs)
    
  }, mc.cores = n.cores)


}

# now we need to decide how to balance the number of sampling events 
# with the number of chemicals out of the 11 important ones

# 28 sampling events with all 11 chems
sapply(activitiesBySet, \(x) length(x))[sapply(fishChems_sets, \(x) length(x) == 11)]

# up to 1152 sampling events with 10 chems if you exclude glyphosate
max10 <- max(sapply(activitiesBySet, \(x) length(x))[sapply(fishChems_sets, \(x) length(x) == 10)])
set10 <- fishChems_sets[which(sapply(activitiesBySet, \(x) length(x) == max10) & sapply(fishChems_sets, \(x) length(x) == 10))][[1]]

# up to 6071 sampling events with 9 chems
max9 <- max(sapply(activitiesBySet, \(x) length(x))[sapply(fishChems_sets, \(x) length(x) == 9)])
set9 <- fishChems_sets[which(sapply(activitiesBySet, \(x) length(x) == max9) & sapply(fishChems_sets, \(x) length(x) == 9))][[1]]

# up to 18089 sampling events with 8 chems
max8 <- max(sapply(activitiesBySet, \(x) length(x))[sapply(fishChems_sets, \(x) length(x) == 8)])
set8 <- fishChems_sets[which(sapply(activitiesBySet, \(x) length(x) == max8) & sapply(fishChems_sets, \(x) length(x) == 8))][[1]]

# up to 22154 sampling events with 7 chems
max7 <- max(sapply(activitiesBySet, \(x) length(x))[sapply(fishChems_sets, \(x) length(x) == 7)])
set7 <- fishChems_sets[which(sapply(activitiesBySet, \(x) length(x) == max7) & sapply(fishChems_sets, \(x) length(x) == 7))][[1]]

# up to 27018 sampling events with 6 chems
max6 <- max(sapply(activitiesBySet, \(x) length(x))[sapply(fishChems_sets, \(x) length(x) == 6)])
set6 <- fishChems_sets[which(sapply(activitiesBySet, \(x) length(x) == max6) & sapply(fishChems_sets, \(x) length(x) == 6))][[1]]

# up to 28171 sampling events with 5 chems
max5 <- max(sapply(activitiesBySet, \(x) length(x))[sapply(fishChems_sets, \(x) length(x) == 5)])
set5 <- fishChems_sets[which(sapply(activitiesBySet, \(x) length(x) == max5) & sapply(fishChems_sets, \(x) length(x) == 5))][[1]]



# since it has the largest jump in sample size let's go with `set8`
activities <- activitiesBySet[which(sapply(activitiesBySet, \(x) length(x) == max8) & sapply(fishChems_sets, \(x) length(x) == 8))][[1]]

# let's look at what other chemicals may also be shared by these sampling events and
# whether any of them are potentially toxic to fish
otherChems <- colnames(chem[chem$ActivityIdentifier %in% activities, 6:ncol(chem)])[which(apply(chem[chem$ActivityIdentifier %in% activities, 6:ncol(chem)], 2, \(x) all(!is.na(x))))]
otherChems_cas <- na.omit(unique(chem_long2$cas[chem_long2$CharacteristicName %in% str_replace(otherChems, "chem_", "")]))

otherChems_tox <- table(lc50_fish_family$cname[lc50_fish_family$cas %in% otherChems_cas])

# these 10 chemicals have >= 50 96-hour LC50s reported, so we can consider including them in dataset
fishChems_full <- names(otherChems_tox)[otherChems_tox >= 50]
fishChems_full_cas <- unique(lc50_fish_family$cas[lc50_fish_family$cname %in% fishChems_full])

# but first let's see how many times they were actually detected
fishChems_full_detection <- table(chem_long2$CharacteristicName[chem_long2$cas %in% fishChems_full_cas], 
                                  chem_long2$ResultDetectionConditionText[chem_long2$cas %in% fishChems_full_cas])[,1]
fishChems_full_detection
# all were detected >100 times so it is fine to include all 10

# reduce chemical data to only include sampling events with those 10 shared chemicals
# remove columns for other chems
colnames_fishChems <- paste("chem_", unique(chem_long2$CharacteristicName[chem_long2$cas %in% fishChems_full_cas]), sep = "")

chem_reduced <- chem %>% 
  filter(ActivityIdentifier %in% activities) %>%
  dplyr::select(do.call(c, flatten(list("ActivityIdentifier", 
                                        "CollectionDate", 
                                        "Latitude_dd", 
                                        "Longitude_dd", 
                                        "SiteNumber",
                                        colnames_fishChems))))

# get fish data
fish <- getFishData(dataType = "abun", agency = "USGS", standardize = "MGMS")

# get shared sites for chem+inverts and chem+fish
commonSites_fish <- intersect(chem_reduced$SiteNumber, fish$SiteNumber)

# reduce both datasets to only include shared sites
chem_fish_sharedSites <- chem_reduced %>% filter(SiteNumber %in% commonSites_fish)
fish_sharedSites <- fish %>% filter(SiteNumber %in% commonSites_fish)

# save filtered data
data.table::fwrite(chem_fish_sharedSites, "../data/chemData_filtered_fishToxBasis_012226.csv")
data.table::fwrite(fish_sharedSites, "../data/fishData_filtered_fishToxBasis_012226.csv")

