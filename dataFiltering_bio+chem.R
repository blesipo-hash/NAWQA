

library(finsyncR)
library(tidyverse)

# get chemical data
chem <- data.table::fread("../data/chemData_wide_071025.csv") %>%
  rename("SiteNumber" = "MonitoringLocationIdentifier",
         "CollectionDate" = "ActivityStartDate",
         "Latitude_dd" = "LatitudeMeasure",
         "Longitude_dd" = "LongitudeMeasure")

# find subset of data-rich chem sites
chem_summary <- apply(chem, 1, \(x) {

  activity <- x["ActivityIdentifier"]
  chems <- names(x[6:length(x)][which(!is.na(x[6:length(x)]))])
  nChems <- length(chems)
  
  return(list(chems = chems, nChems = nChems, activity = activity))
  
})

# visualize the number of sampling e
table(sapply(chem_summary, \(x) x$nChems))
hist(sapply(chem_summary, \(x) x$nChems))

# find sweet spot for balancing chems measured with number of sampling events
# looks like sampling events where >= 210 chemicals measured is good
chem_summary_reduced <- chem_summary[sapply(chem_summary, \(x) x$nChems) >= 210]

# if we look at chemicals measured at all of those sampling events there are 66 in common
# with many pesticides included in that list
# there are two oxygen analogues we will remove leaving 64 chems
shared_chems <- Reduce(intersect, sapply(chem_summary_reduced, 
                                         \(x) x$chems))[!Reduce(intersect, sapply(chem_summary_reduced, \(x) x$chems)) %in%
                                                          c("chem_Azinphos-methyl oxygen analog",
                                                            "chem_Terbufos oxygen analog sulfone")]



# after identifying the 64 shared chemicals, go back to full dataset and get all activities
# with those 64 chems because there may be instances where less than 210 were measured, but
# those 64 were still present
unique_activities <- unique(chem$ActivityIdentifier)

activities <- unique_activities[sapply(unique_activities, \(x) {
  activ <- chem[chem$ActivityIdentifier == x,]
  check_shared <- all(!is.na(activ[,..shared_chems]))
  return(check_shared)
})]


# reduce chemical data to only include sampling events with those 64 shared chemicals
# remove columns for other chems
chem_reduced <- chem %>% 
  filter(ActivityIdentifier %in% activities) %>%
  dplyr::select(do.call(c, flatten(list("ActivityIdentifier", 
                                        "CollectionDate", 
                                        "Latitude_dd", 
                                        "Longitude_dd", 
                                        "SiteNumber",
                                        shared_chems))))

# get macroinvert and fish data
bugs <- getInvertData(dataType = "density", taxonFix = "remove", agency = "USGS", rarefy = F)
fish <- getFishData(dataType = "abun", agency = "USGS", standardize = "MGMS")


# get shared sites for chem+inverts and chem+fish
commonSites_bugs <- intersect(chem_reduced$SiteNumber, bugs$SiteNumber)
commonSites_fish <- intersect(chem_reduced$SiteNumber, fish$SiteNumber)

# reduce both datasets to only include shared sites

# data for sites with both invert and chem data
chem_bugs_sharedSites <- chem_reduced %>% filter(SiteNumber %in% commonSites_bugs)
bugs_sharedSites <- bugs %>% filter(SiteNumber %in% commonSites_bugs)

# data for sites with both fish and chem data
chem_fish_sharedSites <- chem_reduced %>% filter(SiteNumber %in% commonSites_fish)
fish_sharedSites <- fish %>% filter(SiteNumber %in% commonSites_fish)

# save filtered data
data.table::fwrite(chem_bugs_sharedSites, "../data/chemData_filtered_sitesWithInvertData_011326.csv")
data.table::fwrite(bugs_sharedSites, "../data/invertData_filtered_sitesWithChemData_011326.csv")

data.table::fwrite(chem_fish_sharedSites, "../data/chemData_filtered_sitesWithFishData_011326.csv")
data.table::fwrite(fish_sharedSites, "../data/fishData_filtered_sitesWithChemData_011326.csv")

