devtools::install_github("USEPA/finsyncR",
                         build_vignette = TRUE)

library(tidyverse)
library(parallel)

chemData <- data.table::fread("../data/NAWQAdata_052824_covariateColumnsAdded.csv")

# handle units and detection limit: standardize units and assume non-detects = 1/2 detection limit

# convert all units to ug/L
chemData$ResultStd_ugL <- NA
chemData$DetectionStd_ugL <- NA

chemData$ResultStd_ugL[chemData$ResultMeasure.MeasureUnitCode == "ug/l"] <- chemData$ResultMeasureValue[chemData$ResultMeasure.MeasureUnitCode == "ug/l"]
chemData$ResultStd_ugL[chemData$ResultMeasure.MeasureUnitCode == "ng/l"] <- chemData$ResultMeasureValue[chemData$ResultMeasure.MeasureUnitCode == "ng/l"]/1000
chemData$ResultStd_ugL[chemData$ResultMeasure.MeasureUnitCode == "mg/l"] <- chemData$ResultMeasureValue[chemData$ResultMeasure.MeasureUnitCode == "mg/l"]*1000

chemData$DetectionStd_ugL[chemData$DetectionQuantitationLimitMeasure.MeasureUnitCode == "ug/l"] <- chemData$DetectionQuantitationLimitMeasure.MeasureValue[chemData$DetectionQuantitationLimitMeasure.MeasureUnitCode == "ug/l"]
chemData$DetectionStd_ugL[chemData$DetectionQuantitationLimitMeasure.MeasureUnitCode == "ng/l"] <- chemData$DetectionQuantitationLimitMeasure.MeasureValue[chemData$DetectionQuantitationLimitMeasure.MeasureUnitCode == "ng/l"]/1000
chemData$DetectionStd_ugL[chemData$DetectionQuantitationLimitMeasure.MeasureUnitCode == "mg/l"] <- chemData$DetectionQuantitationLimitMeasure.MeasureValue[chemData$DetectionQuantitationLimitMeasure.MeasureUnitCode == "mg/l"]*1000


# handle all detection scenarios: 
# "Not Detected"                       ""                                   "Detected Not Quantified"           
# "Systematic Contamination"           "Present Above Quantification Limit"

# set non-detects to "1/2 detection"
chemData$ResultStd_ugL[chemData$ResultDetectionConditionText == "Not Detected" &
                         chemData$DetectionQuantitationLimitMeasure.MeasureUnitCode %in% c("ug/l","ng/l","mg/l")] <- chemData$DetectionStd_ugL[chemData$ResultDetectionConditionText == "Not Detected" &
                                                                                                                                                 chemData$DetectionQuantitationLimitMeasure.MeasureUnitCode %in% c("ug/l","ng/l","mg/l")]/2
# for samples labeled "Systematic Contamination" set result to deviation from minimum concentration (assumed to be equal to baseline contaminant concentration)                                                          

contaminatedChems <- unique(chemData$CharacteristicName[chemData$ResultDetectionConditionText == "Systematic Contamination"])

for (i in seq_along(contaminatedChems)) {
  
  chem <- contaminatedChems[i]
  dat <- chemData[chemData$ResultDetectionConditionText == "Systematic Contamination" & chemData$CharacteristicName == chem,]
  
  # ensure the ResultStd starts from the uncorrected value (in case code is run multiple times)
  dat$ResultStd_ugL[dat$ResultMeasure.MeasureUnitCode == "ug/l"] <- dat$ResultMeasureValue[dat$ResultMeasure.MeasureUnitCode == "ug/l"]
  dat$ResultStd_ugL[dat$ResultMeasure.MeasureUnitCode == "ng/l"] <- dat$ResultMeasureValue[dat$ResultMeasure.MeasureUnitCode == "ng/l"]/1000
  dat$ResultStd_ugL[dat$ResultMeasure.MeasureUnitCode == "mg/l"] <- dat$ResultMeasureValue[dat$ResultMeasure.MeasureUnitCode == "mg/l"]*1000
  
  bcl <- min(dat$ResultStd_ugL)
  
  correctedResult <- dat$ResultStd_ugL - bcl
  
  chemData$ResultStd_ugL[chemData$ResultDetectionConditionText == "Systematic Contamination" & chemData$CharacteristicName == chem] <- correctedResult
  
}

# for results labeled "Present Above Quantification Limit" set result equal to the upper reporting limit
chemData$ResultStd_ugL[chemData$ResultDetectionConditionText == "Present Above Quantification Limit"] <- chemData$DetectionStd_ugL[chemData$ResultDetectionConditionText == "Present Above Quantification Limit"]

# for results labeled "Detected Not Quantified" set result equal to the detection limit
chemData$ResultStd_ugL[chemData$ResultDetectionConditionText == "Detected Not Quantified"] <- chemData$DetectionStd_ugL[chemData$ResultDetectionConditionText == "Detected Not Quantified"]

# 2 non-detect rows that show NAs because no detection limit is given. setting these to 0
chemData$ResultStd_ugL[is.na(chemData$ResultStd_ugL) & (chemData$ResultMeasure.MeasureUnitCode %in% c("ug/l","ng/l","mg/l") | chemData$DetectionQuantitationLimitMeasure.MeasureUnitCode %in% c("ug/l","ng/l","mg/l"))] <- 0


# get only necessary columns
chemData_abbrev <- chemData[chemData$ResultMeasure.MeasureUnitCode %in% c("ug/l","ng/l","mg/l") | chemData$DetectionQuantitationLimitMeasure.MeasureUnitCode %in% c("ug/l","ng/l","mg/l"),] %>%
  dplyr::select(ActivityIdentifier, ActivityStartDate, MonitoringLocationIdentifier, LatitudeMeasure, LongitudeMeasure, 
                USGSPCode, chemical_group, CharacteristicName, parm_nm, ResultStd_ugL, ResultDetectionConditionText) %>%
  mutate(logResult = log(ResultStd_ugL),
         Filtered = ifelse(str_detect(parm_nm, " filtered"), TRUE, FALSE),
         Micrograms = ifelse(str_detect(parm_nm, "micrograms"), TRUE, FALSE)) %>%
  dplyr::filter(Filtered) %>%
  group_by(ActivityIdentifier) %>%
  mutate(nChems = n_distinct(USGSPCode)) %>% ungroup()

# parse out redundant measurements, e.g. unfiltered vs. filtered samples taken at the same time/place
# and measurements that were repeated in both nanograms and micrograms
table(chemData_abbrev$Filtered)
table(chemData_abbrev$Micrograms)

characteristics <- unique(chemData_abbrev$CharacteristicName)

char_summary <- data.frame(CharacteristicName = characteristics, DetailedNames = "", n_tot = NA, filtered_nm = "", n_filtered = 0, unfiltered_nm = "", n_unfiltered = 0)

for (i in seq_along(characteristics)) {
  
  char_nm <- characteristics[i]
  
  parms <- unique(chemData_abbrev$parm_nm[chemData_abbrev$CharacteristicName == char_nm])
  filtered <- parms[sapply(parms, \(x) str_detect(x, " filtered"))]
  unfiltered <- parms[sapply(parms, \(x) str_detect(x, "unfiltered"))]
  
  char_summary$DetailedNames[char_summary$CharacteristicName == char_nm] <- str_flatten(parms, collapse = " || ")
  char_summary$filtered_nm[char_summary$CharacteristicName == char_nm] <- str_flatten(filtered, collapse = " || ")
  char_summary$unfiltered_nm[char_summary$CharacteristicName == char_nm] <- str_flatten(unfiltered, collapse = " || ")
  
  char_summary$n_tot[char_summary$CharacteristicName == char_nm] <- length(parms)
  char_summary$n_filtered[char_summary$CharacteristicName == char_nm] <- length(filtered)
  char_summary$n_unfiltered[char_summary$CharacteristicName == char_nm] <- length(unfiltered)
  
}

# check for redundant filtered/unfiltered measurements (i.e. the same sample reported as both filtered and unfiltered)
# and/or whether the same sample was reported in both micrograms and nanograms
microgramChems <- unique(chemData_abbrev$CharacteristicName[chemData_abbrev$Micrograms])

# IMPORTANT!!! Disable OpenBLAS multithreading before starting parallel engine
RhpcBLASctl::blas_set_num_threads(1)

if (RhpcBLASctl::blas_get_num_procs() == 1) {
  # specify number of threads to leave free
  free_threads <- 6
  
  # set up parallel computing
  n.cores <- parallel::detectCores() - free_threads
  
  my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
  
  clusterExport(cl = my.cluster, varlist = c("chemData_abbrev", "microgramChems"))

  # run parallel function
  unitCheck <- parLapply(X = microgramChems, fun = \(x) {
    
    dat_index <- which(chemData_abbrev$CharacteristicName == x)
    dat <- chemData_abbrev[chemData_abbrev$CharacteristicName == x,]
    
    samplingEvents <- unique(dat$ActivityIdentifier)
    
    checks <- lapply(samplingEvents,\(y) {
      # check whether nanogram and microgram values were reported for the same sample
      NanoAndMicro <- any(dat$Micrograms[dat$ActivityIdentifier == y]) & any(!dat$Micrograms[dat$ActivityIdentifier == y])
      res <- NanoAndMicro
      
      return(res)
    })
    
    print(paste("Done:", x))
    
    res <- samplingEvents[sapply(checks, \(y) y)]
    
    return(res)
    
  },
  cl = my.cluster)
  
  # always stop the cluster when you are done using it
  parallel::stopCluster(cl = my.cluster)

}

names(unitCheck) <- microgramChems

unitCheck_reduced <- unitCheck[sapply(unitCheck, \(x) length(x) >0)]

unitCheck_selections <- list()

chemData_final <- chemData_abbrev

for (i in seq_along(unitCheck_reduced)) {
  
  char_nm <- names(unitCheck_reduced)[i]
  
  activities <- unitCheck_reduced[[i]]
  
  selectedParms <- c()
  
  for (j in seq_along(activities)) {
    
    parm_nms <- unique(chemData_abbrev$parm_nm[chemData_abbrev$CharacteristicName == char_nm & chemData_abbrev$ActivityIdentifier == activities[j]])
    
    selectedParm <- menu(parm_nms, title = "which parameter should be kept? (the rest will be discarded)")
    
    selectedParms[j] <- parm_nms[selectedParm]
    
    chemData_final <- chemData_final %>%
      dplyr::filter(!(CharacteristicName == char_nm & ActivityIdentifier == activities[j] & parm_nm != selectedParms[j]))
    
  }
  
  unitCheck_selections[[i]] <- selectedParms
  
}

names(unitCheck_selections) <- names(unitCheck_reduced)

# manually check the remaining parameter names for potential removal

unique(chemData_final$parm_nm)

rm_parms <- c("Alachlor, water, filtered (0.7 micron glass fiber filter), enzyme-linked immunosorbent assay, recoverable, micrograms per liter",
              "Glyphosate, water, filtered, immunoassay, unadjusted, recoverable, micrograms per liter",
              "Total anatoxin-a, water, unfiltered, freeze/thaw and then filtered (0.7 micron filter), enzyme-linked immunosorbent assay, recoverable, micrograms per liter",
              "Atrazine, water, filtered, immunoassay, unadjusted, recoverable, micrograms per liter",
              "Cyanazine, water, filtered (0.7 micron glass fiber filter), enzyme-linked immunosorbent assay, recoverable, micrograms per liter",
              "2,4-D, water, filtered (0.7 micron glass fiber filter), enzyme-linked immunosorbent assay, recoverable, micrograms per liter",
              "Metolachlor, water, filtered (0.7 micron glass fiber filter), enzyme-linked immunosorbent assay, recoverable, micrograms per liter")

# final check for duplicate parameter names for a single chemical for a single sampling event (will check these manually and decide which to keep)
chemData_final2 <- chemData_final %>%
  dplyr::filter(!parm_nm %in% rm_parms) %>%
  group_by(CharacteristicName, ActivityIdentifier) %>%
  mutate(nDups = n_distinct(parm_nm)) %>%
  ungroup()

dupList <- unique(paste(chemData_final2$CharacteristicName[chemData_final2$nDups == 2], 
                        chemData_final2$ActivityIdentifier[chemData_final2$nDups == 2], 
                        sep = "||"))

dupParms <- unique(chemData_final2$parm_nm[chemData_final2$nDups == 2])

# only 3 chems have duplicates, remove them manually
chemData_final3 <- chemData_final2 %>%
  dplyr::filter(!(nDups == 2 & 
                    parm_nm %in% c("Terbacil, water, filtered (0.7 micron glass fiber filter), recoverable, micrograms per liter",
                                   "Methyl parathion, water, filtered (0.7 micron glass fiber filter), recoverable, micrograms per liter",
                                   "2-Methyl-4,6-dinitrophenol, water, filtered (0.7 micron glass fiber filter), recoverable, micrograms per liter"))) %>%
  dplyr::select(-nDups, -nChems, -Micrograms, -Filtered) %>%
  group_by(CharacteristicName, ActivityIdentifier) %>% slice(1) %>% ungroup()

# some additional duplicates identified by the pivot_wider() function, just using slice() to remove
# dupCheck <- chemData_final3 |>
#   dplyr::summarise(n = dplyr::n(), .by = c(ActivityIdentifier, CharacteristicName)) |>
#   dplyr::filter(n > 1L) 
#
# unique(dupCheck$CharacteristicName)

# data.table::fwrite(chemData_final3, "../data/chemData_long_071025.csv")

# chemData_final3 <- data.table::fread("../data/chemData_long_071025.csv")

# data frame to join lat long site and date data by activity identifier
activityDat <- chemData_final3 %>%
  dplyr::select(ActivityIdentifier, LatitudeMeasure, LongitudeMeasure, MonitoringLocationIdentifier, ActivityStartDate) %>%
  group_by(ActivityIdentifier) %>% slice(1) %>% ungroup()

# data to wide form
chemData_wide <- chemData_final3 %>%
  dplyr::select(-parm_nm, -ResultDetectionConditionText, -logResult, -chemical_group, -USGSPCode) %>%
  pivot_wider(names_from = CharacteristicName, 
              names_prefix = "chem_",
              values_from = ResultStd_ugL,
              id_cols = ActivityIdentifier,
              values_fill = NA) %>%
  left_join(activityDat, by = "ActivityIdentifier") %>%
  relocate(any_of(c("ActivityIdentifier", 
                    "ActivityStartDate", 
                    "LatitudeMeasure", 
                    "LongitudeMeasure", 
                    "MonitoringLocationIdentifier")), 
           .before = starts_with("chem_"))

# data.table::fwrite(chemData_wide, "../data/chemData_wide_071025.csv")
  

