library(data.table)
library(tidyverse)

##### TEMPLATE: prepare chemical and macroinvertebrate data for risk quotient calculations #####
##### This script is adapted from `riskQuotient_TU_calcs_fish.R` and includes placeholders. #####

# --- INPUT FILES (PLACEHOLDER: update paths/filenames as needed) ---
invert_file <- "../data/invertData_filtered_sitesWithChemData_011326.csv"
chem_file <- "../data/chemData_filtered_sitesWithInvertData_011326.csv"
tox_file <- "../data/standartox_allData_090825.csv"

# PLACEHOLDER: invert taxonomy crosswalk
# Provide a file/table that links your invert abundance column names to taxonomic ranks
# needed for toxicity joins (e.g., Species -> Genus/Family/Order or TaxonCode -> Family).
# Example expected columns might include:
#   TaxonLabel, Genus, Family, Order, Class
invert_taxonomy_file <- "../data/INVERT_TAXONOMY_CROSSWALK_PLACEHOLDER.csv"

# --- OUTPUT FILES (PLACEHOLDER: update naming conventions as desired) ---
out_rq_array <- "invert_riskQuotientArray_TEMPLATE.qs"
out_tu_csv <- "../data/invert_toxicUnits_TEMPLATE.csv"
out_invert_with_tu_csv <- "../data/invert_taxaWithToxicUnits_TEMPLATE.csv"

# read in invert density and chemical data
invert <- fread(invert_file) %>%
  mutate(CollectionDate = as_date(CollectionDate))
chem <- fread(chem_file) %>%
  mutate(CollectionDate = as_date(CollectionDate))

# aggregate chemistry by invert sampling event:
# for each invert sampling event, average chemistry within +/- buff days
buff <- 90

chem_list_invert <- list()

for (i in 1:nrow(invert)) {

  # site of sampling event
  site <- invert[i, ]$SiteNumber

  # event date and date window
  sampDate <- invert[i, ]$CollectionDate
  minDate <- sampDate - buff
  maxDate <- sampDate + buff

  # all chemistry data in window and site
  chemDat <- chem %>%
    filter(CollectionDate >= minDate,
           CollectionDate <= maxDate,
           SiteNumber == site)

  activities <- unique(chemDat$ActivityIdentifier)

  chemDat_avg <- NA

  if (nrow(chemDat) > 0) {

    chemDat_avg <- as.data.frame(matrix(colMeans(chemDat[, 6:ncol(chemDat)]),
                                        ncol = ncol(chemDat) - 5))
    colnames(chemDat_avg) <- colnames(chemDat)[6:ncol(chemDat)]

    chemDat_avg$SiteNumber <- site
    chemDat_avg$CollectionDate <- sampDate
    chemDat_avg$ActivityIdentifiers <- str_flatten(activities, collapse = " | ")

    chemDat_avg <- chemDat_avg %>%
      relocate(any_of(c("SiteNumber", "CollectionDate", "ActivityIdentifiers")),
               .before = starts_with("chem_"))
  }

  chem_list_invert[[i]] <- list(
    site = site,
    sampDate = sampDate,
    activities = activities,
    samplingEvents_chem = nrow(chemDat),
    chemDat_avg = chemDat_avg
  )
}

# number of invert events with chemistry within buff days
length(which(sapply(chem_list_invert, function(x) all(!is.na(x$chemDat_avg)))))

# remove events with no chemistry match
chem_list_invert_abbrev <- chem_list_invert[sapply(chem_list_invert, function(x) all(!is.na(x$chemDat_avg)))]

# average number of chemistry events per invert event
mean(na.omit(sapply(chem_list_invert_abbrev, function(x) x$samplingEvents_chem)))

# filter chemistry and invert data to matched events
chem_final <- do.call(rbind, lapply(chem_list_invert_abbrev, function(x) x$chemDat_avg))
invert_final <- invert %>%
  filter(paste(CollectionDate, SiteNumber) %in%
           paste(chem_final$CollectionDate, chem_final$SiteNumber))

# ----------------------------------------------------------------------------
# PLACEHOLDER: identify invert abundance columns
# ----------------------------------------------------------------------------
# The fish script uses fixed abundance-column start positions.
# For macroinvertebrates, replace the line below with your real column selection.
# Example options:
#   - abundance_cols <- 23:ncol(invert_final)
#   - abundance_cols <- grep("^Taxon_", colnames(invert_final))
#   - abundance_cols <- which(colnames(invert_final) %in% target_taxa)
abundance_cols <- 23:ncol(invert_final)

# optional cleanup for Inf values in abundance columns
if (length(abundance_cols) > 0) {
  invert_final[, abundance_cols] <- apply(invert_final[, abundance_cols], 2, function(x) {
    if (any(is.infinite(x))) {
      meanConc <- mean(x[!is.infinite(x)], na.rm = TRUE)
      x[is.infinite(x)] <- meanConc
      return(x)
    }
    x
  }) %>% as.data.frame()
}

#### Join toxicity information

# taxon labels from abundance columns
invert_taxa <- colnames(invert_final)[abundance_cols]

# load taxonomy crosswalk (PLACEHOLDER: ensure required columns are present)
invert_tax <- fread(invert_taxonomy_file)

# load toxicity data
tox_data <- fread(tox_file)

# ----------------------------------------------------------------------------
# PLACEHOLDER: chosen toxicity endpoint/duration filters for invertebrates
# ----------------------------------------------------------------------------
# Update these based on your invert toxicity decision.
# Examples:
#   endpoint %in% c("LC50", "EC50")
#   duration == "48" or duration == "96"
#   group %in% c("insect", "invertebrate")
tox_filtered <- tox_data %>%
  filter(
    # PLACEHOLDER: invert grouping (e.g., group %in% c("insect", "invertebrate"))
    group %in% c("insect", "invertebrate"),
    # PLACEHOLDER: endpoint choice
    endpoint %in% c("LC50"),
    # PLACEHOLDER: duration choice
    duration_unit == "h",
    duration == "96",
    qualifier == "=",
    concentration_unit %in% c("g/l", "ppb", "g/m3")
  )

# convert toxicity units to ug/L
tox_filtered$concentration_std <- NA

tox_filtered$concentration_std[tox_filtered$concentration_unit == "g/l"] <-
  tox_filtered$concentration[tox_filtered$concentration_unit == "g/l"] * 1e6
tox_filtered$concentration_std[tox_filtered$concentration_unit == "ppb"] <-
  tox_filtered$concentration[tox_filtered$concentration_unit == "ppb"]
tox_filtered$concentration_std[tox_filtered$concentration_unit == "g/m3"] <-
  tox_filtered$concentration[tox_filtered$concentration_unit == "g/m3"] * 1e3

# ----------------------------------------------------------------------------
# PLACEHOLDER: choose taxonomic rank for toxicity aggregation
# ----------------------------------------------------------------------------
# For fish, toxicity is aggregated by family.
# For inverts, choose one rank (e.g., Family, Genus, Order) with best coverage.
agg_rank_col <- "family"  # PLACEHOLDER

# PLACEHOLDER: chemical key mapping your NAWQA chemical names to CAS
# Must produce at least: cas, CharacteristicName
cas_key_file <- "../data/cas_chemname_key.csv"
if (!file.exists(cas_key_file)) {
  stop("Missing CAS key file: ", cas_key_file,
       "\nPLACEHOLDER: provide a chemical key with `cas` and `CharacteristicName` columns.")
}
cas_key <- fread(cas_key_file)

if (!all(c("cas", "CharacteristicName") %in% colnames(cas_key))) {
  stop("CAS key must include columns: cas, CharacteristicName")
}

# keep only chemicals present in chemistry table
chem_names <- str_replace(colnames(chem_final)[4:ncol(chem_final)], "^chem_", "")
cas_key_use <- cas_key %>%
  filter(CharacteristicName %in% chem_names)

# aggregate toxicity by selected invert rank and CAS
tox_agg <- tox_filtered %>%
  filter(cas %in% cas_key_use$cas) %>%
  group_by(cas, .data[[agg_rank_col]]) %>%
  reframe(tox_rank_mean = mean(concentration_std, na.rm = TRUE)) %>%
  left_join(cas_key_use %>% dplyr::select(cas, CharacteristicName) %>% distinct(), by = "cas")

# ----------------------------------------------------------------------------
# PLACEHOLDER: map invert taxa to chosen aggregation rank
# ----------------------------------------------------------------------------
# crosswalk must include:
#   TaxonLabel (matching abundance column names) and selected rank column (agg_rank_col)
if (!"TaxonLabel" %in% colnames(invert_tax)) {
  stop("Invert taxonomy crosswalk must include `TaxonLabel`.")
}
if (!agg_rank_col %in% colnames(invert_tax)) {
  stop("Invert taxonomy crosswalk must include chosen rank column: ", agg_rank_col)
}

invert_tax_use <- invert_tax %>%
  filter(TaxonLabel %in% invert_taxa)

# taxa with toxicity coverage for all selected chemicals (optional strict filter)
tox_table <- table(tox_agg$CharacteristicName, tox_agg[[agg_rank_col]])
ranks_with_full_coverage <- colnames(tox_table)[apply(tox_table, 2, function(x) all(x >= 1))]

invert_taxa_final <- invert_tax_use$TaxonLabel[invert_tax_use[[agg_rank_col]] %in% ranks_with_full_coverage]

# filter invert data to taxa with toxicity coverage
invert_dat_tox <- invert_final %>%
  dplyr::select(colnames(invert_final)[1:(min(abundance_cols) - 1)], invert_taxa_final)

##### calculate risk quotients for each invert taxon based on chosen rank-level toxicity #####

invert_rq_dat <- array(
  data = NA,
  dim = c(length(invert_taxa_final),
          ncol(chem_final) - 3,
          nrow(invert_dat_tox)),
  dimnames = list(
    Taxon = invert_taxa_final,
    Chemical = colnames(chem_final)[4:ncol(chem_final)],
    SamplingEvent = paste(invert_dat_tox$SiteNumber, invert_dat_tox$CollectionDate, sep = "|")
  )
)

for (i in 1:length(invert_rq_dat[1, , 1])) {

  ch <- dimnames(invert_rq_dat)$Chemical[i]
  ch_name <- str_replace(ch, "^chem_", "")

  for (j in 1:length(invert_rq_dat[, 1, 1])) {

    tx <- dimnames(invert_rq_dat)$Taxon[j]
    rank_val <- invert_tax_use[[agg_rank_col]][invert_tax_use$TaxonLabel == tx][1]

    tox_val <- tox_agg$tox_rank_mean[
      tox_agg[[agg_rank_col]] == rank_val & tox_agg$CharacteristicName == ch_name
    ]

    if (length(tox_val) == 0 || is.na(tox_val[1])) {
      next
    }

    ch_RQs <- chem_final[, ch] / tox_val[1]
    invert_rq_dat[j, i, ] <- ch_RQs
  }
}

# sum RQs across chemicals to get TUs by taxon and sampling event
invert_tu_dat <- t(apply(invert_rq_dat, c(1, 3), sum, na.rm = TRUE)) %>% as.data.frame()

# rejoin event info
invert_tu_dat$SiteNumber <- sapply(rownames(invert_tu_dat), function(x) str_split_1(x, "\\|")[1])
invert_tu_dat$CollectionDate <- sapply(rownames(invert_tu_dat), function(x) str_split_1(x, "\\|")[2])

invert_tu_dat2 <- invert_tu_dat %>%
  relocate(any_of(c("SiteNumber", "CollectionDate")))

# save outputs (PLACEHOLDER filenames above)
qs::qsave(invert_rq_dat, out_rq_array)
fwrite(invert_tu_dat2, out_tu_csv)
fwrite(invert_dat_tox, out_invert_with_tu_csv)
