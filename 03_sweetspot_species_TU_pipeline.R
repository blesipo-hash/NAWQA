suppressWarnings({
  library(data.table)
})

today_tag <- format(Sys.Date(), "%Y%m%d")
buff_days <- 90L
min_families_sequence <- c(15L, 10L, 5L)
top_k_sequence <- c(25L, 20L, 15L, 10L)

fish_file <- "fishData_filtered_sitesWithChemData_011326.csv"
chem_file <- "chemData_filtered_sitesWithFishData_011326.csv"
tax_file <- "FishTaxonomy_42024.csv"
key_file <- "cas_chemname_key.csv"
tox_file <- "standartox_allData_090825.csv"

out_eventmeans <- sprintf("fish_chem_eventMeans_sweetspot_%s.csv", today_tag)
out_species <- sprintf("fish_speciesWithToxicUnits_sweetspot_%s.csv", today_tag)
out_tu <- sprintf("fish_toxicUnits_sweetspot_%s.csv", today_tag)
out_qc <- sprintf("fish_TU_QC_sweetspot_%s.txt", today_tag)

standardize_site <- function(x) {
  x_chr <- trimws(as.character(x))
  out <- x_chr
  numericish <- grepl("^[0-9]+$", x_chr)
  usgs <- grepl("^USGS-", x_chr, ignore.case = TRUE)
  out[!usgs & numericish] <- sprintf("USGS-%08d", as.integer(x_chr[!usgs & numericish]))
  usgs2 <- grepl("^USGS-", out, ignore.case = TRUE)
  if (any(usgs2)) {
    suffix <- sub("^USGS-", "", toupper(out[usgs2]))
    digits <- gsub("[^0-9]", "", suffix)
    pad <- nchar(digits) > 0
    suffix[pad] <- sprintf("%08d", as.integer(digits[pad]))
    out[usgs2] <- paste0("USGS-", suffix)
  }
  out
}

convert_to_ugL <- function(x, unit) {
  u <- tolower(trimws(as.character(unit)))
  out <- rep(NA_real_, length(x))
  out[u %in% c("ug/l", "µg/l", "μg/l", "ppb")] <- x[u %in% c("ug/l", "µg/l", "μg/l", "ppb")]
  out[u %in% c("ng/l")] <- x[u %in% c("ng/l")] / 1000
  out[u %in% c("mg/l")] <- x[u %in% c("mg/l")] * 1000
  out[u %in% c("g/l")] <- x[u %in% c("g/l")] * 1e6
  out[u %in% c("g/m3", "g/m^3")] <- x[u %in% c("g/m3", "g/m^3")] * 1e3
  out
}

for (f in c(fish_file, chem_file, tax_file, key_file, tox_file)) {
  if (!file.exists(f)) {
    stop(sprintf("Required file missing: %s", f))
  }
}

fish <- fread(fish_file)
chem <- fread(chem_file)
tax <- fread(tax_file)
key <- fread(key_file)
tox <- fread(tox_file)

for (nm in c("SiteNumber", "CollectionDate")) {
  if (!nm %in% names(fish)) stop(sprintf("fish missing %s", nm))
  if (!nm %in% names(chem)) stop(sprintf("chem missing %s", nm))
}

fish[, CollectionDate := as.Date(CollectionDate)]
chem[, CollectionDate := as.Date(CollectionDate)]
fish[, SiteNumber := standardize_site(SiteNumber)]
chem[, SiteNumber := standardize_site(SiteNumber)]

events <- unique(fish[, .(SiteNumber, CollectionDate)])
chem_cols <- grep("^chem_", names(chem), value = TRUE)
if (!length(chem_cols)) stop("No chem_ columns detected in chem input")

meta_candidates <- c(
  "Agency", "SampleID", "ProjectLabel", "SiteNumber", "StudyReachName", "CollectionDate",
  "CollectionYear", "CollectionMonth", "CollectionDayOfYear", "Latitude_dd", "Longitude_dd",
  "CoordinateDatum", "COMID", "StreamOrder", "WettedWidth", "PredictedWettedWidth_m",
  "NARS_Ecoregion", "SampleTypeCode", "ReachLengthFished_m", "SampleMethod", "MethodEffort",
  "MethodEffort_units"
)
meta_cols <- intersect(names(fish), meta_candidates)
fish_species_cols <- setdiff(names(fish), meta_cols)

setkey(chem, SiteNumber, CollectionDate)
chem_event_means <- events[, {
  lo <- CollectionDate - buff_days
  hi <- CollectionDate + buff_days
  sub <- chem[.(SiteNumber)][CollectionDate >= lo & CollectionDate <= hi]
  if (nrow(sub) == 0) {
    NULL
  } else {
    ans <- sub[, lapply(.SD, function(x) mean(as.numeric(x), na.rm = TRUE)), .SDcols = chem_cols]
    ans[, n_chem_samples_in_window := nrow(sub)]
    ans
  }
}, by = .(SiteNumber, CollectionDate)]

fwrite(chem_event_means, out_eventmeans)

matched_events <- chem_event_means[, .(SiteNumber, CollectionDate)]
setkey(matched_events, SiteNumber, CollectionDate)
setkey(fish, SiteNumber, CollectionDate)
fish_final <- fish[matched_events, nomatch = 0]

for (cc in fish_species_cols) {
  v <- as.numeric(fish_final[[cc]])
  fin <- is.finite(v)
  if (any(!fin) && any(fin)) {
    v[!fin] <- mean(v[fin], na.rm = TRUE)
  }
  fish_final[[cc]] <- v
}

needed_tox <- c("qualifier", "endpoint", "duration_unit", "duration", "concentration", "concentration_unit")
miss_tox <- setdiff(needed_tox, names(tox))
if (length(miss_tox)) stop(sprintf("tox missing columns: %s", paste(miss_tox, collapse = ", ")))

tox[, duration_num := as.numeric(duration)]
tox <- tox[
  qualifier == "=" & endpoint %in% c("LC50", "LC50*") &
    duration_unit == "h" & duration_num == 96
]
tox[, cas_use := if ("cas" %in% names(tox)) as.character(cas) else as.character(casnr)]
tox[, concentration_std_ugL := convert_to_ugL(as.numeric(concentration), concentration_unit)]
tox <- tox[is.finite(concentration_std_ugL) & concentration_std_ugL > 0 & !is.na(cas_use) & nzchar(cas_use)]

if (!"family" %in% names(tox)) stop("tox file must include family column for family-level LC50 aggregation")
tox[, family := trimws(as.character(family))]

key[, CharacteristicName := trimws(as.character(CharacteristicName))]
key[, cas := trimws(as.character(cas))]
chem_map <- data.table(chem_col = chem_cols, CharacteristicName = sub("^chem_", "", chem_cols))
chem_map <- merge(chem_map, unique(key[, .(CharacteristicName, cas)]), by = "CharacteristicName", all.x = TRUE)
chem_map <- chem_map[!is.na(cas) & nzchar(cas)]

lc50_fam <- tox[, .(lc50_fam_ugL = mean(concentration_std_ugL, na.rm = TRUE)), by = .(cas = cas_use, family)]
coverage <- lc50_fam[, .(n_families = uniqueN(family)), by = cas]
chem_map_cov <- merge(chem_map, coverage, by = "cas", all.x = TRUE)
chem_map_cov[is.na(n_families), n_families := 0L]

selection_success <- FALSE
selected <- list()

for (minf in min_families_sequence) {
  cand <- chem_map_cov[n_families >= minf]
  if (!nrow(cand)) next
  setorder(cand, -n_families, chem_col)

  for (k in top_k_sequence) {
    cand_k <- cand[seq_len(min(k, .N))]
    fam_count <- lc50_fam[cas %in% cand_k$cas, .N, by = family]
    families_final <- fam_count[N == nrow(cand_k), family]

    if (length(families_final) >= 3) {
      selection_success <- TRUE
      selected <- list(minf = minf, k = k, cand = cand_k, families = families_final)
      break
    }
  }
  if (selection_success) break
}

qc <- c(
  sprintf("date_tag: %s", today_tag),
  sprintf("events_total: %d", nrow(events)),
  sprintf("events_matched_with_chem_window: %d", nrow(chem_event_means)),
  sprintf("chem_cols_total: %d", length(chem_cols)),
  sprintf("chem_cols_mapped_to_cas: %d", nrow(chem_map))
)

if (!selection_success) {
  qc <- c(qc,
          "coverage_selection: FAILED",
          "No combination of min_families_per_chem (15/10/5) and top_k_chems (25/20/15/10) yielded >=3 families with full chemical coverage.")
  writeLines(qc, out_qc)
  stop("Coverage-aware selection failed; see QC report.")
}

selected_chems <- selected$cand
families_final <- selected$families

# Species -> Family mapping
tax[, Species_period := if ("Species_period" %in% names(tax)) as.character(Species_period) else gsub(" ", ".", as.character(Species))]
tax[, Species_norm := gsub("\\.", " ", Species_period)]
tax[, Family := as.character(Family)]

species_lookup <- data.table(
  species_col = fish_species_cols,
  Species_norm = gsub("\\.", " ", fish_species_cols)
)
species_lookup <- merge(species_lookup, unique(tax[, .(Species_norm, Family)]), by = "Species_norm", all.x = TRUE)
species_lookup <- species_lookup[Family %in% families_final]

species_keep <- intersect(species_lookup$species_col, fish_species_cols)
if (!length(species_keep)) {
  qc <- c(qc,
          sprintf("coverage_selection: SUCCESS (min_families=%d, top_k=%d)", selected$minf, selected$k),
          sprintf("selected_chems: %d", nrow(selected_chems)),
          sprintf("families_final: %d", length(families_final)),
          "species_retained: 0 (no fish species mapped to retained families)")
  writeLines(qc, out_qc)
  stop("No fish species retained after taxonomy + family coverage filters.")
}

# Event chem matrix for selected chemicals
chem_selected_cols <- selected_chems$chem_col
chem_event_sub <- merge(matched_events, chem_event_means[, c("SiteNumber", "CollectionDate", chem_selected_cols), with = FALSE],
                        by = c("SiteNumber", "CollectionDate"), all.x = TRUE)

# family x cas LC50 lookup for selected chems
lc50_lookup <- lc50_fam[cas %in% selected_chems$cas & family %in% families_final]
setkey(lc50_lookup, family, cas)

species_family <- species_lookup[species_col %in% species_keep, .(species_col, family = Family)]
species_family <- unique(species_family)

n_evt <- nrow(chem_event_sub)
n_sp <- nrow(species_family)

tu_mat <- matrix(NA_real_, nrow = n_evt, ncol = n_sp)
nc_mat <- matrix(0L, nrow = n_evt, ncol = n_sp)

for (j in seq_len(n_sp)) {
  fam_j <- species_family$family[j]
  sp_j <- species_family$species_col[j]
  cas_j <- selected_chems$cas
  lc50_j <- lc50_lookup[list(fam_j, cas_j), lc50_fam_ugL]

  rq_cols <- selected_chems$chem_col
  rq_vals <- sapply(seq_along(rq_cols), function(i) {
    env <- as.numeric(chem_event_sub[[rq_cols[i]]])
    lc <- lc50_j[i]
    if (!is.finite(lc) || lc <= 0) {
      rep(NA_real_, length(env))
    } else {
      env / lc
    }
  })

  if (is.vector(rq_vals)) rq_vals <- matrix(rq_vals, ncol = 1)
  tu_mat[, j] <- rowSums(rq_vals, na.rm = TRUE)
  tu_mat[!is.finite(tu_mat[, j]), j] <- NA_real_
  nc_mat[, j] <- rowSums(is.finite(rq_vals), na.rm = TRUE)
}

colnames(tu_mat) <- species_family$species_col
colnames(nc_mat) <- paste0(species_family$species_col, "__nchems")

# Save RQ array as list of event matrices by species x chemical
rq_array <- array(NA_real_, dim = c(n_sp, length(chem_selected_cols), n_evt),
                  dimnames = list(species_family$species_col, chem_selected_cols,
                                  paste(chem_event_sub$SiteNumber, chem_event_sub$CollectionDate, sep = "|")))
for (j in seq_len(n_sp)) {
  fam_j <- species_family$family[j]
  cas_j <- selected_chems$cas
  lc50_j <- lc50_lookup[list(fam_j, cas_j), lc50_fam_ugL]
  for (i in seq_along(chem_selected_cols)) {
    env <- as.numeric(chem_event_sub[[chem_selected_cols[i]]])
    lc <- lc50_j[i]
    rq_array[j, i, ] <- if (is.finite(lc) && lc > 0) env / lc else NA_real_
  }
}

rq_qs_file <- sprintf("fish_riskQuotientArray_sweetspot_%s.qs", today_tag)
rq_rds_file <- sprintf("fish_riskQuotientArray_sweetspot_%s.rds", today_tag)
if (requireNamespace("qs", quietly = TRUE)) {
  qs::qsave(rq_array, rq_qs_file)
  rq_saved <- rq_qs_file
} else {
  saveRDS(rq_array, rq_rds_file)
  rq_saved <- rq_rds_file
}

fish_species_used <- fish_final[, c(meta_cols, species_keep), with = FALSE]
fwrite(fish_species_used, out_species)

tu_dt <- data.table(SiteNumber = chem_event_sub$SiteNumber, CollectionDate = chem_event_sub$CollectionDate)
tu_dt <- cbind(tu_dt, as.data.table(tu_mat))
fwrite(tu_dt, out_tu)

qc <- c(
  qc,
  sprintf("coverage_selection: SUCCESS (min_families=%d, top_k=%d)", selected$minf, selected$k),
  sprintf("selected_chems: %d", length(chem_selected_cols)),
  sprintf("families_final: %d", length(families_final)),
  sprintf("species_total_in_fish: %d", length(fish_species_cols)),
  sprintf("species_retained: %d", length(species_keep)),
  sprintf("pct_species_retained: %.2f", 100 * length(species_keep) / max(1, length(fish_species_cols))),
  sprintf("rq_array_saved: %s", rq_saved),
  "TU_summary_all_species_values:"
)
qc <- c(qc, capture.output(print(summary(as.vector(tu_mat)))))
writeLines(qc, out_qc)

cat(sprintf("Wrote: %s\n", out_eventmeans))
cat(sprintf("Wrote: %s\n", out_species))
cat(sprintf("Wrote: %s\n", out_tu))
cat(sprintf("Wrote: %s\n", rq_saved))
cat(sprintf("Wrote: %s\n", out_qc))
