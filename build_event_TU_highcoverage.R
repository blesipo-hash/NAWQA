#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

today_tag <- format(Sys.Date(), "%Y%m%d")
BUFF_DAYS <- 90L
DENSE_MIN_CHEMS <- 150L
CHEM_COVERAGE_MIN <- 0.85
EVENT_MIN_CHEMS_USED <- 20L
EXCLUDE_CHEMS <- c("chem_Azinphos-methyl oxygen analog", "chem_Terbufos oxygen analog sulfone")
EPS_TU <- 1e-8

event_size_file <- "output/data/fish_event_size_metrics_20260226.csv"
species_event_file <- "output/data/fish_species_event_size_metrics_20260226.csv"
raw_chem_candidates <- c("data/chemData_wide_071025.csv", "../data/chemData_wide_071025.csv", "chemData_wide_071025.csv")
stdtox_candidates <- c("data/standartox_allData_090825.csv", "../data/standartox_allData_090825.csv")
cas_candidates <- c("data/cas_chemname_key.csv", "../data/cas_chemname_key.csv", "cas_chemname_key.csv")

find_one <- function(paths, label) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) {
    stop(sprintf("No %s file found. Checked: %s", label, paste(paths, collapse = ", ")), call. = FALSE)
  }
  hit[1]
}

as_num <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  out <- suppressWarnings(as.numeric(x))
  out[!is.finite(out)] <- NA_real_
  out
}

safe_name <- function(x) {
  tolower(gsub("[^a-z0-9]+", "", x))
}

standardize_site <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "na", "nan")] <- NA_character_
  digits <- gsub("[^0-9]", "", x)
  digits[digits == ""] <- NA_character_
  out <- rep(NA_character_, length(digits))
  ok <- !is.na(digits)
  if (any(ok)) {
    d <- digits[ok]
    d <- ifelse(nchar(d) >= 8, substr(d, nchar(d) - 7, nchar(d)), sprintf("%08s", d))
    d <- gsub(" ", "0", d, fixed = TRUE)
    out[ok] <- paste0("USGS-", d)
  }
  out
}

must_have_cols <- function(dt, cols, label) {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) {
    stop(sprintf("%s missing required columns: %s", label, paste(miss, collapse = ", ")), call. = FALSE)
  }
}

for (f in c(event_size_file, species_event_file)) {
  if (!file.exists(f)) stop(sprintf("Required file not found: %s", f), call. = FALSE)
}

event_size_dt <- fread(event_size_file)
species_event_dt <- fread(species_event_file)

must_have_cols(event_size_dt, c("SiteNumber", "CollectionDate"), "fish_event_size_metrics")
must_have_cols(species_event_dt, c("SiteNumber", "CollectionDate"), "fish_species_event_size_metrics")

event_size_dt[, SiteNumber := standardize_site(SiteNumber)]
event_size_dt[, CollectionDate := as.Date(CollectionDate)]
species_event_dt[, SiteNumber := standardize_site(SiteNumber)]
species_event_dt[, CollectionDate := as.Date(CollectionDate)]

event_size_dt <- event_size_dt[!is.na(SiteNumber) & !is.na(CollectionDate) & grepl("^USGS-[0-9]{8}$", SiteNumber)]
species_event_dt <- species_event_dt[!is.na(SiteNumber) & !is.na(CollectionDate) & grepl("^USGS-[0-9]{8}$", SiteNumber)]

events_all <- unique(rbindlist(list(
  event_size_dt[, .(SiteNumber, CollectionDate)],
  species_event_dt[, .(SiteNumber, CollectionDate)]
), use.names = TRUE))
setorder(events_all, SiteNumber, CollectionDate)
events_all[, event_id__ := .I]
if (!nrow(events_all)) stop("No valid events after standardization/filtering in fish event files.", call. = FALSE)

chem_file <- find_one(raw_chem_candidates, "wide chemistry")
chem <- fread(chem_file)

if (!"SiteNumber" %in% names(chem) && "MonitoringLocationIdentifier" %in% names(chem)) {
  setnames(chem, "MonitoringLocationIdentifier", "SiteNumber")
}
if (!"CollectionDate" %in% names(chem) && "ActivityStartDate" %in% names(chem)) {
  setnames(chem, "ActivityStartDate", "CollectionDate")
}
must_have_cols(chem, c("SiteNumber", "CollectionDate", "ActivityIdentifier"), "chemistry data")

chem[, SiteNumber := standardize_site(SiteNumber)]
chem[, CollectionDate := as.Date(CollectionDate)]
chem <- chem[!is.na(SiteNumber) & !is.na(CollectionDate) & grepl("^USGS-[0-9]{8}$", SiteNumber)]
chem <- chem[SiteNumber %in% unique(events_all$SiteNumber)]

chem_cols_all <- grep("^chem_", names(chem), value = TRUE)
if (length(chem_cols_all) < 10) {
  stop(sprintf("Chemistry table has too few chem_* columns (%d found, need >=10).", length(chem_cols_all)), call. = FALSE)
}

chem[, (chem_cols_all) := lapply(.SD, as_num), .SDcols = chem_cols_all]
chem[, nChems := rowSums(!is.na(.SD)), .SDcols = chem_cols_all]
chem_dense <- chem[nChems >= DENSE_MIN_CHEMS]
if (!nrow(chem_dense)) stop("No dense chemistry rows found at DENSE_MIN_CHEMS threshold.", call. = FALSE)

coverage_vec <- chem_dense[, lapply(.SD, function(x) mean(!is.na(x))), .SDcols = chem_cols_all]
coverage_long <- melt(coverage_vec, measure.vars = names(coverage_vec), variable.name = "chem_col", value.name = "coverage")
panel_chems <- coverage_long[coverage >= CHEM_COVERAGE_MIN, chem_col]
panel_chems <- setdiff(panel_chems, EXCLUDE_CHEMS)
if (length(panel_chems) < 5) {
  stop(sprintf("Too few panel chemicals selected (%d) using CHEM_COVERAGE_MIN=%.3f.", length(panel_chems), CHEM_COVERAGE_MIN), call. = FALSE)
}

stdtox_file <- find_one(stdtox_candidates, "StandardTox")
stdtox <- fread(stdtox_file)
stdtox_lookup <- setNames(names(stdtox), safe_name(names(stdtox)))
req_std <- c("group", "endpoint", "duration_unit", "duration", "qualifier", "concentration", "concentration_unit", "cas")
req_std_safe <- safe_name(req_std)
miss_std <- req_std_safe[!req_std_safe %in% names(stdtox_lookup)]
if (length(miss_std)) {
  stop(sprintf("StandardTox missing required columns (safe_name): %s", paste(miss_std, collapse = ", ")), call. = FALSE)
}

stdtox_sub <- stdtox[, ..stdtox_lookup[req_std_safe]]
setnames(stdtox_sub, names(stdtox_sub), req_std)
stdtox_sub[, group := tolower(trimws(as.character(group)))]
stdtox_sub[, endpoint := toupper(trimws(as.character(endpoint)))]
stdtox_sub[, duration_unit := tolower(trimws(as.character(duration_unit)))]
stdtox_sub[, qualifier := trimws(as.character(qualifier))]
stdtox_sub[, duration := as_num(duration)]
stdtox_sub[, concentration := as_num(concentration)]
stdtox_sub[, concentration_unit := tolower(trimws(as.character(concentration_unit)))]
stdtox_sub[, cas := trimws(as.character(cas))]

stdtox_sub <- stdtox_sub[
  group == "fish" & endpoint %in% c("LC50", "LC50*") & duration_unit == "h" & duration == 96 & qualifier == "="
]

stdtox_sub[grepl("ug/l", concentration_unit, fixed = TRUE), concentration_unit := "ug/l"]
stdtox_sub[, concentration_std_ugL := fifelse(
  concentration_unit == "g/l", concentration * 1e6,
  fifelse(concentration_unit == "mg/l", concentration * 1e3,
    fifelse(concentration_unit == "ug/l", concentration,
      fifelse(concentration_unit == "ng/l", concentration / 1000,
        fifelse(concentration_unit == "ppb", concentration,
          fifelse(concentration_unit == "g/m3", concentration * 1e3, NA_real_)
        )
      )
    )
  )
)]

stdtox_sub <- stdtox_sub[is.finite(concentration_std_ugL) & concentration_std_ugL > 0 & !is.na(cas) & cas != ""]
lc50_bench <- stdtox_sub[, .(lc50_bench_ugL = median(concentration_std_ugL, na.rm = TRUE), tox_n_records = .N), by = cas]
if (!nrow(lc50_bench)) stop("No valid fish 96h LC50 benchmark rows after filtering.", call. = FALSE)

cas_file <- find_one(cas_candidates, "CAS map")
cas_map <- fread(cas_file)
must_have_cols(cas_map, c("cas", "CharacteristicName"), "CAS map")
cas_map[, cas := trimws(as.character(cas))]
cas_map[, CharacteristicName := trimws(as.character(CharacteristicName))]

map_dt <- data.table(chem_col = panel_chems, CharacteristicName = sub("^chem_", "", panel_chems))
map_dt <- merge(map_dt, cas_map[, .(CharacteristicName, cas)], by = "CharacteristicName", all.x = TRUE, sort = FALSE)
map_dt <- merge(map_dt, lc50_bench, by = "cas", all.x = TRUE, sort = FALSE)
map_dt <- map_dt[!is.na(cas) & cas != "" & !is.na(lc50_bench_ugL) & is.finite(lc50_bench_ugL) & lc50_bench_ugL > 0]

sel_chems <- unique(map_dt$chem_col)
if (!length(sel_chems)) stop("No selected chemistry columns successfully mapped to CAS+LC50 benchmark.", call. = FALSE)

events_join <- copy(events_all)
events_join[, `:=`(start = CollectionDate - BUFF_DAYS, end = CollectionDate + BUFF_DAYS)]
chem_sub <- chem[, c("SiteNumber", "CollectionDate", "ActivityIdentifier", sel_chems), with = FALSE]
setkey(chem_sub, SiteNumber, CollectionDate)

cand <- chem_sub[events_join, on = .(SiteNumber, CollectionDate >= start, CollectionDate <= end), allow.cartesian = TRUE, nomatch = 0L]

if (nrow(cand)) {
  cand[, `:=`(event_date = i.CollectionDate, chem_date = CollectionDate, event_id__ = i.event_id__)]
  cand[, lag_days := abs(as.integer(chem_date - event_date))]
  setorder(cand, event_id__, lag_days, chem_date, ActivityIdentifier)
  near_one <- cand[, .SD[1], by = event_id__]
  near_evt <- near_one[, c("event_id__", "SiteNumber", "event_date", "ActivityIdentifier", "lag_days", sel_chems), with = FALSE]
  setnames(near_evt, "event_date", "CollectionDate")
  setnames(near_evt, "lag_days", "chem_lag_days")
} else {
  near_evt <- events_all[, .(event_id__, SiteNumber, CollectionDate)]
  near_evt[, `:=`(ActivityIdentifier = NA_character_, chem_lag_days = NA_integer_)]
  for (cc in sel_chems) near_evt[, (cc) := NA_real_]
}

lc50_by_chem <- setNames(map_dt$lc50_bench_ugL, map_dt$chem_col)
for (cc in sel_chems) {
  near_evt[, (paste0("rq__", cc)) := get(cc) / lc50_by_chem[[cc]]]
}
rq_cols <- paste0("rq__", sel_chems)
near_evt[, n_chems_used := rowSums(!is.na(.SD)), .SDcols = rq_cols]
near_evt[, TU_event := rowSums(.SD, na.rm = TRUE), .SDcols = rq_cols]
near_evt[n_chems_used < EVENT_MIN_CHEMS_USED, TU_event := NA_real_]
near_evt[is.na(TU_event) | !is.finite(TU_event), logTU := NA_real_]
near_evt[is.finite(TU_event), logTU := log10(TU_event + EPS_TU)]

tu_event <- near_evt[, .(SiteNumber, CollectionDate, TU_event, logTU, n_chems_used, chem_lag_days)]

model_dt_event <- merge(event_size_dt, tu_event, by = c("SiteNumber", "CollectionDate"), all.x = TRUE, sort = FALSE)
model_dt_species_event <- merge(species_event_dt, tu_event, by = c("SiteNumber", "CollectionDate"), all.x = TRUE, sort = FALSE)

dir.create("output/data", recursive = TRUE, showWarnings = FALSE)
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)

file_panel <- file.path("output/data", sprintf("chem_panel_for_fishSize_%s.csv", today_tag))
file_tu <- file.path("output/data", sprintf("event_TU_for_size_events_highcoverage_%s.csv", today_tag))
file_model_event <- file.path("output/data", sprintf("modeldata_event_TU_SIZE_highcoverage_%s.csv", today_tag))
file_model_species <- file.path("output/data", sprintf("modeldata_species_event_TU_SIZE_highcoverage_%s.csv", today_tag))
file_qc <- file.path("output/results", sprintf("TU_fishSize_highcoverage_QC_%s.txt", today_tag))

fwrite(map_dt[, .(chem_col, CharacteristicName, cas, tox_n_records, lc50_bench_ugL)], file_panel)
fwrite(tu_event, file_tu)
fwrite(model_dt_event, file_model_event)
fwrite(model_dt_species_event, file_model_species)

qc_lines <- c(
  "TU fish size high-coverage QC",
  sprintf("date_tag: %s", today_tag),
  "",
  "Settings:",
  sprintf("BUFF_DAYS: %d", BUFF_DAYS),
  sprintf("DENSE_MIN_CHEMS: %d", DENSE_MIN_CHEMS),
  sprintf("CHEM_COVERAGE_MIN: %.3f", CHEM_COVERAGE_MIN),
  sprintf("EVENT_MIN_CHEMS_USED: %d", EVENT_MIN_CHEMS_USED),
  sprintf("EXCLUDE_CHEMS: %s", paste(EXCLUDE_CHEMS, collapse = " | ")),
  sprintf("EPS_TU: %g", EPS_TU),
  "",
  "Counts:",
  sprintf("n_event_size_rows: %d", nrow(event_size_dt)),
  sprintf("n_species_event_rows: %d", nrow(species_event_dt)),
  sprintf("n_unique_events_all: %d", nrow(events_all)),
  sprintf("n_events_with_any_chem_match: %d", tu_event[!is.na(chem_lag_days), .N]),
  sprintf("n_events_with_TU: %d", tu_event[!is.na(TU_event), .N]),
  sprintf("n_selected_chems_with_LC50: %d", length(sel_chems)),
  "",
  "Summaries:",
  sprintf("summary(n_chems_used): %s", paste(capture.output(summary(tu_event$n_chems_used)), collapse = " ")),
  sprintf("summary(chem_lag_days): %s", paste(capture.output(summary(tu_event$chem_lag_days)), collapse = " ")),
  sprintf("summary(logTU): %s", paste(capture.output(summary(tu_event$logTU)), collapse = " ")),
  "",
  "Source files used:",
  sprintf("event_size_file: %s", event_size_file),
  sprintf("species_event_file: %s", species_event_file),
  sprintf("chem_file: %s", chem_file),
  sprintf("stdtox_file: %s", stdtox_file),
  sprintf("cas_file: %s", cas_file)
)
writeLines(qc_lines, con = file_qc)

cat("Wrote:\n")
cat(file_panel, "\n", sep = "")
cat(file_tu, "\n", sep = "")
cat(file_model_event, "\n", sep = "")
cat(file_model_species, "\n", sep = "")
cat(file_qc, "\n", sep = "")
