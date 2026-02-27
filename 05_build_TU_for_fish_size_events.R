suppressPackageStartupMessages({
  library(data.table)
})

today_tag <- format(Sys.Date(), "%Y%m%d")
buff_days <- 90L

dir.create("output/data", recursive = TRUE, showWarnings = FALSE)
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)

size_file <- "output/data/fish_event_size_metrics_20260226.csv"
chem_candidates <- c("data/chemData_filtered_fishToxBasis_012226.csv", "../data/chemData_filtered_fishToxBasis_012226.csv", "chemData_filtered_fishToxBasis_012226.csv")
stdtox_candidates <- c("data/standartox_allData_090825.csv", "../data/standartox_allData_090825.csv")
cas_candidates <- c("data/cas_chemname_key.csv", "../data/cas_chemname_key.csv", "cas_chemname_key.csv")

find_one <- function(paths, label) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) stop(sprintf("Missing required input for %s. Checked: %s", label, paste(paths, collapse = ", ")), call. = FALSE)
  hit[1]
}

standardize_site <- function(x) {
  x_chr <- trimws(as.character(x))
  x_chr[x_chr %in% c("", "NA", "NaN")] <- NA_character_
  out <- x_chr

  usgs_now <- !is.na(out) & grepl("^USGS-", out, ignore.case = TRUE)
  if (any(usgs_now)) {
    out[usgs_now] <- gsub("\\D", "", out[usgs_now])
    out[usgs_now][nchar(out[usgs_now]) == 0] <- NA_character_
  }

  numericish <- !is.na(out) & grepl("^[0-9]+$", out)
  other <- !is.na(out) & !usgs_now & !numericish
  if (any(other)) {
    out[other] <- gsub("\\D", "", out[other])
    out[other][nchar(out[other]) == 0] <- NA_character_
  }

  keep <- !is.na(out)
  if (any(keep)) {
    digs <- out[keep]
    digs8 <- ifelse(nchar(digs) >= 8, substr(digs, nchar(digs) - 7, nchar(digs)), sprintf("%08s", digs))
    digs8 <- gsub(" ", "0", digs8)
    out[keep] <- digs8
  }

  ifelse(!is.na(out), paste0("USGS-", out), NA_character_)
}

parse_date <- function(x) as.Date(x)

safe_name <- function(x) tolower(gsub("[^a-z0-9]+", "", x))

as_num <- function(x) {
  y <- suppressWarnings(as.numeric(x))
  y[is.infinite(y)] <- NA_real_
  y
}

size_out <- sprintf("output/data/chem_means_for_size_events_%s.csv", today_tag)
tu_out <- sprintf("output/data/event_TU_for_size_events_%s.csv", today_tag)
model_out <- sprintf("output/data/modeldata_Q1Q2_event_TU_SIZE_%s.csv", today_tag)
qc_out <- sprintf("output/results/TU_size_event_QC_%s.txt", today_tag)

if (!file.exists(size_file)) stop(sprintf("Missing required fish size events file: %s", size_file), call. = FALSE)
chem_file <- find_one(chem_candidates, "chem tox-basis wide table")
stdtox_file <- find_one(stdtox_candidates, "StandardTox CSV")
cas_file <- find_one(cas_candidates, "CAS mapping")

# ---- read inputs ----
size <- fread(size_file)
chem <- fread(chem_file)
stdtox <- fread(stdtox_file)
cas_map <- fread(cas_file)

if (!all(c("SiteNumber", "CollectionDate") %in% names(size))) {
  stop("Fish size file must contain SiteNumber and CollectionDate.", call. = FALSE)
}

chem_required <- c("SiteNumber", "CollectionDate", "ActivityIdentifier")
missing_chem <- setdiff(chem_required, names(chem))
if (length(missing_chem)) {
  stop(sprintf("Chem file missing required columns: %s", paste(missing_chem, collapse = ", ")), call. = FALSE)
}

if (!all(c("cas", "CharacteristicName") %in% names(cas_map))) {
  stop("CAS mapping file must contain columns: cas, CharacteristicName", call. = FALSE)
}

size[, SiteNumber := standardize_site(SiteNumber)]
size[, CollectionDate := parse_date(CollectionDate)]
size <- size[!is.na(SiteNumber) & grepl("^USGS-\\d{8}$", SiteNumber) & !is.na(CollectionDate)]

chem[, SiteNumber := standardize_site(SiteNumber)]
chem[, CollectionDate := parse_date(CollectionDate)]
chem <- chem[!is.na(SiteNumber) & grepl("^USGS-\\d{8}$", SiteNumber) & !is.na(CollectionDate)]

chem_cols <- grep("^chem_", names(chem), value = TRUE)
if (length(chem_cols) < 5) {
  stop(sprintf("Expected at least 5 chem_* columns; found %d.", length(chem_cols)), call. = FALSE)
}
for (cc in chem_cols) set(chem, j = cc, value = as_num(chem[[cc]]))

# ---- non-equi window alignment ----
size_ev <- copy(size)
size_ev[, event_id__ := .I]
size_ev[, start := CollectionDate - buff_days]
size_ev[, end := CollectionDate + buff_days]

setkey(chem, SiteNumber, CollectionDate)
joined <- chem[size_ev, on = .(SiteNumber, CollectionDate >= start, CollectionDate <= end), allow.cartesian = TRUE, nomatch = 0L]

if (nrow(joined) == 0) {
  warning("No chemistry rows matched any fish size event within ±90 days.")
  chem_means <- size_ev[, .(
    SiteNumber,
    CollectionDate,
    n_chem_samples_in_window = 0L,
    chem_lag_mean_days = NA_real_,
    ActivityIdentifiers = NA_character_
  )]
  for (cc in chem_cols) chem_means[, (cc) := NA_real_]
} else {
  joined[, lag_days := abs(as.integer(CollectionDate - i.CollectionDate))]
  joined[, chem_row_date := CollectionDate]

  chem_means <- joined[, c(
    list(
      SiteNumber = first(i.SiteNumber),
      CollectionDate = first(i.CollectionDate),
      n_chem_samples_in_window = uniqueN(ActivityIdentifier),
      chem_lag_mean_days = mean(lag_days, na.rm = TRUE),
      ActivityIdentifiers = paste(unique(ActivityIdentifier), collapse = "|")
    ),
    lapply(.SD, function(v) mean(v, na.rm = TRUE))
  ), by = event_id__, .SDcols = chem_cols]

  for (cc in chem_cols) {
    chem_means[is.nan(get(cc)), (cc) := NA_real_]
  }

  chem_means <- size_ev[chem_means, on = "event_id__"]
  keep <- c("SiteNumber", "CollectionDate", "n_chem_samples_in_window", "chem_lag_mean_days", "ActivityIdentifiers", chem_cols)
  chem_means <- chem_means[, ..keep]
}

setorder(chem_means, SiteNumber, CollectionDate)
fwrite(chem_means, size_out)

# ---- StandardTox filtering + unit standardization ----
nm_std <- names(stdtox)
std_idx <- setNames(nm_std, safe_name(nm_std))
need_std <- c("group", "endpoint", "duration_unit", "duration", "qualifier", "concentration", "concentration_unit", "cas")
missing_std <- need_std[!need_std %in% names(std_idx)]
if (length(missing_std)) {
  stop(sprintf("StandardTox missing required columns (normalized names): %s", paste(missing_std, collapse = ", ")), call. = FALSE)
}

stdtox2 <- copy(stdtox)
setnames(stdtox2, std_idx["group"], "group")
setnames(stdtox2, std_idx["endpoint"], "endpoint")
setnames(stdtox2, std_idx["duration_unit"], "duration_unit")
setnames(stdtox2, std_idx["duration"], "duration")
setnames(stdtox2, std_idx["qualifier"], "qualifier")
setnames(stdtox2, std_idx["concentration"], "concentration")
setnames(stdtox2, std_idx["concentration_unit"], "concentration_unit")
setnames(stdtox2, std_idx["cas"], "cas")

stdtox2[, group := tolower(trimws(as.character(group)))]
stdtox2[, endpoint := toupper(trimws(as.character(endpoint)))]
stdtox2[, duration_unit := tolower(trimws(as.character(duration_unit)))]
stdtox2[, duration := as_num(duration)]
stdtox2[, qualifier := trimws(as.character(qualifier))]
stdtox2[, concentration := as_num(concentration)]
stdtox2[, concentration_unit := tolower(trimws(as.character(concentration_unit)))]
stdtox2[, cas := trimws(as.character(cas))]

stdtox2 <- stdtox2[
  group == "fish" & endpoint %in% c("LC50", "LC50*") &
    duration_unit == "h" & !is.na(duration) & duration == 96 &
    qualifier == "="
]

stdtox2[, concentration_std_ugL := NA_real_]
stdtox2[concentration_unit %in% c("g/l", "g per l"), concentration_std_ugL := concentration * 1e6]
stdtox2[concentration_unit %in% c("mg/l", "mg per l"), concentration_std_ugL := concentration * 1e3]
stdtox2[concentration_unit %in% c("ug/l", "µg/l", "μg/l", "ug per l"), concentration_std_ugL := concentration]
stdtox2[concentration_unit %in% c("ng/l", "ng per l"), concentration_std_ugL := concentration / 1000]
stdtox2[concentration_unit %in% c("ppb"), concentration_std_ugL := concentration]
stdtox2[concentration_unit %in% c("g/m3", "g/m^3"), concentration_std_ugL := concentration * 1e3]
# Conservative handling for l/l or ul/l: drop unless clearly mass-equivalent (not available here)
stdtox2[concentration_unit %in% c("l/l", "ul/l", "µl/l", "μl/l"), concentration_std_ugL := NA_real_]

stdtox2 <- stdtox2[is.finite(concentration_std_ugL) & concentration_std_ugL > 0 & !is.na(cas) & nzchar(cas)]

lc50_bench <- stdtox2[, .(
  lc50_bench_ugL = median(concentration_std_ugL, na.rm = TRUE),
  tox_n_records = .N
), by = cas]

# ---- map chem_* to CAS ----
map_dt <- data.table(chem_col = chem_cols)
map_dt[, CharacteristicName := gsub("^chem_", "", chem_col)]
cas_map2 <- unique(cas_map[, .(cas = trimws(as.character(cas)), CharacteristicName = trimws(as.character(CharacteristicName)))])
map_dt <- merge(map_dt, cas_map2, by = "CharacteristicName", all.x = TRUE)
map_dt <- merge(map_dt, lc50_bench, by = "cas", all.x = TRUE)
map_dt <- map_dt[!is.na(cas) & nzchar(cas) & is.finite(lc50_bench_ugL) & lc50_bench_ugL > 0]

if (!nrow(map_dt)) {
  stop("No chem_* columns successfully mapped to CAS with valid LC50 benchmarks.", call. = FALSE)
}

# ---- compute event TU ----
chem_for_tu <- copy(chem_means)
for (i in seq_len(nrow(map_dt))) {
  cc <- map_dt$chem_col[i]
  rq_col <- paste0("rq__", cc)
  chem_for_tu[, (rq_col) := as_num(get(cc)) / map_dt$lc50_bench_ugL[i]]
}
rq_cols <- grep("^rq__", names(chem_for_tu), value = TRUE)

chem_for_tu[, n_chems_used := rowSums(!is.na(.SD)), .SDcols = rq_cols]
chem_for_tu[, TU_event := rowSums(.SD, na.rm = TRUE), .SDcols = rq_cols]
chem_for_tu[n_chems_used == 0, TU_event := NA_real_]
chem_for_tu[, logTU := ifelse(is.finite(TU_event), log10(TU_event + 1e-8), NA_real_)]

tu_event <- chem_for_tu[, .(SiteNumber, CollectionDate, TU_event, logTU, n_chems_used, n_chem_samples_in_window)]
setorder(tu_event, SiteNumber, CollectionDate)
fwrite(tu_event, tu_out)

# ---- merge to size events ----
model_dt <- merge(size, tu_event, by = c("SiteNumber", "CollectionDate"), all.x = TRUE, sort = FALSE)
setorder(model_dt, SiteNumber, CollectionDate)
fwrite(model_dt, model_out)

# ---- QC report ----
qc_n_size <- nrow(size)
qc_any_chem <- tu_event[!is.na(n_chem_samples_in_window) & n_chem_samples_in_window > 0, .N]
qc_with_tu <- tu_event[is.finite(logTU), .N]

miss_site <- model_dt[!is.finite(logTU), .N, by = SiteNumber][order(-N)][1:min(10, .N)]
map_show <- map_dt[order(chem_col), .(chem_col, cas, tox_n_records, lc50_bench_ugL)]

qc_lines <- c(
  sprintf("n_size_events=%d", qc_n_size),
  sprintf("n_events_with_any_chem_in_window=%d", qc_any_chem),
  sprintf("n_events_with_TU=%d", qc_with_tu),
  "",
  "summary_n_chem_samples_in_window:",
  paste(capture.output(summary(tu_event$n_chem_samples_in_window)), collapse = "\n"),
  "",
  "summary_n_chems_used:",
  paste(capture.output(summary(tu_event$n_chems_used)), collapse = "\n"),
  "",
  "summary_TU_event:",
  paste(capture.output(summary(tu_event$TU_event)), collapse = "\n"),
  "",
  "summary_logTU:",
  paste(capture.output(summary(tu_event$logTU)), collapse = "\n"),
  "",
  "top_10_sites_with_most_missing_TU:",
  if (nrow(miss_site) == 0) "none" else paste(capture.output(print(miss_site)), collapse = "\n"),
  "",
  "mapped_chems_CAS_tox_records_medianLC50:",
  paste(capture.output(print(map_show)), collapse = "\n"),
  "",
  sprintf("source_size_file=%s", size_file),
  sprintf("source_chem_file=%s", chem_file),
  sprintf("source_standardtox_file=%s", stdtox_file),
  sprintf("source_cas_file=%s", cas_file)
)
writeLines(qc_lines, qc_out)

cat(sprintf("Wrote %s\n", size_out))
cat(sprintf("Wrote %s\n", tu_out))
cat(sprintf("Wrote %s\n", model_out))
cat(sprintf("Wrote %s\n", qc_out))
