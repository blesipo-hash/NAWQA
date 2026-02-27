suppressPackageStartupMessages({
  library(data.table)
})

today_tag <- format(Sys.Date(), "%Y%m%d")
BUFF_DAYS <- 90L
MIN_CHEMS_PER_EVENT <- 210L
EXCLUDE_CHEMS <- c("chem_Azinphos-methyl oxygen analog", "chem_Terbufos oxygen analog sulfone")

dir.create("output/data", recursive = TRUE, showWarnings = FALSE)
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)

size_file <- "output/data/fish_event_size_metrics_20260226.csv"
raw_chem_candidates <- c("data/chemData_wide_071025.csv", "../data/chemData_wide_071025.csv", "chemData_wide_071025.csv")
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

as_num <- function(x) {
  y <- suppressWarnings(as.numeric(x))
  y[is.infinite(y)] <- NA_real_
  y
}

safe_name <- function(x) tolower(gsub("[^a-z0-9]+", "", x))

if (!file.exists(size_file)) stop(sprintf("Missing required fish size events file: %s", size_file), call. = FALSE)
raw_chem_file <- find_one(raw_chem_candidates, "raw wide chemistry table")
stdtox_file <- find_one(stdtox_candidates, "StandardTox CSV")
cas_file <- find_one(cas_candidates, "CAS mapping")

size <- fread(size_file)
chem <- fread(raw_chem_file)
stdtox <- fread(stdtox_file)
cas_map <- fread(cas_file)

if (!all(c("SiteNumber", "CollectionDate") %in% names(size))) {
  stop("Fish size file must contain SiteNumber and CollectionDate.", call. = FALSE)
}

# normalize raw chem key columns from common names
if (!("SiteNumber" %in% names(chem))) {
  if ("MonitoringLocationIdentifier" %in% names(chem)) {
    setnames(chem, "MonitoringLocationIdentifier", "SiteNumber")
  } else stop("Raw chem file missing SiteNumber (or MonitoringLocationIdentifier).", call. = FALSE)
}
if (!("CollectionDate" %in% names(chem))) {
  if ("ActivityStartDate" %in% names(chem)) {
    setnames(chem, "ActivityStartDate", "CollectionDate")
  } else stop("Raw chem file missing CollectionDate (or ActivityStartDate).", call. = FALSE)
}
if (!("ActivityIdentifier" %in% names(chem))) stop("Raw chem file must include ActivityIdentifier.", call. = FALSE)

if (!all(c("cas", "CharacteristicName") %in% names(cas_map))) {
  stop("CAS mapping file must contain columns: cas, CharacteristicName", call. = FALSE)
}

size[, SiteNumber := standardize_site(SiteNumber)]
size[, CollectionDate := as.Date(CollectionDate)]
size <- size[!is.na(SiteNumber) & grepl("^USGS-\\d{8}$", SiteNumber) & !is.na(CollectionDate)]

chem[, SiteNumber := standardize_site(SiteNumber)]
chem[, CollectionDate := as.Date(CollectionDate)]
chem <- chem[!is.na(SiteNumber) & grepl("^USGS-\\d{8}$", SiteNumber) & !is.na(CollectionDate)]

chem_cols_all <- grep("^chem_", names(chem), value = TRUE)
if (length(chem_cols_all) < 10) stop(sprintf("Expected broad raw chemistry; found only %d chem_* columns.", length(chem_cols_all)), call. = FALSE)
for (cc in chem_cols_all) set(chem, j = cc, value = as_num(chem[[cc]]))

# ---- choose repeated broad panel (template logic) ----
chem[, nChems := rowSums(!is.na(.SD)), .SDcols = chem_cols_all]
chem_dense <- chem[nChems >= MIN_CHEMS_PER_EVENT]
if (!nrow(chem_dense)) stop(sprintf("No chemistry rows meet MIN_CHEMS_PER_EVENT=%d", MIN_CHEMS_PER_EVENT), call. = FALSE)

presence <- colSums(!is.na(chem_dense[, ..chem_cols_all]))
shared_chems <- names(presence)[presence == nrow(chem_dense)]
shared_chems <- setdiff(shared_chems, EXCLUDE_CHEMS)
if (length(shared_chems) < 5) stop("Too few shared chemicals after dense-event intersection and exclusions.", call. = FALSE)

# retain activities where all shared chems are present
activity_ok <- chem[, .(all_shared_present = all(!is.na(unlist(.SD)))), by = ActivityIdentifier, .SDcols = shared_chems]
keep_activities <- activity_ok[all_shared_present == TRUE, ActivityIdentifier]
chem_reduced <- chem[ActivityIdentifier %in% keep_activities, c("ActivityIdentifier", "CollectionDate", "SiteNumber", shared_chems), with = FALSE]

# restrict to sites with fish size events
common_sites <- intersect(unique(chem_reduced$SiteNumber), unique(size$SiteNumber))
chem_reduced <- chem_reduced[SiteNumber %in% common_sites]

wide_panel_out <- sprintf("output/data/chemData_wide_shared_panel_for_fishSize_%s.csv", today_tag)
fwrite(chem_reduced, wide_panel_out)

# ---- align chemistry to fish-size events (±90d) ----
size_ev <- copy(size)
size_ev[, event_id__ := .I]
size_ev[, start := CollectionDate - BUFF_DAYS]
size_ev[, end := CollectionDate + BUFF_DAYS]

setkey(chem_reduced, SiteNumber, CollectionDate)
joined <- chem_reduced[size_ev, on = .(SiteNumber, CollectionDate >= start, CollectionDate <= end), allow.cartesian = TRUE, nomatch = 0L]

if (!nrow(joined)) stop("No chemistry rows matched fish-size events within BUFF_DAYS window.", call. = FALSE)
joined[, lag_days := abs(as.integer(CollectionDate - i.CollectionDate))]

chem_means <- joined[, c(
  list(
    SiteNumber = first(i.SiteNumber),
    CollectionDate = first(i.CollectionDate),
    n_chem_samples_in_window = .N,
    chem_lag_mean_days = mean(lag_days, na.rm = TRUE),
    ActivityIdentifiers = paste(unique(ActivityIdentifier), collapse = "|")
  ),
  lapply(.SD, function(v) mean(v, na.rm = TRUE))
), by = event_id__, .SDcols = shared_chems]
for (cc in shared_chems) chem_means[is.nan(get(cc)), (cc) := NA_real_]

chem_means <- size_ev[chem_means, on = "event_id__"]
chem_means <- chem_means[, c("SiteNumber", "CollectionDate", "n_chem_samples_in_window", "chem_lag_mean_days", "ActivityIdentifiers", shared_chems), with = FALSE]

chem_means_out <- sprintf("output/data/chem_means_for_size_events_widepanel_%s.csv", today_tag)
fwrite(chem_means, chem_means_out)

# ---- StandardTox pooled fish 96h LC50 ----
nm_std <- names(stdtox)
std_idx <- setNames(nm_std, safe_name(nm_std))
need_std <- c("group", "endpoint", "duration_unit", "duration", "qualifier", "concentration", "concentration_unit", "cas")
missing_std <- need_std[!need_std %in% names(std_idx)]
if (length(missing_std)) stop(sprintf("StandardTox missing required columns: %s", paste(missing_std, collapse = ", ")), call. = FALSE)

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

stdtox2 <- stdtox2[group == "fish" & endpoint %in% c("LC50", "LC50*") & duration_unit == "h" & duration == 96 & qualifier == "="]
stdtox2[, concentration_std_ugL := NA_real_]
stdtox2[concentration_unit %in% c("g/l", "g per l"), concentration_std_ugL := concentration * 1e6]
stdtox2[concentration_unit %in% c("mg/l", "mg per l"), concentration_std_ugL := concentration * 1e3]
stdtox2[concentration_unit %in% c("ug/l", "µg/l", "μg/l", "ug per l"), concentration_std_ugL := concentration]
stdtox2[concentration_unit %in% c("ng/l", "ng per l"), concentration_std_ugL := concentration / 1000]
stdtox2[concentration_unit %in% c("ppb"), concentration_std_ugL := concentration]
stdtox2[concentration_unit %in% c("g/m3", "g/m^3"), concentration_std_ugL := concentration * 1e3]
stdtox2[concentration_unit %in% c("l/l", "ul/l", "µl/l", "μl/l"), concentration_std_ugL := NA_real_]
stdtox2 <- stdtox2[is.finite(concentration_std_ugL) & concentration_std_ugL > 0 & !is.na(cas) & nzchar(cas)]

lc50_bench <- stdtox2[, .(lc50_bench_ugL = median(concentration_std_ugL, na.rm = TRUE), tox_n_records = .N), by = cas]

map_dt <- data.table(chem_col = shared_chems)
map_dt[, CharacteristicName := gsub("^chem_", "", chem_col)]
cas_map2 <- unique(cas_map[, .(cas = trimws(as.character(cas)), CharacteristicName = trimws(as.character(CharacteristicName)))])
map_dt <- merge(map_dt, cas_map2, by = "CharacteristicName", all.x = TRUE)
map_dt <- merge(map_dt, lc50_bench, by = "cas", all.x = TRUE)
map_dt <- map_dt[!is.na(cas) & nzchar(cas) & is.finite(lc50_bench_ugL) & lc50_bench_ugL > 0]
if (!nrow(map_dt)) stop("No shared chemicals mapped to CAS with valid pooled fish LC50 benchmarks.", call. = FALSE)

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

tu_out <- sprintf("output/data/event_TU_for_size_events_widepanel_%s.csv", today_tag)
fwrite(tu_event, tu_out)

model_dt <- merge(size, tu_event, by = c("SiteNumber", "CollectionDate"), all.x = TRUE, sort = FALSE)
setorder(model_dt, SiteNumber, CollectionDate)
model_out <- sprintf("output/data/modeldata_Q1Q2_event_TU_SIZE_widepanel_%s.csv", today_tag)
fwrite(model_dt, model_out)

qc_out <- sprintf("output/results/TU_size_event_widepanel_QC_%s.txt", today_tag)
miss_site <- model_dt[!is.finite(logTU), .N, by = SiteNumber][order(-N)][1:min(10, .N)]
qc_lines <- c(
  sprintf("MIN_CHEMS_PER_EVENT=%d", MIN_CHEMS_PER_EVENT),
  sprintf("BUFF_DAYS=%d", BUFF_DAYS),
  sprintf("raw_chem_cols=%d", length(chem_cols_all)),
  sprintf("shared_chems_selected=%d", length(shared_chems)),
  sprintf("mapped_shared_chems=%d", nrow(map_dt)),
  sprintf("n_size_events=%d", nrow(size)),
  sprintf("n_size_events_with_TU=%d", tu_event[is.finite(logTU), .N]),
  "",
  "shared_chems:",
  paste(shared_chems, collapse = ", "),
  "",
  "mapped_chems_CAS_bench:",
  paste(capture.output(print(map_dt[order(chem_col), .(chem_col, cas, tox_n_records, lc50_bench_ugL)])), collapse = "\n"),
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
  "top_10_sites_with_missing_TU:",
  if (nrow(miss_site) == 0) "none" else paste(capture.output(print(miss_site)), collapse = "\n"),
  "",
  sprintf("source_size=%s", size_file),
  sprintf("source_raw_chem=%s", raw_chem_file),
  sprintf("source_standardtox=%s", stdtox_file),
  sprintf("source_cas=%s", cas_file)
)
writeLines(qc_lines, qc_out)

cat(sprintf("Wrote %s\n", wide_panel_out))
cat(sprintf("Wrote %s\n", chem_means_out))
cat(sprintf("Wrote %s\n", tu_out))
cat(sprintf("Wrote %s\n", model_out))
cat(sprintf("Wrote %s\n", qc_out))
