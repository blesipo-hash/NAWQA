suppressPackageStartupMessages({
  library(data.table)
})

today_tag <- format(Sys.Date(), "%Y%m%d")

dir.create("output/data", recursive = TRUE, showWarnings = FALSE)
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)

abun_event_out <- sprintf("output/data/fish_abundance_event_covariates_%s.csv", today_tag)
model_q1q2_out <- sprintf("output/data/modeldata_Q1Q2_event_%s.csv", today_tag)
model_q3_out   <- sprintf("output/data/modeldata_Q3_species_event_%s.csv", today_tag)
qc_out         <- sprintf("output/results/abundance_covariate_merge_QC_%s.txt", today_tag)

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
    digs8 <- ifelse(
      nchar(digs) >= 8,
      substr(digs, nchar(digs) - 7, nchar(digs)),
      sprintf("%08s", digs)
    )
    digs8 <- gsub(" ", "0", digs8)
    out[keep] <- digs8
  }

  ifelse(!is.na(out), paste0("USGS-", out), NA_character_)
}

find_latest_file <- function(pattern) {
  files <- Sys.glob(pattern)
  if (!length(files)) return(NA_character_)
  info <- file.info(files)
  files[order(info$mtime, decreasing = TRUE)][1]
}

find_existing_file <- function(filename, search_dirs = c(".", "../data", "data", "output/data")) {
  candidates <- file.path(search_dirs, filename)
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) {
    stop(sprintf("Could not find %s in: %s", filename, paste(search_dirs, collapse = ", ")))
  }
  hits[1]
}

rename_first_match <- function(dt, target, candidates) {
  if (target %in% names(dt)) return(dt)
  hit <- candidates[candidates %in% names(dt)][1]
  if (is.na(hit)) {
    stop(sprintf("Missing required column '%s' and no candidate found (%s).",
                 target, paste(candidates, collapse = ", ")))
  }
  setnames(dt, hit, target)
  dt
}

nearest_date_merge <- function(size_dt, abun_dt, exact_threshold = 0.70, max_lag_days = 90) {
  size_dt <- copy(size_dt)
  abun_dt <- copy(abun_dt)

  size_dt[, CollectionDate := as.Date(CollectionDate)]
  abun_dt[, CollectionDate := as.Date(CollectionDate)]

  abun_keep <- unique(abun_dt[, .(SiteNumber, CollectionDate, total_abun_cov, richness_cov)])

  merged_exact <- merge(
    size_dt,
    abun_keep,
    by = c("SiteNumber", "CollectionDate"),
    all.x = TRUE,
    sort = FALSE
  )

  merged_exact[, abun_match_lag_days := as.integer(NA)]
  exact_rows <- merged_exact[!is.na(total_abun_cov) | !is.na(richness_cov), .N]
  exact_rate <- if (nrow(merged_exact) > 0) exact_rows / nrow(merged_exact) else 0

  if (exact_rate >= exact_threshold || nrow(size_dt) == 0) {
    return(list(
      merged = merged_exact,
      used_nearest = FALSE,
      exact_match_count = exact_rows,
      exact_match_rate = exact_rate,
      nearest_match_count = 0L
    ))
  }

  size_with_id <- copy(size_dt)
  size_with_id[, row_id__ := .I]

  candidate <- merge(
    size_with_id[, .(row_id__, SiteNumber, size_date = CollectionDate)],
    abun_keep[, .(SiteNumber, abun_date = CollectionDate, total_abun_cov, richness_cov)],
    by = "SiteNumber",
    allow.cartesian = TRUE
  )

  if (nrow(candidate) > 0) {
    candidate[, lag_days := abs(as.integer(size_date - abun_date))]
    candidate <- candidate[lag_days <= max_lag_days]
  }

  nearest <- candidate[order(row_id__, lag_days, abun_date)][, .SD[1], by = row_id__]

  merged_nearest <- merge(
    size_with_id,
    nearest[, .(row_id__, total_abun_cov, richness_cov, abun_match_lag_days = lag_days)],
    by = "row_id__",
    all.x = TRUE,
    sort = FALSE
  )
  setorder(merged_nearest, row_id__)
  merged_nearest[, row_id__ := NULL]

  nearest_rows <- merged_nearest[!is.na(total_abun_cov) | !is.na(richness_cov), .N]

  list(
    merged = merged_nearest,
    used_nearest = TRUE,
    exact_match_count = exact_rows,
    exact_match_rate = exact_rate,
    nearest_match_count = nearest_rows - exact_rows
  )
}

fish_abun_file <- find_existing_file("fishData_filtered_sitesWithChemData_011326.csv")
size_event_file <- find_latest_file("output/data/fish_event_size_metrics_*.csv")
size_species_file <- find_latest_file("output/data/fish_species_event_size_metrics_*.csv")

if (is.na(size_event_file) || is.na(size_species_file)) {
  stop("Could not find fish size outputs in output/data/. Run Step 1 scripts first.")
}

fish_abun <- fread(fish_abun_file)

fish_abun <- rename_first_match(
  fish_abun,
  target = "SiteNumber",
  candidates = c("MonitoringLocationIdentifier", "site_no", "siteNumber", "SITE_NO")
)
fish_abun <- rename_first_match(
  fish_abun,
  target = "CollectionDate",
  candidates = c("ActivityStartDate", "Date", "SampleDate", "ActivityDate", "COLLECTION_DATE")
)

fish_abun[, CollectionDate := as.Date(CollectionDate)]
fish_abun[, SiteNumber := standardize_site(SiteNumber)]
fish_abun <- fish_abun[!is.na(SiteNumber) & grepl("^USGS-\\d{8}$", SiteNumber) & !is.na(CollectionDate)]

meta_candidates <- c(
  "SiteNumber", "CollectionDate", "Latitude_dd", "Longitude_dd", "ActivityIdentifier",
  "Agency", "ProjectIdentifier", "Sample", "Gear", "DrainageAreaSqKm", "HUC", "Ecoregion"
)
meta_cols <- intersect(meta_candidates, names(fish_abun))
taxon_cols <- setdiff(names(fish_abun), meta_cols)

if (!length(taxon_cols)) {
  stop("No taxon columns detected after removing metadata columns.")
}

for (col in taxon_cols) {
  vals <- suppressWarnings(as.numeric(fish_abun[[col]]))
  vals[is.infinite(vals)] <- NA_real_
  set(fish_abun, j = col, value = vals)
}

fish_abun[, total_abun_cov := rowSums(.SD, na.rm = TRUE), .SDcols = taxon_cols]
fish_abun[, richness_cov := rowSums(.SD > 0, na.rm = TRUE), .SDcols = taxon_cols]

fish_abun_event <- fish_abun[, .(
  total_abun_cov = sum(total_abun_cov, na.rm = TRUE),
  richness_cov = max(richness_cov, na.rm = TRUE)
), by = .(SiteNumber, CollectionDate)]
fish_abun_event[!is.finite(richness_cov), richness_cov := NA_real_]

setorder(fish_abun_event, SiteNumber, CollectionDate)
fwrite(fish_abun_event, abun_event_out)

size_event <- fread(size_event_file)
size_species <- fread(size_species_file)

for (dt_name in c("size_event", "size_species")) {
  dt <- get(dt_name)
  missing_keys <- setdiff(c("SiteNumber", "CollectionDate"), names(dt))
  if (length(missing_keys)) {
    stop(sprintf("%s missing required columns: %s", dt_name, paste(missing_keys, collapse = ", ")))
  }
  dt[, SiteNumber := standardize_site(SiteNumber)]
  dt[, CollectionDate := as.Date(CollectionDate)]
  dt <- dt[!is.na(SiteNumber) & grepl("^USGS-\\d{8}$", SiteNumber) & !is.na(CollectionDate)]
  assign(dt_name, dt)
}

merge_event <- nearest_date_merge(size_event, fish_abun_event, exact_threshold = 0.70, max_lag_days = 90)
model_q1q2 <- merge_event$merged
setorder(model_q1q2, SiteNumber, CollectionDate)
fwrite(model_q1q2, model_q1q2_out)

merge_species <- nearest_date_merge(size_species, fish_abun_event, exact_threshold = 0.70, max_lag_days = 90)
model_q3 <- merge_species$merged
if ("species_id" %in% names(model_q3)) {
  setorder(model_q3, species_id, SiteNumber, CollectionDate)
} else {
  setorder(model_q3, SiteNumber, CollectionDate)
}
fwrite(model_q3, model_q3_out)

unmatched <- model_q1q2[is.na(total_abun_cov) & is.na(richness_cov), .(SiteNumber, CollectionDate)]
unmatched_show <- head(unique(unmatched), 10)

qc_lines <- c(
  sprintf("size_events=%d", nrow(size_event)),
  sprintf("fish_abun_events=%d", nrow(fish_abun_event)),
  sprintf("exact_match_count=%d", merge_event$exact_match_count),
  sprintf("exact_match_rate=%.4f", merge_event$exact_match_rate),
  sprintf("nearest_date_match_count=%d", merge_event$nearest_match_count),
  sprintf("nearest_used=%s", ifelse(merge_event$used_nearest, "yes", "no")),
  "",
  "abun_match_lag_days_summary:",
  paste(capture.output(summary(model_q1q2$abun_match_lag_days)), collapse = "\n"),
  "",
  "total_abun_cov_summary:",
  paste(capture.output(summary(fish_abun_event$total_abun_cov)), collapse = "\n"),
  "",
  "richness_cov_summary:",
  paste(capture.output(summary(fish_abun_event$richness_cov)), collapse = "\n"),
  "",
  "first_10_unmatched_size_events:",
  if (nrow(unmatched_show) == 0) "none" else paste(capture.output(print(unmatched_show)), collapse = "\n"),
  "",
  sprintf("source_fish_abundance_file=%s", fish_abun_file),
  sprintf("source_size_event_file=%s", size_event_file),
  sprintf("source_size_species_file=%s", size_species_file)
)

writeLines(qc_lines, qc_out)

cat(sprintf("Wrote %s\n", abun_event_out))
cat(sprintf("Wrote %s\n", model_q1q2_out))
cat(sprintf("Wrote %s\n", model_q3_out))
cat(sprintf("Wrote %s\n", qc_out))
