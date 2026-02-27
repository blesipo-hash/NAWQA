suppressPackageStartupMessages({
  library(data.table)
  library(finsyncR)
})

today_tag <- format(Sys.Date(), "%Y%m%d")

dir.create("output/data", recursive = TRUE, showWarnings = FALSE)
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)

finsync_event_out <- sprintf("output/data/finsync_fish_abundance_event_MGMS_%s.csv", today_tag)
model_q1q2_out    <- sprintf("output/data/modeldata_Q1Q2_event_%s.csv", today_tag)
model_q3_out      <- sprintf("output/data/modeldata_Q3_species_event_%s.csv", today_tag)
qc_out            <- sprintf("output/results/finsync_abundance_merge_QC_%s.txt", today_tag)

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
  if (!length(files)) {
    stop(sprintf("No files found matching pattern: %s", pattern))
  }
  info <- file.info(files)
  files[order(info$mtime, decreasing = TRUE)][1]
}

nearest_date_merge <- function(size_dt, abun_dt, exact_threshold = 0.70, max_lag_days = 90) {
  size_dt <- copy(size_dt)
  abun_dt <- copy(abun_dt)

  size_dt[, CollectionDate := as.Date(CollectionDate)]
  abun_dt[, CollectionDate := as.Date(CollectionDate)]

  abun_keep <- unique(abun_dt[, .(SiteNumber, CollectionDate, total_abun_MGMS, richness_MGMS)])

  merged_exact <- merge(
    size_dt,
    abun_keep,
    by = c("SiteNumber", "CollectionDate"),
    all.x = TRUE,
    sort = FALSE
  )

  matched_exact <- merged_exact[!is.na(total_abun_MGMS) | !is.na(richness_MGMS), .N]
  match_rate <- if (nrow(merged_exact) > 0) matched_exact / nrow(merged_exact) else 0

  merged_exact[, abun_match_lag_days := as.integer(NA)]

  if (match_rate >= exact_threshold || nrow(size_dt) == 0) {
    return(list(
      merged = merged_exact,
      used_nearest = FALSE,
      exact_match_count = matched_exact,
      nearest_match_count = 0L,
      match_rate = match_rate
    ))
  }

  size_with_id <- copy(size_dt)
  size_with_id[, row_id__ := .I]

  candidate <- merge(
    size_with_id[, .(row_id__, SiteNumber, size_date = CollectionDate)],
    abun_keep[, .(SiteNumber, abun_date = CollectionDate, total_abun_MGMS, richness_MGMS)],
    by = "SiteNumber",
    allow.cartesian = TRUE
  )

  if (nrow(candidate) > 0) {
    candidate[, lag_days := as.integer(abs(as.numeric(size_date - abun_date)))]
    candidate <- candidate[lag_days <= max_lag_days]
  }

  nearest <- candidate[order(row_id__, lag_days, abun_date)][, .SD[1], by = row_id__]

  out <- merge(size_with_id, nearest[, .(row_id__, total_abun_MGMS, richness_MGMS, abun_match_lag_days = lag_days)],
               by = "row_id__", all.x = TRUE, sort = FALSE)
  setorder(out, row_id__)
  out[, row_id__ := NULL]

  exact_key <- merged_exact[!is.na(total_abun_MGMS) | !is.na(richness_MGMS), .(SiteNumber, CollectionDate)]
  exact_unique <- unique(exact_key)
  nearest_key <- out[!is.na(abun_match_lag_days), .(SiteNumber, CollectionDate)]
  nearest_unique <- unique(nearest_key)

  nearest_only <- fsetdiff(nearest_unique, exact_unique)

  list(
    merged = out,
    used_nearest = TRUE,
    exact_match_count = nrow(exact_unique),
    nearest_match_count = nrow(nearest_only),
    match_rate = match_rate
  )
}

# ---- pull finsync abundance ----
fish_abun <- as.data.table(getFishData(dataType = "abun", agency = "USGS", standardize = "MGMS"))

required_abun_cols <- c("SiteNumber", "CollectionDate")
missing_abun_cols <- setdiff(required_abun_cols, names(fish_abun))
if (length(missing_abun_cols)) {
  stop(sprintf("finsyncR fish abundance missing required columns: %s", paste(missing_abun_cols, collapse = ", ")))
}

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
  stop("No taxon columns detected in finsync abundance table after removing metadata columns.")
}

for (col in taxon_cols) {
  set(fish_abun, j = col, value = suppressWarnings(as.numeric(fish_abun[[col]])))
}

fish_abun_event <- fish_abun[, {
  vals <- .SD
  list(
    total_abun_MGMS = rowSums(vals, na.rm = TRUE),
    richness_MGMS = rowSums(vals > 0, na.rm = TRUE)
  )
}, by = .(SiteNumber, CollectionDate), .SDcols = taxon_cols]

setorder(fish_abun_event, SiteNumber, CollectionDate)
fwrite(fish_abun_event, finsync_event_out)

# ---- load most recent fish-size outputs ----
size_event_file <- find_latest_file("output/data/fish_event_size_metrics_*.csv")
size_species_file <- find_latest_file("output/data/fish_species_event_size_metrics_*.csv")

size_event <- fread(size_event_file)
size_species <- fread(size_species_file)

for (dt_name in c("size_event", "size_species")) {
  dt <- get(dt_name)
  miss <- setdiff(c("SiteNumber", "CollectionDate"), names(dt))
  if (length(miss)) {
    stop(sprintf("%s missing required columns: %s", dt_name, paste(miss, collapse = ", ")))
  }
  dt[, CollectionDate := as.Date(CollectionDate)]
  dt[, SiteNumber := standardize_site(SiteNumber)]
  assign(dt_name, dt)
}

# ---- merge Q1/Q2 event model data ----
merge_event <- nearest_date_merge(size_event, fish_abun_event, exact_threshold = 0.70, max_lag_days = 90)
model_q1q2 <- merge_event$merged
setorder(model_q1q2, SiteNumber, CollectionDate)
fwrite(model_q1q2, model_q1q2_out)

# ---- merge Q3 species-event model data ----
merge_q3 <- nearest_date_merge(size_species, fish_abun_event, exact_threshold = 0.70, max_lag_days = 90)
model_q3 <- merge_q3$merged
setorder(model_q3, species_id, SiteNumber, CollectionDate)
fwrite(model_q3, model_q3_out)

# ---- QC ----
unmatched <- model_q1q2[is.na(total_abun_MGMS) & is.na(richness_MGMS), .(SiteNumber, CollectionDate)]
unmatched_show <- head(unique(unmatched), 10)

lag_summary <- summary(model_q1q2$abun_match_lag_days)
abun_summary <- summary(model_q1q2$total_abun_MGMS)
rich_summary <- summary(model_q1q2$richness_MGMS)

qc_lines <- c(
  sprintf("size_events=%d", nrow(size_event)),
  sprintf("finsync_events=%d", nrow(fish_abun_event)),
  sprintf("exact_match_count=%d", merge_event$exact_match_count),
  sprintf("nearest_date_match_count=%d", merge_event$nearest_match_count),
  sprintf("exact_match_rate=%.4f", merge_event$match_rate),
  sprintf("nearest_used=%s", ifelse(merge_event$used_nearest, "yes", "no")),
  "",
  "abun_match_lag_days_summary:",
  paste(capture.output(print(lag_summary)), collapse = "\n"),
  "",
  "total_abun_MGMS_summary:",
  paste(capture.output(print(abun_summary)), collapse = "\n"),
  "",
  "richness_MGMS_summary:",
  paste(capture.output(print(rich_summary)), collapse = "\n"),
  "",
  "first_10_unmatched_events:",
  if (nrow(unmatched_show) == 0) "none" else paste(capture.output(print(unmatched_show)), collapse = "\n"),
  "",
  sprintf("source_size_event_file=%s", size_event_file),
  sprintf("source_size_species_file=%s", size_species_file)
)

writeLines(qc_lines, qc_out)

cat(sprintf("Wrote %s\n", finsync_event_out))
cat(sprintf("Wrote %s\n", model_q1q2_out))
cat(sprintf("Wrote %s\n", model_q3_out))
cat(sprintf("Wrote %s\n", qc_out))
