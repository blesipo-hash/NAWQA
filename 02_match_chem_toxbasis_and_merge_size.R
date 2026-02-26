suppressWarnings({
  library(data.table)
})

today_tag <- format(Sys.Date(), "%Y%m%d")
buff_days <- 90L

fish_file <- "fishData_filtered_fishToxBasis_012226.csv"
chem_file <- "chemData_filtered_fishToxBasis_012226.csv"
cas_key_file <- "cas_chemname_key.csv"
standartox_file <- "standartox_allData_090825.csv"
existing_tu_file <- "fish_toxicUnits_012326.csv"

out_chem_means <- sprintf("chem_event_means_fishToxBasis_%s.csv", today_tag)
out_merged <- sprintf("fish_size_plus_chem_fishToxBasis_%s.csv", today_tag)
out_tu <- sprintf("fish_size_plus_TUevent_fishToxBasis_%s.csv", today_tag)
out_qc <- sprintf("fish_TU_QC_fishToxBasis_%s.txt", today_tag)

standardize_site <- function(x) {
  x_chr <- trimws(as.character(x))
  already_usgs <- grepl("^USGS-", x_chr, ignore.case = TRUE)
  out <- x_chr
  numericish <- grepl("^[0-9]+$", x_chr)
  need_convert <- !already_usgs & numericish
  if (any(need_convert)) {
    out[need_convert] <- sprintf("USGS-%08d", as.integer(x_chr[need_convert]))
  }
  usgs_now <- grepl("^USGS-", out, ignore.case = TRUE)
  if (any(usgs_now)) {
    suffix <- sub("^USGS-", "", toupper(out[usgs_now]))
    digits <- gsub("[^0-9]", "", suffix)
    use_digits <- nchar(digits) > 0
    suffix[use_digits] <- sprintf("%08d", as.integer(digits[use_digits]))
    out[usgs_now] <- paste0("USGS-", suffix)
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

if (!file.exists(fish_file) || !file.exists(chem_file)) {
  stop("Required fish/chem files for fishToxBasis matching are missing.")
}

fish <- fread(fish_file)
chem <- fread(chem_file)

if (!all(c("SiteNumber", "CollectionDate") %in% names(fish))) {
  stop("fish file missing SiteNumber or CollectionDate")
}
if (!all(c("SiteNumber", "CollectionDate") %in% names(chem))) {
  stop("chem file missing SiteNumber or CollectionDate")
}

fish[, CollectionDate := as.Date(CollectionDate)]
chem[, CollectionDate := as.Date(CollectionDate)]
fish[, SiteNumber := standardize_site(SiteNumber)]
chem[, SiteNumber := standardize_site(SiteNumber)]

events <- unique(fish[, .(SiteNumber, CollectionDate)])
chem_cols <- grep("^chem_", names(chem), value = TRUE)
if (!length(chem_cols)) {
  stop("No chem_ columns found in chem data")
}

setkey(chem, SiteNumber, CollectionDate)
chem_event_means <- events[, {
  lo <- CollectionDate - buff_days
  hi <- CollectionDate + buff_days
  sub <- chem[.(SiteNumber)][CollectionDate >= lo & CollectionDate <= hi]
  n_rows <- nrow(sub)
  if (n_rows > 0) {
    vals <- sub[, lapply(.SD, function(x) mean(as.numeric(x), na.rm = TRUE)), .SDcols = chem_cols]
    vals[, n_chem_samples_in_window := n_rows]
    vals
  } else {
    NULL
  }
}, by = .(SiteNumber, CollectionDate)]

fwrite(chem_event_means, out_chem_means)

size_files <- Sys.glob("fish_event_size_metrics_*.csv")
if (!length(size_files)) {
  stop("No fish_event_size_metrics_*.csv found; run script 01 first.")
}
latest_size_file <- size_files[which.max(file.info(size_files)$mtime)]
fish_size <- fread(latest_size_file)
fish_size[, CollectionDate := as.Date(CollectionDate)]
fish_size[, SiteNumber := standardize_site(SiteNumber)]

fish_size_plus_chem <- merge(
  fish_size,
  chem_event_means,
  by = c("SiteNumber", "CollectionDate"),
  all.x = TRUE
)
setorder(fish_size_plus_chem, SiteNumber, CollectionDate)
fwrite(fish_size_plus_chem, out_merged)

qc_lines <- c(
  sprintf("date_tag: %s", today_tag),
  sprintf("events_total_in_fish_file: %d", nrow(events)),
  sprintf("events_with_chem_window: %d", nrow(chem_event_means)),
  sprintf("chem_cols_detected: %d", length(chem_cols)),
  sprintf("fish_size_metrics_file_used: %s", latest_size_file),
  "n_chem_samples_in_window_summary:"
)
qc_lines <- c(qc_lines, capture.output(print(summary(chem_event_means$n_chem_samples_in_window))))

# TU event summary
if (file.exists(existing_tu_file)) {
  tu_dt <- fread(existing_tu_file)
  tu_dt[, CollectionDate := as.Date(CollectionDate)]
  tu_dt[, SiteNumber := standardize_site(SiteNumber)]
  species_tu_cols <- setdiff(names(tu_dt), c("SiteNumber", "CollectionDate"))

  merged_tu <- merge(
    fish_size_plus_chem,
    tu_dt,
    by = c("SiteNumber", "CollectionDate"),
    all.x = TRUE
  )

  if (length(species_tu_cols)) {
    tu_mat <- as.matrix(merged_tu[, ..species_tu_cols])
    suppressWarnings(storage.mode(tu_mat) <- "numeric")
    merged_tu[, TU_event_mean := rowMeans(tu_mat, na.rm = TRUE)]
    merged_tu[, TU_event_median := apply(tu_mat, 1, function(z) {
      z <- z[is.finite(z)]
      if (!length(z)) NA_real_ else median(z)
    })]
    merged_tu[!is.finite(TU_event_mean), TU_event_mean := NA_real_]
  } else {
    merged_tu[, `:=`(TU_event_mean = NA_real_, TU_event_median = NA_real_)]
  }

  fwrite(merged_tu, out_tu)
  qc_lines <- c(qc_lines,
                sprintf("TU_mode: existing_file (%s)", existing_tu_file),
                sprintf("species_tu_cols: %d", length(species_tu_cols)),
                "TU_event_mean_summary:",
                capture.output(print(summary(merged_tu$TU_event_mean))))

} else if (file.exists(standartox_file) && file.exists(cas_key_file)) {
  std <- fread(standartox_file)
  key <- fread(cas_key_file)

  needed <- c("qualifier", "endpoint", "duration_unit", "duration", "concentration", "concentration_unit")
  missing_needed <- setdiff(needed, names(std))
  if (length(missing_needed)) {
    stop(sprintf("standartox file missing columns: %s", paste(missing_needed, collapse = ", ")))
  }

  std[, duration_num := as.numeric(duration)]
  std <- std[
    qualifier == "=" & endpoint %in% c("LC50", "LC50*") &
      duration_unit == "h" & duration_num == 96
  ]

  if ("group" %in% names(std)) {
    std <- std[grepl("fish", as.character(group), ignore.case = TRUE)]
  }

  std[, cas_use := if ("cas" %in% names(std)) as.character(cas) else as.character(casnr)]
  std[, concentration_std_ugL := convert_to_ugL(as.numeric(concentration), concentration_unit)]
  std <- std[is.finite(concentration_std_ugL) & concentration_std_ugL > 0 & !is.na(cas_use) & nzchar(cas_use)]

  tox_bench <- std[, .(tox_benchmark_ugL = median(concentration_std_ugL, na.rm = TRUE)), by = cas_use]
  setnames(tox_bench, "cas_use", "cas")

  key[, CharacteristicName := trimws(as.character(CharacteristicName))]
  key[, cas := trimws(as.character(cas))]
  map_dt <- data.table(
    chem_col = chem_cols,
    CharacteristicName = sub("^chem_", "", chem_cols)
  )
  map_dt <- merge(map_dt, unique(key[, .(CharacteristicName, cas)]), by = "CharacteristicName", all.x = TRUE)
  map_dt <- merge(map_dt, tox_bench, by = "cas", all.x = TRUE)

  usable_map <- map_dt[!is.na(cas) & is.finite(tox_benchmark_ugL) & tox_benchmark_ugL > 0]
  mapped_cols <- usable_map$chem_col

  out_dt <- copy(fish_size_plus_chem)
  if (length(mapped_cols)) {
    tu_matrix <- sapply(seq_len(nrow(usable_map)), function(i) {
      cc <- usable_map$chem_col[i]
      bench <- usable_map$tox_benchmark_ugL[i]
      as.numeric(out_dt[[cc]]) / bench
    })
    if (is.vector(tu_matrix)) tu_matrix <- matrix(tu_matrix, ncol = 1)

    out_dt[, TU_total_event := rowSums(tu_matrix, na.rm = TRUE)]
    out_dt[!is.finite(TU_total_event), TU_total_event := NA_real_]

    mapped_non_na <- out_dt[, sapply(mapped_cols, function(cc) as.integer(is.finite(as.numeric(get(cc)))))]
    if (is.vector(mapped_non_na)) mapped_non_na <- matrix(mapped_non_na, ncol = 1)
    out_dt[, n_chems_mapped := rowSums(mapped_non_na, na.rm = TRUE)]
    out_dt[, prop_chems_mapped := if (length(mapped_cols)) n_chems_mapped / length(mapped_cols) else NA_real_]
  } else {
    out_dt[, `:=`(TU_total_event = NA_real_, n_chems_mapped = 0L, prop_chems_mapped = NA_real_)]
  }

  fwrite(out_dt, out_tu)
  qc_lines <- c(qc_lines,
                sprintf("TU_mode: from_scratch (%s)", standartox_file),
                sprintf("chem_cols_mapped_to_cas: %d", sum(!is.na(map_dt$cas))),
                sprintf("chem_cols_mapped_to_tox_benchmark: %d", nrow(usable_map)),
                "TU_total_event_summary:",
                capture.output(print(summary(out_dt$TU_total_event))))
} else {
  fwrite(fish_size_plus_chem, out_tu)
  qc_lines <- c(qc_lines,
                "TU_mode: skipped (no existing TU file and no standartox/key available)")
}

writeLines(qc_lines, out_qc)
cat(sprintf("Wrote: %s\n", out_chem_means))
cat(sprintf("Wrote: %s\n", out_merged))
cat(sprintf("Wrote: %s\n", out_tu))
cat(sprintf("Wrote: %s\n", out_qc))
