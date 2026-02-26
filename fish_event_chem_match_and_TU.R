suppressPackageStartupMessages({
  library(data.table)
})

buff_days <- 90
tox_duration_hours <- 48
endpoint_set <- c("LC50", "EC50", "LC50*", "EC50*", "XX50")
min_required_chems <- 5

pick_col <- function(nms, patterns) {
  for (p in patterns) {
    idx <- grep(p, nms, ignore.case = TRUE)
    if (length(idx) > 0) return(nms[idx[1]])
  }
  NA_character_
}

to_ugL <- function(value, unit) {
  u <- tolower(trimws(as.character(unit)))
  mult <- fifelse(u %in% c("ug/l", "µg/l", "μg/l", "ppb", "ug\u00b7l-1", "ug l-1"), 1,
           fifelse(u %in% c("mg/l", "mg\u00b7l-1", "mg l-1"), 1000,
           fifelse(u %in% c("ng/l", "ng\u00b7l-1", "ng l-1"), 0.001,
           fifelse(u %in% c("g/l", "g\u00b7l-1", "g l-1"), 1e6, NA_real_))))
  as.numeric(value) * mult
}

fish_file <- "fishData_filtered_fishToxBasis_012226.csv"
chem_file <- "chemData_filtered_fishToxBasis_012226.csv"

if (!file.exists(fish_file)) stop(sprintf("Missing required file: %s", fish_file))
if (!file.exists(chem_file)) stop(sprintf("Missing required file: %s", chem_file))

fish_events <- unique(fread(fish_file, select = c("SiteNumber", "CollectionDate")))
fish_events[, CollectionDate := as.Date(CollectionDate)]
setorder(fish_events, SiteNumber, CollectionDate)

chem <- fread(chem_file)
if (!all(c("SiteNumber", "CollectionDate") %in% names(chem))) {
  stop("chemData file must include SiteNumber and CollectionDate columns")
}
chem[, CollectionDate := as.Date(CollectionDate)]
chem_cols <- grep("^chem_", names(chem), value = TRUE)
if (!length(chem_cols)) stop("No chem_* columns found in chemData file")
if (length(chem_cols) < min_required_chems) {
  stop(sprintf("Too few chem columns found (%d < min_required_chems=%d)", length(chem_cols), min_required_chems))
}

matched_list <- vector("list", nrow(fish_events))
for (i in seq_len(nrow(fish_events))) {
  s <- fish_events$SiteNumber[i]
  d <- fish_events$CollectionDate[i]
  in_window <- chem[
    SiteNumber == s & CollectionDate >= (d - buff_days) & CollectionDate <= (d + buff_days)
  ]
  if (nrow(in_window) == 0) next

  row <- in_window[, lapply(.SD, function(x) mean(as.numeric(x), na.rm = TRUE)), .SDcols = chem_cols]
  row[, `:=`(
    SiteNumber = s,
    CollectionDate = d,
    n_chem_samples_in_window = nrow(in_window)
  )]
  setcolorder(row, c("SiteNumber", "CollectionDate", "n_chem_samples_in_window", chem_cols))
  matched_list[[i]] <- row
}

chem_event_means <- rbindlist(matched_list, fill = TRUE)
if (nrow(chem_event_means) == 0) {
  stop("No fish events had chem samples within the specified buffer window")
}

chem_out <- sprintf("chem_event_means_fishToxBasis_%s.csv", format(Sys.Date(), "%Y%m%d"))
fwrite(chem_event_means, chem_out)

size_metrics_files <- Sys.glob("fish_event_size_metrics_*.csv")
if (!length(size_metrics_files)) {
  stop("fish_event_size_metrics output not found. Run fish_size_event_metrics.R first.")
}
size_metrics_file <- size_metrics_files[which.max(file.info(size_metrics_files)$mtime)]
size_metrics <- fread(size_metrics_file)
size_metrics[, CollectionDate := as.Date(CollectionDate)]

fish_size_chem <- merge(
  size_metrics,
  chem_event_means,
  by = c("SiteNumber", "CollectionDate"),
  all = FALSE
)
size_chem_out <- sprintf("fish_size_plus_chem_fishToxBasis_%s.csv", format(Sys.Date(), "%Y%m%d"))
fwrite(fish_size_chem, size_chem_out)

qc_lines <- character()
qc_lines <- c(qc_lines, sprintf("fish events total: %d", nrow(fish_events)))
qc_lines <- c(qc_lines, sprintf("fish events with chem matched: %d", nrow(chem_event_means)))
qc_lines <- c(qc_lines, "n_chem_samples_in_window summary:")
qc_lines <- c(qc_lines, capture.output(print(summary(chem_event_means$n_chem_samples_in_window))))
qc_lines <- c(qc_lines, sprintf("# chem columns available: %d", length(chem_cols)))

tox_files <- Sys.glob("*")
tox_files <- tox_files[grepl("^(standartox|standardtox).*(csv|rds)$", basename(tox_files), ignore.case = TRUE)]

tu_done <- FALSE
if (length(tox_files) > 0) {
  tox_file <- tox_files[which.max(file.info(tox_files)$mtime)]
  tox <- if (grepl("\\.rds$", tox_file, ignore.case = TRUE)) readRDS(tox_file) else fread(tox_file)
  tox <- as.data.table(tox)

  nms <- names(tox)
  qual_col <- pick_col(nms, c("^qualifier$", "resultqualifier", "qual"))
  endpoint_col <- pick_col(nms, c("endpoint"))
  duration_col <- pick_col(nms, c("duration", "exposure"))
  toxval_col <- pick_col(nms, c("tox.*value", "value", "result", "conc", "effect"))
  unit_col <- pick_col(nms, c("unit"))
  cas_col <- pick_col(nms, c("^cas$", "casrn", "cas_number", "chemical.*cas"))
  chem_name_col <- pick_col(nms, c("chemical.*name", "chemname", "substance", "compound", "chemical"))

  if (!is.na(qual_col)) tox <- tox[get(qual_col) == "=" | is.na(get(qual_col))]
  if (!is.na(endpoint_col)) tox <- tox[toupper(trimws(get(endpoint_col))) %in% toupper(endpoint_set)]

  if (!is.na(duration_col)) {
    dur_num <- suppressWarnings(as.numeric(gsub("[^0-9.]+", "", as.character(tox[[duration_col]]))))
    tox <- tox[!is.na(dur_num) & dur_num == tox_duration_hours]
  } else {
    warning("No duration column found in toxicity file; skipping duration filter.")
    qc_lines <- c(qc_lines, "WARNING: No duration column found in toxicity file; duration filter skipped.")
  }

  if (!is.na(toxval_col) && !is.na(unit_col)) {
    tox[, tox_value_ugL := to_ugL(get(toxval_col), get(unit_col))]
    tox <- tox[is.finite(tox_value_ugL) & tox_value_ugL > 0]

    if (!is.na(cas_col)) {
      tox[, chem_id := trimws(tolower(as.character(get(cas_col))))]
    } else if (!is.na(chem_name_col)) {
      tox[, chem_id := trimws(tolower(as.character(get(chem_name_col))))]
    } else {
      tox[, chem_id := NA_character_]
    }
    tox <- tox[!is.na(chem_id) & nzchar(chem_id)]

    bench <- tox[, .(tox_benchmark_ugL = median(tox_value_ugL, na.rm = TRUE)), by = chem_id]

    chem_map <- data.table(
      chem_col = chem_cols,
      chem_token = tolower(sub("^chem_", "", chem_cols))
    )

    if (!is.na(cas_col)) {
      chem_map[, chem_id := fifelse(chem_token %in% bench$chem_id, chem_token, NA_character_)]
    } else {
      chem_map[, chem_id := NA_character_]
    }

    if (!is.na(chem_name_col)) {
      name_bench <- unique(tox[, .(chem_name_id = trimws(tolower(as.character(get(chem_name_col))))), by = chem_id])
      setnames(name_bench, "chem_name_id", "chem_token")
      chem_map <- merge(chem_map, name_bench, by = "chem_token", all.x = TRUE)
      chem_map[is.na(chem_id.x), chem_id.x := chem_id.y]
      chem_map[, `:=`(chem_id = chem_id.x, chem_id.x = NULL, chem_id.y = NULL)]
    }

    chem_map <- merge(chem_map, bench, by = "chem_id", all.x = TRUE)
    mapped <- chem_map[!is.na(tox_benchmark_ugL)]

    qc_lines <- c(qc_lines, sprintf("# chem columns mapped to tox benchmarks: %d", nrow(mapped)))

    if (nrow(mapped) > 0) {
      tu_dt <- copy(fish_size_chem)
      for (j in seq_len(nrow(mapped))) {
        cc <- mapped$chem_col[j]
        b <- mapped$tox_benchmark_ugL[j]
        tu_name <- paste0("TU_", cc)
        tu_dt[, (tu_name) := as.numeric(get(cc)) / b]
      }
      tu_cols <- grep("^TU_chem_", names(tu_dt), value = TRUE)
      tu_dt[, TU_total_event := rowSums(.SD, na.rm = TRUE), .SDcols = tu_cols]

      tu_out <- sprintf("fish_size_plus_TU_fishToxBasis_%s.csv", format(Sys.Date(), "%Y%m%d"))
      fwrite(tu_dt, tu_out)
      tu_done <- TRUE

      qc_lines <- c(qc_lines, "TU_total_event summary:")
      qc_lines <- c(qc_lines, capture.output(print(summary(tu_dt$TU_total_event))))
      qc_lines <- c(qc_lines, sprintf("TU output file: %s", tu_out))
    } else {
      qc_lines <- c(qc_lines, "No chem columns could be mapped to toxicity benchmarks; TU not computed.")
    }
  } else {
    qc_lines <- c(qc_lines, "Toxicity file found but missing tox value/unit columns; TU not computed.")
  }
} else {
  qc_lines <- c(qc_lines, "No toxicity file matching standartox|standardtox(.csv|.rds) found; TU calculation skipped.")
}

qc_out <- sprintf("fish_TU_QC_%s.txt", format(Sys.Date(), "%Y%m%d"))
writeLines(qc_lines, qc_out)

cat(paste(qc_lines, collapse = "\n"), "\n")
message(sprintf("Wrote: %s", chem_out))
message(sprintf("Wrote: %s", size_chem_out))
message(sprintf("Wrote: %s", qc_out))
if (!tu_done) message("TU output was not written.")
