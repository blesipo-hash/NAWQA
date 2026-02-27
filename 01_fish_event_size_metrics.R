suppressPackageStartupMessages({
  library(data.table)
})

today_tag  <- format(Sys.Date(), "%Y%m%d")
input_file <- "SizeData/fish_size.csv"

output_data_file <- sprintf("output/data/fish_event_size_metrics_%s.csv", today_tag)
qc_file          <- sprintf("output/results/fish_event_size_metrics_QC_%s.txt", today_tag)

dir.create("output/data", recursive = TRUE, showWarnings = FALSE)
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop(sprintf("Required input file not found: %s", input_file))
}

# ---- helpers ----

weighted_quantile <- function(x, w, probs = c(0.5), na.rm = TRUE) {
  if (length(x) != length(w)) stop("x and w must be the same length")
  x <- suppressWarnings(as.numeric(x))
  w <- suppressWarnings(as.numeric(w))

  if (na.rm) {
    keep <- is.finite(x) & is.finite(w) & w > 0
    x <- x[keep]
    w <- w[keep]
  }
  if (!length(x)) return(rep(NA_real_, length(probs)))

  ord <- order(x)
  x <- x[ord]
  w <- w[ord]

  cw <- cumsum(w)
  total <- sum(w)
  if (!is.finite(total) || total <= 0) return(rep(NA_real_, length(probs)))

  vapply(probs, function(p) {
    p2 <- min(max(p, 0), 1)
    idx <- which(cw >= p2 * total)[1]
    if (is.na(idx)) NA_real_ else x[idx]
  }, numeric(1))
}

# Standardize SiteNumber robustly:
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
      substr(digs, nchar(digs) - 7, nchar(digs)),  # last 8 digits
      sprintf("%08s", digs)                        # pad left if shorter
    )
    digs8 <- gsub(" ", "0", digs8)
    out[keep] <- digs8
  }

  out <- ifelse(!is.na(out), paste0("USGS-", out), NA_character_)
  out
}

# ---- read + validate ----

fish <- fread(input_file)

required_cols <- c("SiteNumber", "CollectionDate", "Abundance", "TotalLength_mm", "Weight_g")
missing_cols <- setdiff(required_cols, names(fish))
if (length(missing_cols)) {
  stop(sprintf("Missing required columns in %s: %s", input_file, paste(missing_cols, collapse = ", ")))
}

# Parse + standardize keys
fish[, CollectionDate := as.Date(CollectionDate)]
fish[, SiteNumber := standardize_site(SiteNumber)]

# Keep only valid site IDs + dates
fish <- fish[!is.na(SiteNumber) & grepl("^USGS-\\d{8}$", SiteNumber)]
fish <- fish[!is.na(CollectionDate)]

# Coerce numerics
fish[, Abundance := suppressWarnings(as.numeric(Abundance))]
fish[, TotalLength_mm := suppressWarnings(as.numeric(TotalLength_mm))]
fish[, Weight_g := suppressWarnings(as.numeric(Weight_g))]

# Replace missing/invalid abundance with 0
fish[is.na(Abundance) | Abundance < 0, Abundance := 0]

# Richness taxon (prefer full binomial-like fields)
if ("ScientificName" %in% names(fish)) {
  richness_source <- "ScientificName"
  fish[, richness_taxon := trimws(as.character(ScientificName))]
} else if ("PublishedTaxonName" %in% names(fish)) {
  richness_source <- "PublishedTaxonName"
  fish[, richness_taxon := trimws(as.character(PublishedTaxonName))]
} else if ("BenchTaxonName" %in% names(fish)) {
  richness_source <- "BenchTaxonName"
  fish[, richness_taxon := trimws(as.character(BenchTaxonName))]
} else if (all(c("Genus", "Species") %in% names(fish))) {
  richness_source <- "Genus+Species"
  fish[, richness_taxon := trimws(paste(Genus, Species))]
} else {
  richness_source <- "none"
  fish[, richness_taxon := NA_character_]
}

# Condition factor K (Fulton's K; mm->cm conversion)
fish[Weight_g > 0 & TotalLength_mm > 0,
     K := 100 * Weight_g / ((TotalLength_mm / 10)^3)]
fish[!(Weight_g > 0 & TotalLength_mm > 0), K := NA_real_]

# ---- event-level aggregation ----

event_metrics <- fish[, {
  # validity masks for length, weight, and K (only count fish with positive abundance)
  len_valid <- is.finite(TotalLength_mm) & TotalLength_mm > 0 & is.finite(Abundance) & Abundance > 0
  wt_valid  <- is.finite(Weight_g)       & Weight_g > 0       & is.finite(Abundance) & Abundance > 0
  k_valid   <- is.finite(K)              & K > 0              & is.finite(Abundance) & Abundance > 0

  len_q <- weighted_quantile(TotalLength_mm[len_valid], Abundance[len_valid], probs = c(0.5, 0.9))
  wt_q  <- weighted_quantile(Weight_g[wt_valid],        Abundance[wt_valid],  probs = c(0.5, 0.9))
  k_q   <- weighted_quantile(K[k_valid],                Abundance[k_valid],   probs = c(0.5, 0.9))

  list(
    n_total = sum(Abundance, na.rm = TRUE),
    n_len   = sum(Abundance[len_valid], na.rm = TRUE),
    n_wt    = sum(Abundance[wt_valid],  na.rm = TRUE),

    mean_len   = if (sum(Abundance[len_valid]) > 0) weighted.mean(TotalLength_mm[len_valid], Abundance[len_valid]) else NA_real_,
    median_len = len_q[1],
    p90_len    = len_q[2],
    max_len    = if (any(len_valid)) max(TotalLength_mm[len_valid]) else NA_real_,

    mean_wt    = if (sum(Abundance[wt_valid]) > 0) weighted.mean(Weight_g[wt_valid], Abundance[wt_valid]) else NA_real_,
    median_wt  = wt_q[1],
    p90_wt     = wt_q[2],
    max_wt     = if (any(wt_valid)) max(Weight_g[wt_valid]) else NA_real_,

    mean_K     = if (sum(Abundance[k_valid]) > 0) weighted.mean(K[k_valid], Abundance[k_valid]) else NA_real_,
    median_K   = k_q[1],
    p90_K      = k_q[2],

    richness_event = if (richness_source == "none") NA_integer_
    else uniqueN(richness_taxon[is.finite(Abundance) & Abundance > 0 &
                                  !is.na(richness_taxon) & nzchar(richness_taxon)])
  )
}, by = .(SiteNumber, CollectionDate)]

setorder(event_metrics, SiteNumber, CollectionDate)
fwrite(event_metrics, output_data_file)

# ---- QC ----
len_missing_prop <- fish[, mean(!(is.finite(TotalLength_mm) & TotalLength_mm > 0))]
wt_missing_prop  <- fish[, mean(!(is.finite(Weight_g) & Weight_g > 0))]

qc_lines <- c(
  sprintf("rows_used=%d", nrow(fish)),
  sprintf("events=%d", nrow(event_metrics)),
  sprintf("sites=%d", uniqueN(event_metrics$SiteNumber)),
  sprintf("date_min=%s", as.character(min(event_metrics$CollectionDate, na.rm = TRUE))),
  sprintf("date_max=%s", as.character(max(event_metrics$CollectionDate, na.rm = TRUE))),
  sprintf("richness_source=%s", richness_source),
  sprintf("missing_prop_length=%0.6f", len_missing_prop),
  sprintf("missing_prop_weight=%0.6f", wt_missing_prop)
)
writeLines(qc_lines, qc_file)

cat(sprintf("Wrote %s\n", output_data_file))
cat(sprintf("Wrote %s\n", qc_file))
