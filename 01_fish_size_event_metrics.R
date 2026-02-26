suppressWarnings({
  library(data.table)
})

today_tag <- format(Sys.Date(), "%Y%m%d")
input_file <- "fish_size.csv"
output_file <- sprintf("fish_event_size_metrics_%s.csv", today_tag)

if (!file.exists(input_file)) {
  stop(sprintf("Required input file not found: %s", input_file))
}

weighted_quantile <- function(x, w, probs = 0.5) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  if (!length(x)) {
    return(rep(NA_real_, length(probs)))
  }
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w)
  total <- sum(w)
  if (!is.finite(total) || total <= 0) {
    return(rep(NA_real_, length(probs)))
  }
  sapply(probs, function(p) {
    p2 <- min(max(p, 0), 1)
    idx <- which(cw >= p2 * total)[1]
    x[idx]
  })
}

standardize_site <- function(x) {
  x_chr <- trimws(as.character(x))
  already_usgs <- grepl("^USGS-", x_chr, ignore.case = TRUE)
  out <- x_chr

  numericish <- grepl("^[0-9]+$", x_chr)
  need_convert <- !already_usgs & numericish
  if (any(need_convert)) {
    out[need_convert] <- sprintf("USGS-%08d", as.integer(x_chr[need_convert]))
  }

  # Normalize USGS prefix casing and pad numeric portion if needed
  usgs_now <- grepl("^USGS-", out, ignore.case = TRUE)
  if (any(usgs_now)) {
    suffix <- sub("^USGS-", "", toupper(out[usgs_now]))
    suffix_digits <- gsub("[^0-9]", "", suffix)
    pad_ok <- nchar(suffix_digits) > 0
    suffix_final <- suffix
    suffix_final[pad_ok] <- sprintf("%08d", as.integer(suffix_digits[pad_ok]))
    out[usgs_now] <- paste0("USGS-", suffix_final)
  }

  out
}

fish <- fread(input_file)
required_cols <- c("SiteNumber", "CollectionDate", "Abundance", "TotalLength_mm", "Weight_g")
missing_cols <- setdiff(required_cols, names(fish))
if (length(missing_cols)) {
  stop(sprintf("Missing required columns in %s: %s", input_file, paste(missing_cols, collapse = ", ")))
}

fish[, CollectionDate := as.Date(CollectionDate)]
fish[, SiteNumber := standardize_site(SiteNumber)]
fish[, Abundance := as.numeric(Abundance)]
fish[, TotalLength_mm := as.numeric(TotalLength_mm)]
fish[, Weight_g := as.numeric(Weight_g)]

if ("Species" %in% names(fish)) {
  richness_source <- "Species"
  fish[, richness_taxon := trimws(as.character(Species))]
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

fish[, K := ifelse(
  is.finite(Weight_g) & is.finite(TotalLength_mm) & TotalLength_mm > 0,
  100 * Weight_g / ((TotalLength_mm / 10)^3),
  NA_real_
)]

event_metrics <- fish[, {
  len_ok <- is.finite(TotalLength_mm)
  wt_ok <- is.finite(Weight_g)
  k_ok <- is.finite(K)

  list(
    n_total = sum(ifelse(is.finite(Abundance), Abundance, 0), na.rm = TRUE),
    n_len = sum(ifelse(len_ok & is.finite(Abundance), Abundance, 0), na.rm = TRUE),
    n_wt = sum(ifelse(wt_ok & is.finite(Abundance), Abundance, 0), na.rm = TRUE),
    mean_len = if (any(len_ok & is.finite(Abundance) & Abundance > 0)) {
      weighted.mean(TotalLength_mm[len_ok], Abundance[len_ok], na.rm = TRUE)
    } else NA_real_,
    median_len = weighted_quantile(TotalLength_mm, Abundance, 0.5),
    p90_len = weighted_quantile(TotalLength_mm, Abundance, 0.9),
    max_len = if (any(len_ok)) max(TotalLength_mm[len_ok], na.rm = TRUE) else NA_real_,
    mean_wt = if (any(wt_ok & is.finite(Abundance) & Abundance > 0)) {
      weighted.mean(Weight_g[wt_ok], Abundance[wt_ok], na.rm = TRUE)
    } else NA_real_,
    median_wt = weighted_quantile(Weight_g, Abundance, 0.5),
    p90_wt = weighted_quantile(Weight_g, Abundance, 0.9),
    max_wt = if (any(wt_ok)) max(Weight_g[wt_ok], na.rm = TRUE) else NA_real_,
    mean_K = if (any(k_ok & is.finite(Abundance) & Abundance > 0)) {
      weighted.mean(K[k_ok], Abundance[k_ok], na.rm = TRUE)
    } else NA_real_,
    median_K = weighted_quantile(K, Abundance, 0.5),
    p90_K = weighted_quantile(K, Abundance, 0.9),
    richness_event = uniqueN(richness_taxon[!is.na(richness_taxon) & nzchar(richness_taxon)])
  )
}, by = .(SiteNumber, CollectionDate)]

setorder(event_metrics, SiteNumber, CollectionDate)
fwrite(event_metrics, output_file)

cat(sprintf("Wrote %s\n", output_file))
cat(sprintf("QC: events=%d, sites=%d, date_range=[%s,%s], richness_source=%s\n",
            nrow(event_metrics), uniqueN(event_metrics$SiteNumber),
            as.character(min(event_metrics$CollectionDate, na.rm = TRUE)),
            as.character(max(event_metrics$CollectionDate, na.rm = TRUE)),
            richness_source))
