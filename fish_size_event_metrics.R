suppressPackageStartupMessages({
  library(data.table)
})

weighted_quantile <- function(x, w, probs = c(0.5), na.rm = TRUE) {
  if (length(x) != length(w)) stop("x and w must be same length")
  if (na.rm) {
    keep <- !(is.na(x) | is.na(w))
    x <- x[keep]
    w <- w[keep]
  }
  if (!length(x)) return(rep(NA_real_, length(probs)))
  w <- as.numeric(w)
  x <- as.numeric(x)
  w[is.na(w) | w < 0] <- 0
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  if (!length(x) || sum(w) <= 0) return(rep(NA_real_, length(probs)))

  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)

  vapply(probs, function(p) {
    p <- max(min(p, 1), 0)
    idx <- which(cw >= p)[1]
    if (is.na(idx)) x[length(x)] else x[idx]
  }, numeric(1))
}

infile <- "fish_size.csv"
if (!file.exists(infile)) {
  stop(sprintf("Required input file not found: %s", infile))
}

fish_size <- fread(infile)
required_cols <- c("SiteNumber", "CollectionDate", "Abundance", "TotalLength_mm", "Weight_g")
missing_cols <- setdiff(required_cols, names(fish_size))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required columns in %s: %s", infile, paste(missing_cols, collapse = ", ")))
}

fish_size[, CollectionDate := as.Date(CollectionDate)]
fish_size <- fish_size[!is.na(SiteNumber) & !is.na(CollectionDate)]

fish_size[, Abundance := suppressWarnings(as.numeric(Abundance))]
fish_size[, Abundance := fifelse(is.na(Abundance), 0, Abundance)]
fish_size[, TotalLength_mm := suppressWarnings(as.numeric(TotalLength_mm))]
fish_size[, Weight_g := suppressWarnings(as.numeric(Weight_g))]

if ("Species" %in% names(fish_size)) {
  fish_size[, taxon_for_richness := as.character(Species)]
} else if ("BenchTaxonName" %in% names(fish_size)) {
  fish_size[, taxon_for_richness := as.character(BenchTaxonName)]
} else if (all(c("Genus", "Species") %in% names(fish_size))) {
  fish_size[, taxon_for_richness := trimws(paste(Genus, Species))]
} else {
  fish_size[, taxon_for_richness := NA_character_]
}

fish_size[
  Weight_g > 0 & TotalLength_mm > 0,
  K := 100 * Weight_g / ((TotalLength_mm / 10)^3)
]

event_metrics <- fish_size[, {
  len_valid <- !is.na(TotalLength_mm) & TotalLength_mm > 0
  wt_valid <- !is.na(Weight_g) & Weight_g > 0
  k_valid <- !is.na(K) & is.finite(K) & K > 0

  len_q <- weighted_quantile(TotalLength_mm[len_valid], Abundance[len_valid], probs = c(0.5, 0.9))
  wt_q <- weighted_quantile(Weight_g[wt_valid], Abundance[wt_valid], probs = c(0.5, 0.9))
  k_q  <- weighted_quantile(K[k_valid], Abundance[k_valid], probs = c(0.5, 0.9))

  list(
    n_total = sum(Abundance, na.rm = TRUE),
    n_len = sum(Abundance[len_valid], na.rm = TRUE),
    n_wt = sum(Abundance[wt_valid], na.rm = TRUE),

    mean_len = if (sum(Abundance[len_valid], na.rm = TRUE) > 0) {
      weighted.mean(TotalLength_mm[len_valid], w = Abundance[len_valid], na.rm = TRUE)
    } else NA_real_,
    median_len = len_q[1],
    p90_len = len_q[2],
    max_len = if (any(len_valid)) max(TotalLength_mm[len_valid], na.rm = TRUE) else NA_real_,

    mean_wt = if (sum(Abundance[wt_valid], na.rm = TRUE) > 0) {
      weighted.mean(Weight_g[wt_valid], w = Abundance[wt_valid], na.rm = TRUE)
    } else NA_real_,
    median_wt = wt_q[1],
    p90_wt = wt_q[2],
    max_wt = if (any(wt_valid)) max(Weight_g[wt_valid], na.rm = TRUE) else NA_real_,

    mean_K = if (sum(Abundance[k_valid], na.rm = TRUE) > 0) {
      weighted.mean(K[k_valid], w = Abundance[k_valid], na.rm = TRUE)
    } else NA_real_,
    median_K = k_q[1],
    p90_K = k_q[2],

    richness_event = uniqueN(taxon_for_richness[Abundance > 0 & !is.na(taxon_for_richness) & nzchar(taxon_for_richness)])
  )
}, by = .(SiteNumber, CollectionDate)]

outfile <- sprintf("fish_event_size_metrics_%s.csv", format(Sys.Date(), "%Y%m%d"))
fwrite(event_metrics, outfile)

message(sprintf("Wrote event-level fish size metrics: %s", outfile))
