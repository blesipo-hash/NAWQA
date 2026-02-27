# =============================================================================
# 04_Q3_within_species_only.R
# Question 3 only: Q3A/Q3B/Q3C within-species size vs event TU
# - Inputs: output/data/fish_species_event_size_metrics_*.csv
#           output/data/event_TU_for_size_events_*.csv
# - Outputs: output/models/Q3A_*.rds, Q3B_*.rds, Q3C_*.rds
#           output/results/Q3_within_species_summary_*.txt
#           output/figures/Q3A_*, Q3B_*, Q3C_* (2 plots each)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("Package 'lme4' is required. Install with install.packages('lme4')", call. = FALSE)
}

today_tag <- format(Sys.Date(), "%Y%m%d")

dir.create("output/data", recursive = TRUE, showWarnings = FALSE)
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/models", recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# helpers
# ---------------------------
find_latest <- function(pattern) {
  files <- Sys.glob(pattern)
  if (!length(files)) {
    stop(sprintf("No files found matching pattern: %s", pattern), call. = FALSE)
  }
  info <- file.info(files)
  files[order(info$mtime, decreasing = TRUE)][1]
}

standardize_site <- function(x) {
  x_chr <- trimws(as.character(x))
  x_chr[x_chr %in% c("", "NA", "NaN")] <- NA_character_
  out <- x_chr

  has_usgs <- !is.na(out) & grepl("^USGS-", out, ignore.case = TRUE)
  if (any(has_usgs)) {
    out[has_usgs] <- gsub("\\D", "", out[has_usgs])
    out[has_usgs][nchar(out[has_usgs]) == 0] <- NA_character_
  }

  numericish <- !is.na(out) & grepl("^[0-9]+$", out)
  other <- !is.na(out) & !has_usgs & !numericish
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

as_num <- function(x) {
  y <- suppressWarnings(as.numeric(x))
  y[is.infinite(y)] <- NA_real_
  y
}

pick_col_optional <- function(dt, candidates) {
  hit <- candidates[candidates %in% names(dt)][1]
  if (is.na(hit)) {
    return(NA_character_)
  }
  hit
}

pick_col_required <- function(dt, candidates, label) {
  hit <- pick_col_optional(dt, candidates)
  if (is.na(hit)) {
    stop(
      sprintf("Missing required %s column. Tried: %s", label, paste(candidates, collapse = ", ")),
      call. = FALSE
    )
  }
  hit
}

# nearest-date TU fill within site (optional)
nearest_fill_tu <- function(spec_dt, tu_dt, max_lag_days = 90L) {
  spec_dt <- copy(spec_dt)
  tu_dt <- copy(tu_dt)

  spec_dt[, row_id__ := .I]
  cand <- merge(
    spec_dt[, .(row_id__, SiteNumber, d_size = CollectionDate)],
    tu_dt[, .(SiteNumber, d_tu = CollectionDate, TU_event, logTU, n_chems_used, n_chem_samples_in_window)],
    by = "SiteNumber",
    allow.cartesian = TRUE
  )

  if (nrow(cand)) {
    cand[, lag_days := abs(as.integer(d_size - d_tu))]
    cand <- cand[lag_days <= max_lag_days]
  }

  nearest <- cand[order(row_id__, lag_days, d_tu)][, .SD[1], by = row_id__]

  out <- merge(
    spec_dt,
    nearest[, .(
      row_id__,
      TU_event_near = TU_event,
      logTU_near = logTU,
      n_chems_used_near = n_chems_used,
      n_chem_samples_in_window_near = n_chem_samples_in_window,
      tu_lag_days = lag_days
    )],
    by = "row_id__",
    all.x = TRUE,
    sort = FALSE
  )

  # fill missing TU values from nearest-date candidates without creating .x/.y column conflicts
  if ("TU_event" %in% names(out)) {
    out[is.na(TU_event), TU_event := TU_event_near]
  } else {
    out[, TU_event := TU_event_near]
  }
  if ("logTU" %in% names(out)) {
    out[is.na(logTU), logTU := logTU_near]
  } else {
    out[, logTU := logTU_near]
  }
  if ("n_chems_used" %in% names(out)) {
    out[is.na(n_chems_used), n_chems_used := n_chems_used_near]
  } else {
    out[, n_chems_used := n_chems_used_near]
  }
  if ("n_chem_samples_in_window" %in% names(out)) {
    out[is.na(n_chem_samples_in_window), n_chem_samples_in_window := n_chem_samples_in_window_near]
  } else {
    out[, n_chem_samples_in_window := n_chem_samples_in_window_near]
  }

  out[, c("TU_event_near", "logTU_near", "n_chems_used_near", "n_chem_samples_in_window_near") := NULL]
  setorder(out, row_id__)
  out[, row_id__ := NULL]
  out
}

# effect CI plot
plot_effect_ci <- function(beta, lo, hi, title, ylab) {
  d <- data.frame(term = "logTU", beta = beta, lo = lo, hi = hi)
  ggplot(d, aes(x = term, y = beta, ymin = lo, ymax = hi)) +
    geom_pointrange() +
    geom_hline(yintercept = 0, linetype = 2) +
    theme_bw() +
    labs(title = title, x = "", y = ylab)
}

# fixed-effect prediction curve for lmer (back-transform outside)
pred_curve_lmer <- function(m, dat, xvar = "logTU", grid_n = 200) {
  x <- dat[[xvar]]
  xr <- quantile(x, probs = c(0.05, 0.95), na.rm = TRUE)
  xseq <- seq(xr[1], xr[2], length.out = grid_n)

  nd <- data.frame(logTU = xseq)

  # build needed covariates
  mf <- model.frame(m)
  if ("date_sc" %in% names(mf)) {
    nd$date_sc <- 0
  }
  if ("log_n_total" %in% names(mf)) {
    nd$log_n_total <- median(mf$log_n_total, na.rm = TRUE)
  }

  # placeholders for random-effect grouping vars (not used when re.form=NA, but needed for predict)
  if ("species_id_model" %in% names(mf)) {
    nd$species_id_model <- mf$species_id_model[1]
  }
  if ("SiteNumber" %in% names(mf)) {
    nd$SiteNumber <- mf$SiteNumber[1]
  }

  fit <- as.numeric(predict(m, newdata = nd, re.form = NA, allow.new.levels = TRUE))

  # Wald SE using fixed-effect vcov
  X <- model.matrix(lme4::nobars(formula(m)), nd)
  b <- lme4::fixef(m)
  V <- as.matrix(vcov(m))

  keep <- intersect(colnames(X), names(b))
  X <- X[, keep, drop = FALSE]
  b <- b[keep]
  V <- V[keep, keep, drop = FALSE]

  eta <- as.numeric(X %*% b)
  se <- sqrt(pmax(0, diag(X %*% V %*% t(X))))

  data.table(logTU = xseq, eta = eta, se = se, fit_link = fit)
}

# ---------------------------
# locate inputs
# ---------------------------
species_file <- find_latest("output/data/fish_species_event_size_metrics_*.csv")
tu_file <- find_latest("output/data/event_TU_for_size_events_*.csv")

cat("Using:\n")
cat("  species-event size:", species_file, "\n")
cat("  event TU          :", tu_file, "\n")

sp <- fread(species_file)
tu <- fread(tu_file)

# ---------------------------
# standardize keys
# ---------------------------
if (!all(c("SiteNumber", "CollectionDate") %in% names(sp))) {
  stop("Species-event size file must contain SiteNumber and CollectionDate.", call. = FALSE)
}
sp[, SiteNumber := standardize_site(SiteNumber)]
sp[, CollectionDate := as.Date(CollectionDate)]
sp <- sp[!is.na(SiteNumber) & grepl("^USGS-\\d{8}$", SiteNumber) & !is.na(CollectionDate)]

# detect species id
species_col <- pick_col_required(
  sp,
  c("species_id", "species_id_model", "Species", "ScientificName", "PublishedTaxonName", "BenchTaxonName"),
  "species id"
)
sp[, species_id_model := as.factor(as.character(get(species_col)))]
sp[is.na(species_id_model) | !nzchar(as.character(species_id_model)), species_id_model := factor("unknown_species")]

# detect outcomes
mean_len_col <- pick_col_required(sp, c("mean_len_sp", "mean_len", "mean_len_model", "MeanLength_mm"), "mean length (species-event)")
p90_len_col <- pick_col_optional(sp, c("p90_len_sp", "p90_len", "p90_len_species"))
max_len_col <- pick_col_optional(sp, c("max_len_sp", "max_len", "max_len_species"))
K_col <- pick_col_optional(sp, c("mean_K_sp", "mean_K", "median_K_sp", "median_K", "K_mean_sp"))

# TU file must have SiteNumber + CollectionDate + TU_event + logTU
if (!all(c("SiteNumber", "CollectionDate") %in% names(tu))) {
  stop("Event TU file must contain SiteNumber and CollectionDate.", call. = FALSE)
}
tu[, SiteNumber := standardize_site(SiteNumber)]
tu[, CollectionDate := as.Date(CollectionDate)]
tu <- tu[!is.na(SiteNumber) & grepl("^USGS-\\d{8}$", SiteNumber) & !is.na(CollectionDate)]

# Ensure TU columns exist
if (!("logTU" %in% names(tu))) {
  if ("TU_event" %in% names(tu)) {
    tu[, logTU := ifelse(is.finite(as_num(TU_event)), log10(as_num(TU_event) + 1e-8), NA_real_)]
  } else {
    stop("TU file missing both logTU and TU_event.", call. = FALSE)
  }
}
if (!("TU_event" %in% names(tu))) {
  tu[, TU_event := NA_real_]
}

# force numeric
tu[, TU_event := as_num(TU_event)]
tu[, logTU := as_num(logTU)]
if ("n_chems_used" %in% names(tu)) {
  tu[, n_chems_used := as_num(n_chems_used)]
}
if ("n_chem_samples_in_window" %in% names(tu)) {
  tu[, n_chem_samples_in_window := as_num(n_chem_samples_in_window)]
}

tu_keep <- unique(
  tu[, .(
    SiteNumber,
    CollectionDate,
    TU_event,
    logTU,
    n_chems_used = if ("n_chems_used" %in% names(tu)) n_chems_used else NA_real_,
    n_chem_samples_in_window = if ("n_chem_samples_in_window" %in% names(tu)) n_chem_samples_in_window else NA_real_
  )]
)

setkey(tu_keep, SiteNumber, CollectionDate)

# remove TU-like columns from species data so merge never creates logTU.x/logTU.y collisions
tu_like_cols <- intersect(c("TU_event", "logTU", "n_chems_used", "n_chem_samples_in_window"), names(sp))
if (length(tu_like_cols)) {
  sp[, (tu_like_cols) := NULL]
}

# ---------------------------
# merge TU onto species-event data
# ---------------------------
sp2 <- merge(sp, tu_keep, by = c("SiteNumber", "CollectionDate"), all.x = TRUE, sort = FALSE)

# Optional nearest-date fill (within site) for missing TU
use_nearest <- TRUE
max_lag_days <- 90L

if (use_nearest) {
  miss_before <- sp2[is.na(logTU), .N]
  sp2_filled <- nearest_fill_tu(sp2, tu_keep, max_lag_days = max_lag_days)
  miss_after <- sp2_filled[is.na(logTU), .N]
  sp2 <- sp2_filled
} else {
  miss_before <- sp2[is.na(logTU), .N]
  miss_after <- miss_before
}

# ---------------------------
# add covariates
# ---------------------------
sp2[, SiteNumber := as.factor(SiteNumber)]
sp2[, date_num := as.numeric(CollectionDate)]
sp2[, date_sc := as.numeric(scale(date_num))]

# log_n_total: if present, use; else try common alternatives; else omit
ln_total_col <- pick_col_optional(sp2, c("log_n_total", "log_n_total_cov", "log_total_abun", "log_total_abun_MGMS"))
if (!is.na(ln_total_col)) {
  sp2[, log_n_total := as_num(get(ln_total_col))]
} else {
  sp2[, log_n_total := NA_real_]
}

# outcomes -> log scale
sp2[, mean_len_y := as_num(get(mean_len_col))]
sp2[, log_mean_len := ifelse(is.finite(mean_len_y) & mean_len_y > 0, log(mean_len_y), NA_real_)]

if (!is.na(p90_len_col)) {
  sp2[, p90_len_y := as_num(get(p90_len_col))]
  sp2[, log_p90_len := ifelse(is.finite(p90_len_y) & p90_len_y > 0, log(p90_len_y), NA_real_)]
}

if (!is.na(max_len_col)) {
  sp2[, max_len_y := as_num(get(max_len_col))]
  sp2[, log_max_len := ifelse(is.finite(max_len_y) & max_len_y > 0, log(max_len_y), NA_real_)]
}

if (!is.na(K_col)) {
  sp2[, K_y := as_num(get(K_col))]
  sp2[, log_K := ifelse(is.finite(K_y) & K_y > 0, log(K_y), NA_real_)]
}

# ---------------------------
# model fitting (lmer)
# ---------------------------
fit_lmer <- function(dt, response, label, include_ln_total = TRUE) {
  d <- dt[is.finite(get(response)) & is.finite(logTU) & is.finite(date_sc)]
  if (include_ln_total) {
    d <- d[is.finite(log_n_total)]
  }

  if (nrow(d) < 100) {
    stop(sprintf("%s: too few rows after filtering (%d)", label, nrow(d)), call. = FALSE)
  }

  if (include_ln_total) {
    form <- as.formula(
      paste0(response, " ~ logTU + date_sc + I(date_sc^2) + log_n_total + (1|species_id_model) + (1|SiteNumber)")
    )
  } else {
    form <- as.formula(
      paste0(response, " ~ logTU + date_sc + I(date_sc^2) + (1|species_id_model) + (1|SiteNumber)")
    )
  }

  m <- lme4::lmer(form, data = as.data.frame(d))
  list(model = m, data = d, formula = form)
}

include_ln_total <- any(is.finite(sp2$log_n_total))

fits <- list()

# Q3A
fits$Q3A <- fit_lmer(sp2, "log_mean_len", "Q3A (within-species mean length)", include_ln_total = include_ln_total)

# Q3B (p90 preferred, else max)
if (!is.na(p90_len_col)) {
  fits$Q3B <- fit_lmer(sp2, "log_p90_len", "Q3B (within-species p90 length)", include_ln_total = include_ln_total)
} else if (!is.na(max_len_col)) {
  fits$Q3B <- fit_lmer(sp2, "log_max_len", "Q3B (within-species max length)", include_ln_total = include_ln_total)
}

# Q3C
if (!is.na(K_col)) {
  fits$Q3C <- fit_lmer(sp2, "log_K", "Q3C (within-species condition K)", include_ln_total = include_ln_total)
}

# ---------------------------
# outputs: models, plots, summary
# ---------------------------
summary_file <- sprintf("output/results/Q3_within_species_summary_%s.txt", today_tag)
lines <- c(
  "Q3 within-species fish size analysis (Q3A/Q3B/Q3C)",
  sprintf("date_tag=%s", today_tag),
  sprintf("species_event_file=%s", species_file),
  sprintf("event_TU_file=%s", tu_file),
  sprintf("use_nearest=%s max_lag_days=%d", ifelse(use_nearest, "yes", "no"), max_lag_days),
  sprintf("missing_TU_before=%d missing_TU_after=%d", miss_before, miss_after),
  sprintf("include_log_n_total=%s", ifelse(include_ln_total, "yes", "no")),
  ""
)

for (nm in names(fits)) {
  obj <- fits[[nm]]
  m <- obj$model
  d <- as.data.table(obj$data)

  beta <- lme4::fixef(m)["logTU"]
  ci <- suppressMessages(confint(m, method = "Wald"))
  ci_logTU <- ci["logTU", ]

  pct <- (exp(beta) - 1) * 100

  model_path <- sprintf("output/models/%s_lmer_%s.rds", nm, today_tag)
  saveRDS(m, model_path)

  # Plot 1: effect CI
  p1 <- plot_effect_ci(
    beta,
    ci_logTU[1],
    ci_logTU[2],
    title = sprintf("%s: logTU effect (95%% CI)", nm),
    ylab = "Effect on log(outcome)"
  )
  fig1 <- sprintf("output/figures/%s_effect_CI_%s.png", nm, today_tag)
  ggsave(fig1, p1, width = 6.5, height = 4.2, dpi = 320)

  # Plot 2: prediction curve (fixed effects only), back-transformed
  pr <- pred_curve_lmer(m, as.data.frame(d), xvar = "logTU", grid_n = 200)
  pr[, fit_bt := exp(eta)]
  pr[, lo_bt := exp(eta - 1.96 * se)]
  pr[, hi_bt := exp(eta + 1.96 * se)]

  ylab2 <- if (nm == "Q3C") {
    "Predicted K (back-transformed)"
  } else {
    "Predicted length (mm; back-transformed)"
  }

  p2 <- ggplot(pr, aes(x = logTU, y = fit_bt)) +
    geom_ribbon(aes(ymin = lo_bt, ymax = hi_bt), alpha = 0.25) +
    geom_line(linewidth = 1.1) +
    theme_bw() +
    labs(
      title = sprintf("%s: fixed-effect prediction", nm),
      x = "Event toxic units (log10 scale)",
      y = ylab2
    )
  fig2 <- sprintf("output/figures/%s_pred_curve_%s.png", nm, today_tag)
  ggsave(fig2, p2, width = 7.5, height = 4.8, dpi = 320)

  lines <- c(
    lines,
    sprintf("---- %s ----", nm),
    sprintf("n_rows=%d", nrow(d)),
    sprintf("n_sites=%d", uniqueN(d$SiteNumber)),
    sprintf("n_species=%d", uniqueN(d$species_id_model)),
    sprintf("beta_logTU=%.6f", beta),
    sprintf("CI_logTU=[%.6f, %.6f]", ci_logTU[1], ci_logTU[2]),
    sprintf("pct_change_per_1log10TU=%.2f%%", pct),
    sprintf("model_rds=%s", model_path),
    sprintf("plot_effect=%s", fig1),
    sprintf("plot_curve=%s", fig2),
    ""
  )
}

if (!length(fits)) {
  lines <- c(lines, "No models were fit (missing eligible outcome columns for Q3B/Q3C).")
}

writeLines(lines, summary_file)
cat(sprintf("Wrote %s\n", summary_file))
cat("Done. Figures are in output/figures and models in output/models.\n")
