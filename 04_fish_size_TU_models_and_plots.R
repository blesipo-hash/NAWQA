suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(mgcv)
})

today_tag <- format(Sys.Date(), "%Y%m%d")

dir.create("output/data", recursive = TRUE, showWarnings = FALSE)
dir.create("output/results", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# ---- required inputs ----
q1q2_base_file <- "output/data/modeldata_Q1Q2_event_20260226.csv"
q3_base_file   <- "output/data/modeldata_Q3_species_event_20260226.csv"
fish_event_size_file   <- "output/data/fish_event_size_metrics_20260226.csv"
fish_species_size_file <- "output/data/fish_species_event_size_metrics_20260226.csv"

required_files <- c(q1q2_base_file, q3_base_file, fish_event_size_file, fish_species_size_file)
missing_required <- required_files[!file.exists(required_files)]
if (length(missing_required)) {
  stop(sprintf(
    "Missing required Step 3 input file(s): %s",
    paste(missing_required, collapse = ", ")
  ), call. = FALSE)
}

find_tu_file <- function() {
  candidates <- c(
    "data/fish_toxicUnits_012326.csv",
    "../data/fish_toxicUnits_012326.csv",
    "fish_toxicUnits_012326.csv"
  )
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) {
    stop("Could not find TU file 'fish_toxicUnits_012326.csv' in data/ or ../data/.", call. = FALSE)
  }
  hits[1]
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

safe_date <- function(dt) {
  dt[, CollectionDate := as.Date(CollectionDate)]
  dt
}

nearest_date_merge_tu <- function(size_dt, tu_dt, exact_threshold = 0.70, max_lag_days = 90) {
  size_dt <- copy(size_dt)
  tu_dt <- copy(tu_dt)

  size_dt[, CollectionDate := as.Date(CollectionDate)]
  tu_dt[, CollectionDate := as.Date(CollectionDate)]

  tu_keep <- unique(tu_dt[, .(
    SiteNumber,
    CollectionDate,
    TU_event_mean,
    TU_event_median,
    n_species_TU_used,
    logTU
  )])

  merged_exact <- merge(size_dt, tu_keep, by = c("SiteNumber", "CollectionDate"), all.x = TRUE, sort = FALSE)
  merged_exact[, tu_lag_days := as.integer(0)]

  exact_rows <- merged_exact[!is.na(logTU), .N]
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

  size_idx <- copy(size_dt)
  size_idx[, row_id__ := .I]

  cand <- merge(
    size_idx[, .(row_id__, SiteNumber, size_date = CollectionDate)],
    tu_keep[, .(SiteNumber, tu_date = CollectionDate, TU_event_mean, TU_event_median, n_species_TU_used, logTU)],
    by = "SiteNumber",
    allow.cartesian = TRUE
  )

  if (nrow(cand) > 0) {
    cand[, lag_days := abs(as.integer(size_date - tu_date))]
    cand <- cand[lag_days <= max_lag_days]
  }

  nearest <- cand[order(row_id__, lag_days, tu_date)][, .SD[1], by = row_id__]

  merged_nearest <- merge(
    size_idx,
    nearest[, .(row_id__, TU_event_mean, TU_event_median, n_species_TU_used, logTU, tu_lag_days = lag_days)],
    by = "row_id__",
    all.x = TRUE,
    sort = FALSE
  )
  setorder(merged_nearest, row_id__)
  merged_nearest[, row_id__ := NULL]

  newly_filled <- merge(
    merged_exact[, .(SiteNumber, CollectionDate, logTU_exact = logTU)],
    merged_nearest[, .(SiteNumber, CollectionDate, logTU_nearest = logTU)],
    by = c("SiteNumber", "CollectionDate"),
    all = FALSE
  )[is.na(logTU_exact) & !is.na(logTU_nearest), .N]

  nearest_rows_total <- merged_nearest[!is.na(logTU), .N]

  list(
    merged = merged_nearest,
    used_nearest = TRUE,
    exact_match_count = exact_rows,
    exact_match_rate = exact_rate,
    nearest_match_count = newly_filled,
    nearest_total_nonNA = nearest_rows_total
  )
}

get_species_col <- function(dt) {
  cand <- c("species_id_model", "species_id", "Species", "ScientificName")
  hit <- cand[cand %in% names(dt)][1]
  if (is.na(hit)) stop("Could not detect species identifier column.", call. = FALSE)
  hit
}

get_mean_len_col <- function(dt) {
  cand <- c("mean_len_sp", "mean_len", "MeanLength_mm", "LengthMean_mm")
  hit <- cand[cand %in% names(dt)][1]
  if (is.na(hit)) stop("Could not detect mean length column.", call. = FALSE)
  hit
}

get_count_col <- function(dt) {
  cand <- c("n_total", "n", "n_fish", "n_len", "n_individual", "count")
  hit <- cand[cand %in% names(dt)][1]
  if (is.na(hit)) {
    stop(
      sprintf(
        "Could not find an event count column. Expected one of: %s. Available columns include: %s",
        paste(cand, collapse = ", "),
        paste(names(dt), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  hit
}

summarize_decile <- function(dt, y_col) {
  dt <- as.data.table(dt)
  if (!(y_col %in% names(dt))) stop(sprintf("summarize_decile(): missing column %s", y_col), call. = FALSE)

  keep <- dt[is.finite(logTU) & is.finite(get(y_col))]
  if (!nrow(keep)) return(data.table())

  breaks <- unique(quantile(keep$logTU, probs = seq(0, 1, 0.1), na.rm = TRUE, names = FALSE))
  if (length(breaks) < 3) {
    keep[, decile := factor("all")]
  } else {
    keep[, decile := cut(logTU, breaks = breaks, include.lowest = TRUE, ordered_result = TRUE)]
  }

  out <- keep[, .(
    mean_logTU = mean(logTU, na.rm = TRUE),
    mean_y = mean(get(y_col), na.rm = TRUE),
    se_y = sd(get(y_col), na.rm = TRUE) / sqrt(.N),
    n = .N
  ), by = decile][order(mean_logTU)]

  out
}

safe_direction_gam <- function(model, base_dat) {
  base_dat <- as.data.frame(base_dat)
  if (!("logTU" %in% names(base_dat))) return("undetermined")
  x <- base_dat$logTU
  x <- x[is.finite(x)]
  if (length(x) < 10) return("undetermined")

  xseq <- seq(quantile(x, 0.1, na.rm = TRUE), quantile(x, 0.9, na.rm = TRUE), length.out = 60)

  newd <- data.frame(
    logTU = xseq,
    date_num = if ("date_num" %in% names(base_dat)) median(base_dat$date_num, na.rm = TRUE) else 0,
    log_n_total = if ("log_n_total" %in% names(base_dat)) median(base_dat$log_n_total, na.rm = TRUE) else 0
  )

  p <- as.numeric(predict(model, newdata = newd, type = "response"))
  slope <- coef(lm(p ~ xseq))[2]
  if (!is.finite(slope)) return("undetermined")
  if (abs(slope) < 1e-10) return("flat")
  if (slope > 0) "increasing" else "decreasing"
}

make_gam_pred_df <- function(m_gam, dat, grid_n = 200) {
  xseq <- seq(min(dat$logTU, na.rm = TRUE), max(dat$logTU, na.rm = TRUE), length.out = grid_n)
  nd <- data.frame(
    logTU = xseq,
    date_num = median(dat$date_num, na.rm = TRUE),
    log_n_total = median(dat$log_n_total, na.rm = TRUE)
  )
  pr <- predict(m_gam, newdata = nd, se.fit = TRUE, type = "response")
  if (!is.list(pr) || !all(c("fit", "se.fit") %in% names(pr))) {
    stop("predict.gam did not return fit/se.fit as expected.", call. = FALSE)
  }
  out <- cbind(nd, fit = as.numeric(pr$fit), se = as.numeric(pr$se.fit))
  out[order(out$logTU), , drop = FALSE]
}

strip_random_terms <- function(form) {
  ftxt <- paste(deparse(form), collapse = "")
  ftxt <- gsub("\\([^\\)]*\\|[^\\)]*\\)", "", ftxt)
  ftxt <- gsub("\\+\\s*\\+", "+", ftxt)
  ftxt <- gsub("~\\s*\\+", "~", ftxt)
  ftxt <- gsub("\\s+", " ", ftxt)
  as.formula(ftxt)
}

get_q3_coef_table <- function(m_q3, engine) {
  sm <- summary(m_q3)
  if (engine == "lme4") return(sm$coefficients)
  if (is.list(sm$coefficients) && !is.null(sm$coefficients$cond)) return(sm$coefficients$cond)
  stop("Unexpected coefficient structure for Q3 model.", call. = FALSE)
}

get_q3_fixef <- function(m_q3, engine) {
  if (engine == "lme4") return(lme4::fixef(m_q3))
  glmmTMB::fixef(m_q3)$cond
}

get_q3_vcov <- function(m_q3, engine) {
  if (engine == "lme4") return(as.matrix(vcov(m_q3)))
  as.matrix(vcov(m_q3)$cond)
}

# ---- 1) TU summary ----
tu_file <- find_tu_file()
tu <- fread(tu_file)

if (!all(c("SiteNumber", "CollectionDate") %in% names(tu))) {
  stop("TU file must contain SiteNumber and CollectionDate.", call. = FALSE)
}

tu[, SiteNumber := standardize_site(SiteNumber)]
tu <- safe_date(tu)
tu <- tu[!is.na(SiteNumber) & grepl("^USGS-\\d{8}$", SiteNumber) & !is.na(CollectionDate)]

TU_cols <- setdiff(names(tu), c("SiteNumber", "CollectionDate"))
if (!length(TU_cols)) stop("No TU columns detected after excluding SiteNumber and CollectionDate.", call. = FALSE)

for (col in TU_cols) {
  vals <- suppressWarnings(as.numeric(tu[[col]]))
  vals[is.infinite(vals)] <- NA_real_
  set(tu, j = col, value = vals)
}

tu[, TU_event_mean := rowMeans(.SD, na.rm = TRUE), .SDcols = TU_cols]
tu[, TU_event_median := apply(.SD, 1, median, na.rm = TRUE), .SDcols = TU_cols]
tu[, n_species_TU_used := rowSums(!is.na(.SD)), .SDcols = TU_cols]
tu[, TU_event_mean := fifelse(is.nan(TU_event_mean), NA_real_, TU_event_mean)]
tu[, TU_event_median := fifelse(is.nan(TU_event_median), NA_real_, TU_event_median)]
tu[, logTU := log10(TU_event_mean + 1e-8)]

tu_event <- unique(tu[, .(SiteNumber, CollectionDate, TU_event_mean, TU_event_median, n_species_TU_used, logTU)])
setorder(tu_event, SiteNumber, CollectionDate)

tu_event_out <- sprintf("output/data/event_TU_summary_%s.csv", today_tag)
fwrite(tu_event, tu_event_out)

# ---- 2) merge TU into Q1/Q2 ----
q1q2 <- fread(q1q2_base_file)
if (!all(c("SiteNumber", "CollectionDate") %in% names(q1q2))) {
  stop("Q1/Q2 modeldata must contain SiteNumber and CollectionDate.", call. = FALSE)
}
q1q2[, SiteNumber := standardize_site(SiteNumber)]
q1q2 <- safe_date(q1q2)
q1q2 <- q1q2[!is.na(SiteNumber) & grepl("^USGS-\\d{8}$", SiteNumber) & !is.na(CollectionDate)]

merge_q1q2 <- nearest_date_merge_tu(q1q2, tu_event, exact_threshold = 0.70, max_lag_days = 90)
q1q2_tu <- merge_q1q2$merged
q1q2_tu[, date_num := as.numeric(CollectionDate)]

count_col <- get_count_col(q1q2_tu)
q1q2_tu[, log_n_total := log(pmax(get(count_col), 0) + 1)]

if (!("p90_len" %in% names(q1q2_tu))) stop("Q1/Q2 data missing p90_len.", call. = FALSE)
if (!("mean_K" %in% names(q1q2_tu))) stop("Q1/Q2 data missing mean_K.", call. = FALSE)

q1q2_tu_out <- sprintf("output/data/modeldata_Q1Q2_event_TU_%s.csv", today_tag)
setorder(q1q2_tu, SiteNumber, CollectionDate)
fwrite(q1q2_tu, q1q2_tu_out)

# ---- 3) merge TU into Q3 ----
q3 <- fread(q3_base_file)
if (!all(c("SiteNumber", "CollectionDate") %in% names(q3))) {
  stop("Q3 modeldata must contain SiteNumber and CollectionDate.", call. = FALSE)
}
q3[, SiteNumber := standardize_site(SiteNumber)]
q3 <- safe_date(q3)
q3 <- q3[!is.na(SiteNumber) & grepl("^USGS-\\d{8}$", SiteNumber) & !is.na(CollectionDate)]

q3_tu <- merge(
  q3,
  tu_event,
  by = c("SiteNumber", "CollectionDate"),
  all.x = TRUE,
  sort = FALSE
)

species_col <- get_species_col(q3_tu)
len_col <- get_mean_len_col(q3_tu)

q3_tu[, species_id_model := as.character(get(species_col))]
q3_tu[is.na(species_id_model) | !nzchar(species_id_model), species_id_model := "unknown_species"]

q3_tu[, mean_len_model := suppressWarnings(as.numeric(get(len_col)))]
q3_tu[, date_num := as.numeric(CollectionDate)]
q3_tu[, log_mean_len := ifelse(is.finite(mean_len_model) & mean_len_model > 0, log(mean_len_model), NA_real_)]

if (!("log_n_total" %in% names(q3_tu))) {
  event_nt <- unique(q1q2_tu[, .(SiteNumber, CollectionDate, log_n_total)])
  q3_tu <- merge(q3_tu, event_nt, by = c("SiteNumber", "CollectionDate"), all.x = TRUE, sort = FALSE)
}
if (!("log_n_total" %in% names(q3_tu))) stop("Could not populate log_n_total in Q3 data.", call. = FALSE)

q3_tu_out <- sprintf("output/data/modeldata_Q3_species_event_TU_%s.csv", today_tag)
setorder(q3_tu, SiteNumber, CollectionDate)
fwrite(q3_tu, q3_tu_out)

# ---- 4) fit models ----
q1_dat <- q1q2_tu[is.finite(p90_len) & is.finite(logTU) & is.finite(date_num) & is.finite(log_n_total)]
q2_dat <- q1q2_tu[is.finite(mean_K) & is.finite(logTU) & is.finite(date_num) & is.finite(log_n_total)]
q3_dat <- q3_tu[is.finite(log_mean_len) & is.finite(logTU) & is.finite(date_num) & is.finite(log_n_total)]

if (nrow(q1_dat) < 30) stop("Insufficient rows for Q1 model after filtering.", call. = FALSE)
if (nrow(q2_dat) < 30) stop("Insufficient rows for Q2 model after filtering.", call. = FALSE)
if (nrow(q3_dat) < 30) stop("Insufficient rows for Q3 model after filtering.", call. = FALSE)

m_q1 <- gam(p90_len ~ s(logTU) + s(date_num) + log_n_total, data = q1_dat, method = "REML")
m_q2 <- gam(mean_K  ~ s(logTU) + s(date_num) + log_n_total, data = q2_dat, method = "REML")

lmer_available <- requireNamespace("lme4", quietly = TRUE)
glmmTMB_available <- requireNamespace("glmmTMB", quietly = TRUE)

q3_dat[, species_id_model := as.factor(species_id_model)]
q3_dat[, SiteNumber := as.factor(SiteNumber)]
q3_dat[, date_num_sc := as.numeric(scale(date_num))]

q3_formula <- log_mean_len ~ logTU + date_num_sc + I(date_num_sc^2) + log_n_total +
  (1 | species_id_model) + (1 | SiteNumber)

if (lmer_available) {
  m_q3 <- lme4::lmer(q3_formula, data = q3_dat)
  q3_model_engine <- "lme4"
} else if (glmmTMB_available) {
  m_q3 <- glmmTMB::glmmTMB(q3_formula, data = q3_dat)
  q3_model_engine <- "glmmTMB"
} else {
  stop("Neither lme4 nor glmmTMB is available; cannot fit Q3 mixed model.", call. = FALSE)
}

# ---- 5) plots (2 per question) ----
fig_q1a <- sprintf("output/figures/Q1a_p90_len_vs_logTU_%s.png", today_tag)
fig_q1b <- sprintf("output/figures/Q1b_p90_len_deciles_%s.png", today_tag)
fig_q2a <- sprintf("output/figures/Q2a_meanK_vs_logTU_%s.png", today_tag)
fig_q2b <- sprintf("output/figures/Q2b_meanK_deciles_%s.png", today_tag)
fig_q3a <- sprintf("output/figures/Q3a_pred_length_vs_logTU_%s.png", today_tag)
fig_q3b <- sprintf("output/figures/Q3b_species_random_effects_%s.png", today_tag)

plot_q1 <- make_gam_pred_df(m_q1, q1_dat, grid_n = 200)
p_q1a <- ggplot(q1_dat, aes(x = logTU, y = p90_len)) +
  geom_point(alpha = 0.35, size = 1.2) +
  geom_ribbon(data = plot_q1, aes(x = logTU, ymin = fit - 1.96 * se, ymax = fit + 1.96 * se), inherit.aes = FALSE, alpha = 0.25, fill = "steelblue") +
  geom_line(data = plot_q1, aes(x = logTU, y = fit), inherit.aes = FALSE, linewidth = 1.1, color = "steelblue4") +
  theme_bw() +
  labs(x = "Event toxic units (log10 scale)", y = "90th percentile fish length (mm)", title = "Q1a: p90 length vs TU")
ggsave(fig_q1a, p_q1a, width = 8, height = 5, dpi = 320)

q1_dec <- summarize_decile(q1_dat, "p90_len")
p_q1b <- ggplot(q1_dec, aes(x = mean_logTU, y = mean_y)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_y - 1.96 * se_y, ymax = mean_y + 1.96 * se_y), width = 0.03) +
  theme_bw() +
  labs(x = "Event toxic units (log10 decile mean)", y = "Mean p90 length (mm)", title = "Q1b: p90 length by TU deciles")
ggsave(fig_q1b, p_q1b, width = 8, height = 5, dpi = 320)

plot_q2 <- make_gam_pred_df(m_q2, q2_dat, grid_n = 200)
p_q2a <- ggplot(q2_dat, aes(x = logTU, y = mean_K)) +
  geom_point(alpha = 0.35, size = 1.2) +
  geom_ribbon(data = plot_q2, aes(x = logTU, ymin = fit - 1.96 * se, ymax = fit + 1.96 * se), inherit.aes = FALSE, alpha = 0.25, fill = "darkorange") +
  geom_line(data = plot_q2, aes(x = logTU, y = fit), inherit.aes = FALSE, linewidth = 1.1, color = "darkorange4") +
  theme_bw() +
  labs(x = "Event toxic units (log10 scale)", y = "Mean Fulton condition factor (K)", title = "Q2a: mean K vs TU")
ggsave(fig_q2a, p_q2a, width = 8, height = 5, dpi = 320)

q2_dec <- summarize_decile(q2_dat, "mean_K")
p_q2b <- ggplot(q2_dec, aes(x = mean_logTU, y = mean_y)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_y - 1.96 * se_y, ymax = mean_y + 1.96 * se_y), width = 0.03) +
  theme_bw() +
  labs(x = "Event toxic units (log10 decile mean)", y = "Mean K", title = "Q2b: mean K by TU deciles")
ggsave(fig_q2b, p_q2b, width = 8, height = 5, dpi = 320)

coef_tab_q3 <- get_q3_coef_table(m_q3, q3_model_engine)
if (!("logTU" %in% rownames(coef_tab_q3))) stop("Q3 model summary does not contain logTU coefficient.", call. = FALSE)

xseq_q3 <- seq(min(q3_dat$logTU, na.rm = TRUE), max(q3_dat$logTU, na.rm = TRUE), length.out = 200)
new_q3 <- data.frame(
  logTU = xseq_q3,
  date_num_sc = 0,
  log_n_total = median(q3_dat$log_n_total, na.rm = TRUE),
  species_id_model = q3_dat$species_id_model[1],
  SiteNumber = q3_dat$SiteNumber[1]
)

fixed_form <- strip_random_terms(q3_formula)
X <- model.matrix(fixed_form, new_q3)
b <- get_q3_fixef(m_q3, q3_model_engine)
V <- get_q3_vcov(m_q3, q3_model_engine)

keep <- intersect(colnames(X), names(b))
X <- X[, keep, drop = FALSE]
b <- b[keep]
V <- V[keep, keep, drop = FALSE]

eta <- as.numeric(X %*% b)
se  <- sqrt(pmax(0, diag(X %*% V %*% t(X))))

plot_q3a <- data.table(logTU = xseq_q3, fit = eta, se = se)
plot_q3a[, fit_len := exp(fit)]
plot_q3a[, lo_len  := exp(fit - 1.96 * se)]
plot_q3a[, hi_len  := exp(fit + 1.96 * se)]

p_q3a <- ggplot(plot_q3a, aes(x = logTU, y = fit_len)) +
  geom_ribbon(aes(ymin = lo_len, ymax = hi_len), alpha = 0.25, fill = "purple") +
  geom_line(linewidth = 1.1, color = "purple4") +
  theme_bw() +
  labs(x = "Event toxic units (log10 scale)", y = "Predicted mean fish length (mm)", title = "Q3a: fixed-effect TU relationship (back-transformed)")
ggsave(fig_q3a, p_q3a, width = 8, height = 5, dpi = 320)

if (q3_model_engine == "lme4") {
  re_obj <- lme4::ranef(m_q3, condVar = TRUE)$species_id_model
  re_species <- as.data.table(re_obj, keep.rownames = "species_id")
  setnames(re_species, old = names(re_species)[2], new = "effect")
  pv <- attr(re_obj, "postVar")
  re_species[, se := if (!is.null(pv)) sqrt(pv[1, 1, ]) else NA_real_]
} else {
  re <- glmmTMB::ranef(m_q3)$cond$species_id_model
  re_species <- as.data.table(re, keep.rownames = "species_id")
  setnames(re_species, old = names(re_species)[2], new = "effect")
  re_species[, se := NA_real_]
}

re_species[, abs_effect := abs(effect)]
re_top <- re_species[order(-abs_effect)][1:min(30, .N)]
re_top[, species_id := factor(species_id, levels = re_top[order(effect)]$species_id)]

p_q3b <- ggplot(re_top, aes(x = species_id, y = effect)) +
  geom_point() +
  geom_errorbar(aes(ymin = effect - 1.96 * se, ymax = effect + 1.96 * se), width = 0.2, na.rm = TRUE) +
  coord_flip() +
  theme_bw() +
  labs(x = "Species (top 30 by |random intercept|)", y = "Random intercept estimate", title = "Q3b: species random effects")
ggsave(fig_q3b, p_q3b, width = 9, height = 7, dpi = 320)

# ---- 6) results summary text ----
q1_sm <- summary(m_q1)$s.table
q2_sm <- summary(m_q2)$s.table

get_smooth_p <- function(s_table, term = "s(logTU)") {
  if (is.null(s_table) || nrow(s_table) == 0) return(NA_real_)
  if (term %in% rownames(s_table)) return(as.numeric(s_table[term, "p-value"]))
  idx <- grep("logTU", rownames(s_table))
  if (length(idx)) return(as.numeric(s_table[idx[1], "p-value"]))
  NA_real_
}

q1_p  <- get_smooth_p(q1_sm, "s(logTU)")
q2_p  <- get_smooth_p(q2_sm, "s(logTU)")
q1_dir <- safe_direction_gam(m_q1, q1_dat)
q2_dir <- safe_direction_gam(m_q2, q2_dat)

coef_tab_q3 <- get_q3_coef_table(m_q3, q3_model_engine)
q3_est <- coef_tab_q3["logTU", "Estimate"]
q3_se  <- coef_tab_q3["logTU", "Std. Error"]
pcol   <- intersect(colnames(coef_tab_q3), c("Pr(>|t|)", "Pr(>|z|)"))
q3_p   <- if (length(pcol)) coef_tab_q3["logTU", pcol[1]] else NA_real_

result_file <- sprintf("output/results/fish_size_TU_results_%s.txt", today_tag)

res_lines <- c(
  "Fish size + TU modeling summary",
  sprintf("date_tag=%s", today_tag),
  "",
  sprintf("TU file used: %s", tu_file),
  sprintf("n_events_q1q2_after_merge=%d", nrow(q1q2_tu)),
  sprintf("n_species_events_q3_after_merge=%d", nrow(q3_tu)),
  sprintf("q1q2_exact_match_count=%d", merge_q1q2$exact_match_count),
  sprintf("q1q2_exact_match_rate=%.4f", merge_q1q2$exact_match_rate),
  sprintf("q1q2_nearest_newly_filled=%d", merge_q1q2$nearest_match_count),
  "",
  "TU summary (event mean):",
  paste(capture.output(summary(tu_event$TU_event_mean)), collapse = "\n"),
  "",
  "n_species_TU_used distribution:",
  paste(capture.output(summary(tu_event$n_species_TU_used)), collapse = "\n"),
  "",
  sprintf("Q1 smooth p-value (s(logTU))=%s", ifelse(is.na(q1_p), "NA", format(q1_p, digits = 6))),
  sprintf("Q1 direction=%s", q1_dir),
  "- Q1 interpretation:",
  sprintf("  * TU smooth term is %s at conventional thresholds (p=%s).", ifelse(is.finite(q1_p) && q1_p < 0.05, "statistically significant", "not statistically significant"), ifelse(is.na(q1_p), "NA", format(q1_p, digits = 3))),
  sprintf("  * Fitted TU relationship is %s over central TU range.", q1_dir),
  "  * Controls (time trend, log_n_total) are included in the GAM.",
  "",
  sprintf("Q2 smooth p-value (s(logTU))=%s", ifelse(is.na(q2_p), "NA", format(q2_p, digits = 6))),
  sprintf("Q2 direction=%s", q2_dir),
  "- Q2 interpretation:",
  sprintf("  * TU smooth term is %s at conventional thresholds (p=%s).", ifelse(is.finite(q2_p) && q2_p < 0.05, "statistically significant", "not statistically significant"), ifelse(is.na(q2_p), "NA", format(q2_p, digits = 3))),
  sprintf("  * Fitted TU relationship is %s over central TU range.", q2_dir),
  "  * Mean K is modeled with the same controls for comparability with Q1.",
  "",
  sprintf("Q3 engine=%s", q3_model_engine),
  sprintf("Q3 logTU estimate=%.6g", q3_est),
  sprintf("Q3 logTU SE=%.6g", q3_se),
  sprintf("Q3 logTU p-value=%s", ifelse(is.na(q3_p), "NA", format(q3_p, digits = 5))),
  "- Q3 interpretation:",
  sprintf("  * The logTU fixed effect on log-mean length is %.4f (SE %.4f).", q3_est, q3_se),
  "  * Positive => larger within-species mean length with increasing TU; negative => truncation.",
  "  * Random intercepts account for baseline differences among species and sites.",
  "",
  "Figure outputs:",
  fig_q1a,
  fig_q1b,
  fig_q2a,
  fig_q2b,
  fig_q3a,
  fig_q3b
)

writeLines(res_lines, result_file)

cat(sprintf("Wrote %s\n", tu_event_out))
cat(sprintf("Wrote %s\n", q1q2_tu_out))
cat(sprintf("Wrote %s\n", q3_tu_out))
cat(sprintf("Wrote %s\n", result_file))
cat("Wrote figures in output/figures\n")
