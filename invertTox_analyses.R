# TEMPLATE: macroinvertebrate toxicity analyses
# Adapted from `fishTox_analyses.R` with explicit placeholders.

library(RhpcBLASctl)
RhpcBLASctl::blas_set_num_threads(parallel::detectCores() - 4)

library(tidyverse)
library(glmmTMB)

# --- INPUTS (PLACEHOLDER: update filenames to your final outputs) ---
invert_abundance_file <- "../data/invert_taxaWithToxicUnits_TEMPLATE.csv"
invert_tu_file <- "../data/invert_toxicUnits_TEMPLATE.csv"
invert_rq_array_file <- "invert_riskQuotientArray_TEMPLATE.qs"

# PLACEHOLDER: invert taxonomy crosswalk used to define grouping variable(s)
invert_taxonomy_file <- "../data/INVERT_TAXONOMY_CROSSWALK_PLACEHOLDER.csv"

# --- OUTPUT FIGURES (PLACEHOLDER: adjust names as desired) ---
fig_family_slope <- "../figures/InvertGroupToxDeclines_TEMPLATE.png"
fig_taxon_chem <- "../figures/InvertChemicalRelativeToxicity_TEMPLATE.png"
fig_group_trend <- "../figures/InvertGroupTrendToxicUnits_TEMPLATE.png"
fig_rq_scatter <- "../figures/InvertChemRQsByGroup_TEMPLATE.png"

# read abundance and TU data
invert <- data.table::fread(invert_abundance_file) %>%
  as.data.frame() %>%
  pivot_longer(cols = 23:ncol(.), names_to = "Taxon", values_to = "Density")

tox <- data.table::fread(invert_tu_file) %>%
  as.data.frame() %>%
  pivot_longer(cols = 3:ncol(.), names_to = "Taxon", values_to = "ToxicUnits")

# load taxonomy crosswalk
invert_tax <- data.table::fread(invert_taxonomy_file)

# ----------------------------------------------------------------------------
# PLACEHOLDER: choose grouping variable for mixed model
# ----------------------------------------------------------------------------
# Examples: Family, Order, FunctionalFeedingGroup, Habit
# The selected column must exist in invert_tax.
group_var <- "Family"  # PLACEHOLDER

if (!"TaxonLabel" %in% colnames(invert_tax)) {
  stop("Invert taxonomy crosswalk must include `TaxonLabel` for joining.")
}
if (!group_var %in% colnames(invert_tax)) {
  stop("Chosen group_var not found in invert taxonomy crosswalk: ", group_var)
}

# join density, TU, taxonomy
dat <- invert %>%
  left_join(tox, by = c("SiteNumber", "CollectionDate", "Taxon")) %>%
  left_join(invert_tax %>%
              mutate(Taxon = TaxonLabel) %>%
              dplyr::select(Taxon, all_of(group_var)),
            by = "Taxon")

# numeric/scaled predictors
# PLACEHOLDER: if CollectionDate is not directly coercible, parse to Date first.
dat$CollectionDate_scaled <- scale(as.numeric(dat$CollectionDate))
dat$ToxicUnits_scaled <- scale(dat$ToxicUnits)

# factors
dat$SiteNumber <- factor(dat$SiteNumber)
dat[[group_var]] <- factor(dat[[group_var]])
dat$Taxon <- factor(dat$Taxon)

# keep taxa with non-zero totals
taxon_totals <- dat %>%
  group_by(Taxon) %>%
  summarise(total_density = sum(Density, na.rm = TRUE), .groups = "drop")

dat <- dat %>%
  filter(Taxon %in% taxon_totals$Taxon[taxon_totals$total_density > 0])

# keep sampling events with total abundance > 0
sample_totals <- dat %>%
  group_by(SiteNumber, CollectionDate) %>%
  summarise(total_density = sum(Density, na.rm = TRUE), .groups = "drop")

keep_samples <- sample_totals %>%
  filter(total_density > 0) %>%
  dplyr::select(SiteNumber, CollectionDate)

dat_clean <- dat %>%
  semi_join(keep_samples, by = c("SiteNumber", "CollectionDate"))

# retain top abundant taxa (PLACEHOLDER threshold)
top_n_taxa <- 50

top_taxa <- dat_clean %>%
  group_by(Taxon) %>%
  summarise(total_density = sum(Density, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_density)) %>%
  slice(1:top_n_taxa) %>%
  pull(Taxon)

dat_subset <- dat_clean %>%
  filter(Taxon %in% top_taxa) %>%
  droplevels() %>%
  mutate(logDensity = log(Density + 0.1),
         logTU = log(ToxicUnits))

# read risk quotient array
rq_array <- qs::qread(invert_rq_array_file)

# array -> long
rq_long <- as.data.frame.table(
  rq_array,
  responseName = "RiskQuotient"
)

# Expect array dimnames: Taxon, Chemical, SamplingEvent
# harmonize expected names from as.data.frame.table
colnames(rq_long)[1:3] <- c("Taxon", "Chemical", "SamplingEvent")

rq_long <- rq_long %>%
  filter(Taxon %in% dat_subset$Taxon,
         SamplingEvent %in% paste(dat_subset$SiteNumber, dat_subset$CollectionDate, sep = "|"))

# rejoin event info
rq_long$SiteNumber <- sapply(as.character(rq_long$SamplingEvent), function(x) str_split_1(x, "\\|")[1])
rq_long$CollectionDate <- sapply(as.character(rq_long$SamplingEvent), function(x) str_split_1(x, "\\|")[2])

rq_long2 <- rq_long %>%
  left_join(dat_subset %>%
              dplyr::select(SiteNumber,
                            CollectionDate,
                            CollectionDate_scaled,
                            Density,
                            logDensity,
                            Taxon,
                            ToxicUnits,
                            ToxicUnits_scaled,
                            logTU,
                            all_of(group_var)) %>%
              mutate(CollectionDate = as.character(CollectionDate)),
            by = c("CollectionDate", "SiteNumber", "Taxon"))

# remove zeros for site/taxon combos where taxon is never detected at that site
sites <- unique(rq_long2$SiteNumber)
zeroesToRemove <- list()

for (i in seq_along(sites)) {
  site <- sites[i]
  taxa <- unique(rq_long2$Taxon[rq_long2$SiteNumber == site])

  for (j in seq_along(taxa)) {
    if (all(rq_long2$Density[rq_long2$Taxon == taxa[j] & rq_long2$SiteNumber == site] == 0, na.rm = TRUE)) {
      zeroesToRemove[[length(zeroesToRemove) + 1]] <-
        rownames(rq_long2[rq_long2$Taxon == taxa[j] & rq_long2$SiteNumber == site, ])
    }
  }
}

rowsToRemove <- suppressWarnings(as.numeric(unique(do.call(c, zeroesToRemove))))
rowsToRemove <- rowsToRemove[!is.na(rowsToRemove)]

rq_long3 <- rq_long2
if (length(rowsToRemove) > 0) {
  rq_long3 <- rq_long2[-rowsToRemove, ]
}

rq_long3 <- rq_long3 %>%
  group_by(SiteNumber, CollectionDate) %>%
  mutate(rq_mean = mean(RiskQuotient, na.rm = TRUE),
         rq_resid = RiskQuotient - rq_mean) %>%
  ungroup() %>%
  mutate(rq_resid_scaled = scale(rq_resid))

summary(rq_long3$rq_resid_scaled)

# ----------------------------------------------------------------------------
# Mixed model (PLACEHOLDER: update grouping structure if needed)
# ----------------------------------------------------------------------------
# Mirrors fish model:
#   Density ~ ToxicUnits_scaled*Group + CollectionDate_scaled +
#   (1|Taxon) + (0 + rq_resid_scaled|Taxon/Chemical)
mform <- as.formula(
  paste0("Density ~ ToxicUnits_scaled*", group_var,
         " + CollectionDate_scaled + (1|Taxon) + (0 + rq_resid_scaled|Taxon/Chemical)")
)

m2 <- glmmTMB(
  mform,
  data = rq_long3,
  family = tweedie
)

summary(m2)
re_m2 <- ranef(m2, condVar = TRUE) %>% as.data.frame()

library(emmeans)

# group-level toxicity slopes
grp_tox <- emtrends(m2, as.formula(paste0("~", group_var)), var = "ToxicUnits_scaled") %>%
  as.data.frame() %>%
  rename(
    est = ToxicUnits_scaled.trend,
    ci.lower = asymp.LCL,
    ci.upper = asymp.UCL
  ) %>%
  mutate(
    tox_resid = est - mean(est),
    tox_resid.lower = ci.lower - mean(est),
    tox_resid.upper = ci.upper - mean(est)
  )

# plotting helpers
groups <- sort(unique(as.character(grp_tox[[group_var]])))
grp_tox[[group_var]] <- factor(grp_tox[[group_var]], levels = rev(groups))

palette_cols <- scales::hue_pal()(max(length(groups), 3))
shapes <- rep(21:25, length.out = length(groups))

# fig1: relative group responses to TU
fig1 <- ggplot(
  grp_tox,
  aes_string(x = group_var, y = "tox_resid", ymax = "tox_resid.upper", ymin = "tox_resid.lower",
             color = group_var, shape = group_var, fill = group_var)
) +
  geom_point() +
  geom_errorbar(width = 0.25) +
  scale_color_manual(values = palette_cols) +
  scale_fill_manual(values = palette_cols) +
  scale_shape_manual(values = shapes) +
  theme_bw() +
  geom_hline(yintercept = 0, lty = 2) +
  theme(text = element_text(size = 14), legend.position = "top") +
  ylab("Relative Effect of Mixture Toxicity [as ln(Toxic Units)] on Macroinvertebrate Density")

png(fig_family_slope, units = "in", res = 1200, width = 8.5, height = 7)
fig1
dev.off()

# taxon-by-chemical random effects
taxonChemEffects <- re_m2 %>%
  filter(grpvar == "Chemical:Taxon")

taxonChemEffects$Taxon <- sapply(as.character(taxonChemEffects$grp), function(x) str_split_1(x, ":")[2])
taxonChemEffects$Chemical <- sapply(as.character(taxonChemEffects$grp), function(x) str_split_1(x, ":")[1])

taxonChemEffects2 <- taxonChemEffects %>%
  left_join(invert_tax %>% mutate(Taxon = TaxonLabel) %>% dplyr::select(Taxon, all_of(group_var)), by = "Taxon") %>%
  mutate(
    ci.lower = condval - 1.96 * condsd,
    ci.upper = condval + 1.96 * condsd,
    Chemical = str_replace(Chemical, "chem_", "")
  )

taxonChemEffects2[[group_var]] <- factor(taxonChemEffects2[[group_var]], levels = groups)
taxonChemEffects2$Taxon <- factor(
  taxonChemEffects2$Taxon,
  levels = unique(taxonChemEffects2$Taxon[order(taxonChemEffects2[[group_var]])])
)

fig2 <- ggplot(
  taxonChemEffects2,
  aes_string(x = "condval", y = "Taxon", xmin = "ci.lower", xmax = "ci.upper",
             color = group_var, shape = group_var, fill = group_var)
) +
  facet_wrap(~Chemical, nrow = 2) +
  geom_point() +
  geom_errorbar() +
  theme_bw() +
  geom_vline(xintercept = 0, lty = 2) +
  xlab("Relative Effect of Chemical on Taxon Density (controlling for overall mixture toxicity)") +
  scale_color_manual(values = palette_cols) +
  scale_fill_manual(values = palette_cols) +
  scale_shape_manual(values = shapes) +
  theme(legend.position = "top")

png(fig_taxon_chem, units = "in", res = 1200, width = 13, height = 11)
fig2
dev.off()

# predictions along toxic-unit gradient
rqResidZeroLoc <- pnorm(
  0,
  mean = mean(rq_long3$rq_resid, na.rm = TRUE),
  sd = sd(rq_long3$rq_resid, na.rm = TRUE),
  lower.tail = FALSE
)
rqResidZero <- qnorm(
  rqResidZeroLoc,
  mean = mean(rq_long3$rq_resid_scaled, na.rm = TRUE),
  sd = sd(rq_long3$rq_resid_scaled, na.rm = TRUE),
  lower.tail = FALSE
)

newdat <- data.frame(
  grp = rep(groups, each = 1000),
  ToxicUnits_scaled = seq(min(rq_long3$ToxicUnits_scaled, na.rm = TRUE),
                          max(rq_long3$ToxicUnits_scaled, na.rm = TRUE),
                          length.out = 1000),
  ToxicUnits = seq(min(rq_long3$ToxicUnits, na.rm = TRUE),
                   max(rq_long3$ToxicUnits, na.rm = TRUE),
                   length.out = 1000),
  CollectionDate_scaled = mean(rq_long3$CollectionDate_scaled, na.rm = TRUE),
  rq_resid_scaled = rqResidZero,
  Taxon = NA,
  Chemical = NA
)

newdat[[group_var]] <- factor(newdat$grp, levels = groups)
newdat$grp <- NULL

preds <- predict(m2, newdata = newdat, se.fit = TRUE, allow.new.levels = TRUE, type = "link")

newdat$est_link <- preds$fit
newdat$ci.lower_link <- preds$fit - 1.96 * preds$se.fit
newdat$ci.upper_link <- preds$fit + 1.96 * preds$se.fit

inv_link <- family(m2)$linkinv

newdat$est <- inv_link(newdat$est_link)
newdat$ci.lower <- inv_link(newdat$ci.lower_link)
newdat$ci.upper <- inv_link(newdat$ci.upper_link)

fig3 <- ggplot(
  newdat,
  aes_string(y = "est", ymax = "ci.upper", ymin = "ci.lower",
             color = group_var, fill = group_var, linetype = group_var,
             x = "ToxicUnits")
) +
  geom_line() +
  geom_ribbon(alpha = 0.2, color = "transparent") +
  theme_bw() +
  scale_color_manual(values = palette_cols) +
  scale_fill_manual(values = palette_cols) +
  ylab("Density") +
  xlab("Toxic Units")

png(fig_group_trend, units = "in", res = 1200, width = 6, height = 5)
fig3
dev.off()

# RQ summary by event/group/chemical
rq_summary <- rq_long3 %>%
  group_by(SiteNumber, CollectionDate, .data[[group_var]], Chemical) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(Chemical = str_replace(Chemical, "chem_", ""))

# simple faceted RQ scatter by group
fig4 <- ggplot(rq_summary,
               aes(x = Chemical, y = RiskQuotient, color = .data[[group_var]])) +
  facet_wrap(as.formula(paste("~", group_var)), scales = "free_y") +
  geom_point(position = position_jitter(width = 0.2, height = 0), alpha = 0.7) +
  theme_bw() +
  xlab("") +
  ylab("Risk Quotient") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")

png(fig_rq_scatter, units = "in", res = 1200, width = 7, height = 8)
fig4
dev.off()
