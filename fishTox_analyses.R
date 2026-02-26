
# set openblas to use 36 threads, openMP to 1 to avoid conflicts
library(RhpcBLASctl)
RhpcBLASctl::blas_set_num_threads(parallel::detectCores() - 4)


library(tidyverse)
library(glmmTMB)

# read in density and tox data
fish <- data.table::fread("../data/fish_speciesWithToxicUnits_012326.csv") %>%
  as.data.frame() %>%
  pivot_longer(cols = 23:ncol(.), names_to = "Species", values_to = "Density")

tox <- data.table::fread("../data/fish_toxicUnits_012326.csv") %>%
  as.data.frame() %>%
  pivot_longer(cols = 3:ncol(.), names_to = "Species", values_to = "ToxicUnits")

# fish taxonomy
fish_tax <- data.table::fread("../data/FishTaxonomy_42024.csv") %>%
  filter(Species_period %in% unique(fish$Species))


# join density, taxonomic, and tox data
dat <- fish  %>%
  left_join(tox, by = c("SiteNumber", "CollectionDate", "Species")) %>%
  left_join(fish_tax %>% mutate(Species = Species_period) %>% select(-Species_period), by = "Species")

dat$CollectionDate_scaled <- scale(as.numeric(dat$CollectionDate))
dat$ToxicUnits_scaled     <- scale(dat$ToxicUnits)

dat$SiteNumber <- factor(dat$SiteNumber)
dat$Family     <- factor(dat$Family)
dat$Species    <- factor(dat$Species)

sp_totals <- dat %>%
  group_by(Species) %>%
  summarise(total_density = sum(Density))

dat <- dat %>%
  filter(Species %in% sp_totals$Species[sp_totals$total_density > 0])

# total density per sample
sample_totals <- dat %>%
  group_by(SiteNumber, CollectionDate) %>%
  summarise(total_density = sum(Density), .groups = "drop")

# keep only samples with at least one individual
keep_samples <- sample_totals %>%
  filter(total_density > 0) %>%
  select(SiteNumber, CollectionDate)

dat_clean <- dat %>%
  semi_join(keep_samples, by = c("SiteNumber", "CollectionDate"))

# Calculate total density per species; retain only top 50 most abundant species
top_species <- dat_clean %>%
  filter(!Family %in% c("Salmonidae", "Poeciliidae")) %>% # drop fams w few species
  group_by(Species) %>%
  summarise(total_density = sum(Density)) %>%
  arrange(desc(total_density)) %>%
  slice(1:50) %>%
  pull(Species)

dat_subset <- dat_clean %>%
  filter(Species %in% top_species) %>%
  droplevels() %>%
  mutate(logDensity = log(Density + 0.1),
         logTU = log(ToxicUnits))


# read in array of risk quotients
rq_array <- qs::qread("fish_riskQuotientArray_012326.qs")

# convert array to long form for analysis and filter based on `data_subset` from TU analyses
rq_long <- as.data.frame.table(
  rq_array,
  responseName = "RiskQuotient"
)  %>%
  left_join(fish_tax %>% mutate(Species = Species_period) %>% select(Species, Family), by = "Species") %>%
  filter(Species %in% dat_subset$Species,
         SamplingEvent %in% paste(dat_subset$SiteNumber, dat_subset$CollectionDate, sep = "|"))

# rejoin sampling event info
rq_long$SiteNumber <- sapply(as.character(rq_long$SamplingEvent), \(x) str_split_1(x,"\\|")[1])
rq_long$CollectionDate <- sapply(as.character(rq_long$SamplingEvent), \(x) str_split_1(x,"\\|")[2])

rq_long2 <- rq_long %>%
  left_join(dat_subset %>% dplyr::select(SiteNumber, 
                                         CollectionDate, 
                                         CollectionDate_scaled, 
                                         Density, 
                                         logDensity, 
                                         Species,
                                         ToxicUnits,
                                         ToxicUnits_scaled,
                                         logTU) %>%
              mutate(CollectionDate = as.character(CollectionDate)),
            by = c("CollectionDate", "SiteNumber", "Species"))

# remove zeros for sites where species was never detected
sites <- unique(rq_long2$SiteNumber)

zeroesToRemove <- list()

for (i in seq_along(sites)) {
  
  site <- sites[i]
  spp <- unique(rq_long2$Species[rq_long2$SiteNumber == site])
  
  for (j in seq_along(spp)) {
    
    if (all(rq_long2$Density[rq_long2$Species == spp[j] & rq_long2$SiteNumber == site] == 0)) {
      
      zeroesToRemove[[length(zeroesToRemove) + 1]] <- rownames(rq_long2[rq_long2$Species == spp[j] & rq_long2$SiteNumber == site,])
      
    }
    
  }
  
}

rowsToRemove <- as.numeric(unique(do.call(c, zeroesToRemove)))


rq_long3 <- rq_long2[-rowsToRemove,] %>%
  group_by(SiteNumber, CollectionDate) %>%
  mutate(rq_mean = mean(RiskQuotient),
         rq_resid = RiskQuotient - rq_mean) %>%
  ungroup() %>%
  mutate(rq_resid_scaled = scale(rq_resid))

summary(rq_long3$rq_resid_scaled)


mform <- as.formula("Density ~ ToxicUnits_scaled*Family + CollectionDate_scaled + (1|Species) + (0 + rq_resid_scaled|Species/Chemical)")

m2 <- glmmTMB(mform, 
              data = rq_long3,
              family = tweedie)

summary(m2)

re_m2 <- ranef(m2, condVar = T) %>% as.data.frame()

library(emmeans)

famTox <- emtrends(m2, ~Family, var = "ToxicUnits_scaled") %>% as.data.frame() %>% 
  rename("est" = "ToxicUnits_scaled.trend",
         "ci.lower" = "asymp.LCL",
         "ci.upper" = "asymp.UCL") %>%
  mutate(tox_resid = est - mean(est), # get residual toxicity to each family by subtracting mean toxicity
         tox_resid.lower = ci.lower - mean(est),
         tox_resid.upper = ci.upper - mean(est))

fams <- c("Centrarchidae", "Cyprinidae", "Ictaluridae")

famTox$Family <- factor(famTox$Family, levels = rev(fams))

famColors <- c("#E69F00", "#56B4E9", "#009E73")
famShapes <- c(21:23)

# fig1 shows relative effects of mixture toxicity (TU) on families 
# (i.e. deviations from overall toxicity, aka the main effect of toxic units in the model)
fig1 <- ggplot(famTox, aes(x = Family, y = tox_resid, ymax = tox_resid.upper, 
                         ymin = tox_resid.lower, color = Family, shape = Family, fill = Family)) +
  geom_point() +
  geom_errorbar(width = 0.25) +
  scale_color_manual(values = rev(famColors)) +
  scale_fill_manual(values = rev(famColors)) +
  scale_shape_manual(values = rev(famShapes)) +
  theme(legend.position = "top", legend) +
  theme_bw() +
  geom_hline(yintercept = 0, lty = 2) +
  theme(text = element_text(size = 14)) +
  ylab("Relative Effect of Mixture Toxicity [as ln(Toxic Units)] on Fish Density")

png("../figures/FamilyToxDeclines.png", units = "in", res = 1200, width = 8.5, height = 7)

fig1

dev.off()


speciesChemEffects <- re_m2 %>%
  filter(grpvar == "Chemical:Species")

speciesChemEffects$Species <- sapply(as.character(speciesChemEffects$grp), \(x) str_split_1(x,":")[2])
speciesChemEffects$Chemical <- sapply(as.character(speciesChemEffects$grp), \(x) str_split_1(x,":")[1])

speciesChemEffects2 <- speciesChemEffects %>%
  left_join(fish_tax %>% mutate(Species = Species_period) %>% select(Species, Family), by = "Species") %>%
  mutate(ci.lower = condval - 1.96*condsd,
         ci.upper = condval + 1.96*condsd, 
         Chemical = str_replace(Chemical, "chem_", ""),
         Species = str_replace(Species, "\\.", " "))

speciesChemEffects2$Family <- factor(speciesChemEffects2$Family, levels = fams)
speciesChemEffects2$Species <- factor(speciesChemEffects2$Species, levels = unique(speciesChemEffects2$Species[order(speciesChemEffects2$Family)]))



# fig2 plots relative effects of each chemical on each species (i.e. when this chemical dominates the mixture, there are more density declines)
fig2 <- ggplot(speciesChemEffects2, aes(x = condval, y = Species, 
                              xmin = ci.lower, xmax = ci.upper,
                              color = Family,
                              shape = Family,
                              fill = Family)) +
  facet_wrap(~Chemical, nrow = 2) +
  geom_point() +
  geom_errorbar() +
  theme_bw() +
  geom_vline(xintercept = 0, lty = 2) +
  xlab("Relative Effect of Chemical on Population Density Controlling for Overall Mixture Toxicity") +
  scale_color_manual(values = famColors) +
  scale_fill_manual(values = famColors) +
  scale_shape_manual(values = famShapes) +
  theme(legend.position = "top", legend)


png("../figures/ChemicalRelativeToxicity.png", units = "in", res = 1200, width = 13, height = 11)

fig2

dev.off()

# generate predictions of family responses to gradient of toxic units
# controlling for overall change over time and residual risk quotient (relative contribution of chemicals to mixture toxicity)
rqResidZeroLoc <- pnorm(0, mean = mean(rq_long3$rq_resid), sd = sd(rq_long3$rq_resid), 
                     lower.tail = F)
rqResidZero <- qnorm(rqResidZeroLoc, mean = mean(rq_long3$rq_resid_scaled), sd = sd(rq_long3$rq_resid_scaled), 
                  lower.tail = F)


newdat <- data.frame(
  Family = rep(fams, each = 1000),
  ToxicUnits_scaled = seq(min(rq_long3$ToxicUnits_scaled), max(rq_long3$ToxicUnits_scaled), by = (max(rq_long3$ToxicUnits_scaled) - min(rq_long3$ToxicUnits_scaled))/999),
  ToxicUnits = seq(min(rq_long3$ToxicUnits), max(rq_long3$ToxicUnits), by = (max(rq_long3$ToxicUnits) - min(rq_long3$ToxicUnits))/999),
  CollectionDate_scaled = mean(rq_long3$CollectionDate_scaled),
  rq_resid_scaled = rqResidZero,
  Species = NA,
  Chemical = NA
)

newdat$Family <- factor(newdat$Family, levels = fams)

preds <- predict(m2, newdata = newdat, se.fit = T, allow.new.levels = T, type = "link")

newdat$est_link <- preds$fit
newdat$ci.lower_link <- preds$fit - 1.96*preds$se.fit
newdat$ci.upper_link <- preds$fit + 1.96*preds$se.fit

# transform CI to conform with tweedie family
inv_link <- family(m2)$linkinv

newdat$est <- inv_link(newdat$est_link)
newdat$ci.lower <- inv_link(newdat$ci.lower_link)
newdat$ci.upper <- inv_link(newdat$ci.upper_link)


summary(newdat[,c("est", "ci.lower", "ci.upper", "ToxicUnits")])


# fig3 plots predicted responses of each family to increasing toxic units
# useful for estimating at what level of toxicity each family is likely to be extirpated
# cyprinids hit zero first, followed by ictalurids then centrarchids
fig3 <- ggplot(newdat, aes(y = est, ymax = ci.upper, ymin = ci.lower, 
                         color = Family, fill = Family, linetype = Family,
                         x = ToxicUnits)) +
  geom_line() + geom_ribbon(alpha = 0.2, color = "transparent") +
  theme_bw() + xlim(0.00025,0.01) + ylim(0, 0.051) +
  scale_color_manual(values = famColors) +
  scale_fill_manual(values = famColors) +
  ylab("Density") +
  xlab("Toxic Units")

png("../figures/FamilyTrendToxicUnits.png", units = "in", res = 1200, width = 6, height = 5)

fig3

dev.off()

# get all unique risk quotient values for each chemical at each sampling event
# RQs are based on family-level toxicity so we will need to slice by family to
# avoid pseudoreplicating
rq_summary <- rq_long3 %>%
  group_by(SiteNumber, CollectionDate, Family, Chemical) %>% 
  slice(1) %>% ungroup() %>%
  mutate(Chemical = str_replace(Chemical, "chem_", ""))

# plots for each family, note that the black lines (dense collection of points) generally
# represent risk quotients calculated for non-detects, i.e. where we assumed the concentration
# is equal to 1/2 detection limit

# cyprinids in general had lower risk quotients than the other two families
cypRQplot <- ggplot(rq_summary[rq_summary$Family == "Cyprinidae",], aes(x = Chemical, y = RiskQuotient)) +
  facet_wrap(~Family) +
  geom_point(position = position_jitter()) +
  theme_bw() + ylab("Risk Quotient") + xlab("") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.00001))

cypRQplot

# centrarchids seem to have a fewer chemicals in general contributing to estimated toxicity
# i.e. the distribution of risk quotients seems to favor a smaller handful of chems
cenRQplot <- ggplot(rq_summary[rq_summary$Family == "Centrarchidae",], aes(x = Chemical, y = RiskQuotient)) +
  facet_wrap(~Family) +
  geom_point(position = position_jitter()) +
  theme_bw() + ylab("") + xlab("") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.00001))

cenRQplot

# ictalurids toxicity is mostly attributed to one chemical: dieldrin with minor contributions 
# from a handlful of other chems
ictRQplot <- ggplot(rq_summary[rq_summary$Family == "Ictaluridae",], aes(x = Chemical, y = RiskQuotient)) +
  facet_wrap(~Family) +
  geom_point(position = position_jitter()) +
  theme_bw() + ylab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.00001))

ictRQplot


# fig 4 shows relative contribution of each chemical to site mixture toxicity
# using a scatterplot of the risk quotients at each sampling event
fig4 <- gridExtra::grid.arrange(cenRQplot, cypRQplot, ictRQplot)


png("../figures/ChemRQsByFamily.png", units = "in", res = 1200, width = 5, height = 9)

gridExtra::grid.arrange(cenRQplot, cypRQplot, ictRQplot)

dev.off()




