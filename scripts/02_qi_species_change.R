# =============================================================================
# 02_qi_species_change.R
# Methods: Question i
#   - "Defining site-level species changes"
#   - "Defining initial commonness and rarity"
#
# Input : data/ac_raw_data.csv, data/SupplementaryDataS4_390species.csv
# Output: output/intermediate/qi_data_complete.rds           zero-filled counts
#         output/intermediate/qi_species_change.csv          one row per site x species
#         output/intermediate/qi_species_change_initial_ra.csv
#                                                            rows with a finite initial
#                                                            relative abundance + 10% bin
# =============================================================================

source("R/config.R")
source("R/functions_species_change.R")
library(dplyr)
library(tidyr)

require_files(path_data(FILES$qi_raw_counts), path_data(FILES$link_sauer))

# ---- 1. load counts: January-February counts, 1990-2018 -----------------------
raw  <- read.csv(path_data(FILES$qi_raw_counts))
data <- raw %>%
  filter(Season.x %in% SEASON_KEEP, Year >= YEAR_MIN, Year <= YEAR_MAX)

# Link & Sauer (2002) effort parameters (muP = p, muB = B)
mu <- read.csv(path_data(FILES$link_sauer)) %>%
  rename(Species = BLsciname) %>%
  dplyr::select(Species, muP, muB)

# ---- 2. generate zero counts -------------------------------------------------------
# Methods: "compiling a species list for each site across the full survey period
# and treating any listed species absent from a particular count as zero".
data_complete <- data %>%
  group_by(SiteCode) %>%
  tidyr::complete(
    Year    = unique(Year),       # survey years at this site
    Species = unique(Species),    # species ever recorded at this site
    fill = list(Count = 0)
  ) %>%
  arrange(SiteCode, Year, Species) %>%
  left_join(mu, by = "Species") %>%
  ungroup()

saveRDS(data_complete, path_inter("qi_data_complete.rds"))

# ---- 3. site-level species change ----------------------------------------------------
# gain / loss / increasing / declining / no significant change
results_all <- analyze_species_trends_cbc(data_complete)
write.csv(results_all, path_inter("qi_species_change.csv"), row.names = FALSE)

cat("\nSite x species series by category:\n")
print(table(results_all$Category))
cat("\nModel used:\n")
print(table(results_all$model_type, useNA = "ifany"))
cat("\nEffort offset applied (TRUE) or not (FALSE):\n")
print(table(results_all$used_offset))

# Diagnostic: series whose model could not be fitted
n_fallback_lost <- sum(is.na(results_all$model_type) &
                         results_all$Category == "No significant change species")
cat("\nSeries labelled 'No significant change' because no model could be fitted:",
    n_fallback_lost, "\n")

# ---- 4. initial relative abundance --------------------------------------------------
initial_ra <- calculate_relative_abundance_first_year(data_complete) %>%
  dplyr::select(SiteCode, Species, relative_abundance)

results_ra <- merge(results_all, initial_ra, by = c("SiteCode", "Species"))

# species absent from the first count have -Inf (log10(0)) and are excluded
results_ra_finite <- results_ra[is.finite(results_ra$relative_abundance), ]
results_ra_finite$pSAD2 <- bin_relative_abundance(results_ra_finite$relative_abundance)

write.csv(results_ra_finite, path_inter("qi_species_change_initial_ra.csv"), row.names = FALSE)

# ---- 5. check: why is relative abundance NA? ---------------------------------------------
# (a site whose largest first-count abundance is 1 gives log10(1) = 0 in the denominator)
diag_input <- data %>%
  filter(Count > 0) %>%
  dplyr::select(SiteCode, Species, Year, Count)
ra_diag <- diagnose_relative_abundance_na(diag_input)
cat("\nReasons for NA initial relative abundance (counts > 0 only):\n")
print(table(ra_diag$na_reason, useNA = "ifany"))
cat("Sites:", length(unique(ra_diag$SiteCode)), "\n")

save_session_info("02_qi_species_change")
