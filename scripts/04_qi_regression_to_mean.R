# =============================================================================
# 04_qi_regression_to_mean.R
# Methods: Question i -> "Testing for regression-to-the-mean"  (Fig. S4)
#
# At the sites with >= 3 counts in 1990-1999 (data/rtm_site.csv), the species-level
# designations, initial-commonness classes and abundance-change models are
# repeated using (0) all years and versions left-censored by (1), (2) and (3)
# years, i.e. with the earliest 1, 2 or 3 counts of every site removed.
#
# Input : output/intermediate/qi_data_complete.rds, qi_species_change.csv (script 02)
#         data/rtm_site.csv
# Output: output/intermediate/rtm_species_change_<scenario>.csv
#         output/tables/rtm_models.csv
#         output/figures/rtm_fig_S4_v1.png, rtm_fig_S4_v2.png
#
# =============================================================================

source("R/config.R")
source("R/functions_species_change.R")
source("R/functions_plotting.R")
library(dplyr)
library(lme4)
library(effects)
library(ggplot2)
library(ggpubr)

require_files(path_inter("qi_data_complete.rds"), path_inter("qi_species_change.csv"),
              path_data(FILES$rtm_sites))

data_complete <- readRDS(path_inter("qi_data_complete.rds"))
results_all   <- read.csv(path_inter("qi_species_change.csv"))
rtm_sites     <- read.csv(path_data(FILES$rtm_sites))$x

# ---- check the site list against the Methods rule ---------------------------------------
# Methods: "sites that met the inclusion criteria for Question i and had at least
# three counts from 1990-1999" (a count = one site-year here)
n_early <- data_complete %>%
  filter(Year >= 1990, Year <= 1999) %>%
  distinct(SiteCode, Year) %>%
  count(SiteCode)
cat("Sites in rtm_site.csv:", length(rtm_sites), "\n")
cat("Sites with >= 3 survey years in 1990-1999 in the data:",
    sum(n_early$n >= 3), "\n")

rtm_data <- data_complete %>% filter(SiteCode %in% rtm_sites)

# ---- run the four scenarios ---------------------------------------------------------------
scenarios <- c(all_years = 0, censor_1yr = 1, censor_2yr = 2, censor_3yr = 3)
scenario_labels <- c(all_years  = "All years",
                     censor_1yr = "First count removed",
                     censor_2yr = "First two counts removed",
                     censor_3yr = "First three counts removed")

res <- list()      # per scenario: $b = finite initial abundance, $c = + fitted rate
fits <- list()     # per scenario: fitted model, effect data frame
current <- rtm_data

for (s in names(scenarios)) {
  k <- scenarios[[s]]
  message("Scenario: ", s)

  if (k == 0) {
    # the change analysis is run separately for each site x species series, so the
    # all-years result for these sites is a subset of the full result
    res_s <- results_all %>% filter(SiteCode %in% rtm_sites)
  } else {
    current <- drop_first_year(current)
    res_s   <- analyze_species_trends_cbc(current)
  }
  write.csv(res_s, path_inter(paste0("rtm_species_change_", s, ".csv")), row.names = FALSE)

  # initial relative abundance from the first remaining count of each site
  ra_s <- calculate_relative_abundance_first_year(current) %>%
    dplyr::select(SiteCode, Species, relative_abundance)

  merged <- merge(res_s, ra_s, by = c("SiteCode", "Species"))
  merged_b <- merged[is.finite(merged$relative_abundance), ]
  merged_b$pSAD2 <- bin_relative_abundance(merged_b$relative_abundance)
  merged_c <- merged_b[!is.na(merged_b$estimate), ]

  res[[s]] <- list(b = merged_b, c = merged_c)

  # rate of change ~ initial relative abundance, species as random factor
  dat_c <- merged_c
  dat_c$Sp2 <- as.factor(dat_c$Species)
  m <- lmer(estimate ~ relative_abundance + (1 | Sp2), data = dat_c)
  fits[[s]] <- list(model = m,
                    fit = as.data.frame(effects::effect(term = "relative_abundance", mod = m)))
}

# ---- model table ----------------------------------------------------------------------------
rtm_table <- do.call(rbind, lapply(names(fits), function(s) {
  m  <- fits[[s]]$model
  cf <- coef(summary(m))["relative_abundance", ]
  an <- car::Anova(m)
  data.frame(scenario = s,
             n_sites = n_distinct(res[[s]]$c$SiteCode),
             n_series = nrow(res[[s]]$c),
             estimate = cf[["Estimate"]], std_error = cf[["Std. Error"]],
             t_value = cf[["t value"]],
             chisq = an["relative_abundance", "Chisq"],
             p_value = an["relative_abundance", "Pr(>Chisq)"])
}))
print(rtm_table)
write.csv(rtm_table, path_table("rtm_models.csv"), row.names = FALSE)

# ---- Fig. S4 -----------------------------------------------------------------------------------
bars <- lapply(names(res), function(s) {
  d <- res[[s]]$b %>% group_by(pSAD2, Category) %>% summarise(sp_no = n(), .groups = "drop")
  plot_category_bar(d) + labs(title = scenario_labels[[s]]) +
    theme(axis.text.x = element_text(size = 7))
})
points_full <- lapply(names(res), function(s) {
  plot_estimate_scatter(res[[s]]$c, fit = fits[[s]]$fit) + ylim(-3, 3)
})
points_zoom <- lapply(names(res), function(s) {
  plot_estimate_scatter(res[[s]]$c, fit = fits[[s]]$fit) + ylim(-0.5, 0.5)
})

fig_s4_v1 <- ggarrange(plotlist = c(bars, points_full), ncol = 4, nrow = 2,
                       common.legend = TRUE, legend = "bottom")
ggsave(path_fig("rtm_fig_S4_v1.png"), fig_s4_v1, width = 15, height = 6, dpi = 500)

fig_s4_v2 <- ggarrange(plotlist = c(bars, points_zoom, points_full), ncol = 4, nrow = 3,
                       common.legend = TRUE, legend = "bottom")
ggsave(path_fig("rtm_fig_S4_v2.png"), fig_s4_v2, width = 15, height = 9, dpi = 500)

save_session_info("04_qi_regression_to_mean")
