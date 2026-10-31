# =============================================================================
# 05_qi_null_model.R
# Methods: Question i -> "Null model for abundance change"  (Fig. S1)
#
# 1,000 simulated assemblages of 100 species; log10 initial abundance drawn from
# N(1.44, 1.00); annual log rates of change drawn from N(0, sd) with
# sd = 0.1, 0.2 and 0.3; 28 time points. Each population is classified and its
# rate of change estimated, then rate of change is regressed on initial
# relative abundance.
#
# Output: output/intermediate/null_classifications_sd<sd>.rds
#         output/tables/null_models.csv
#         output/figures/null_model_fig_S1.png
# =============================================================================

source("R/config.R")
source("R/functions_species_change.R")    # bin_relative_abundance(), PSAD_LABELS
source("R/functions_null_model.R")
source("R/functions_plotting.R")
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)

set.seed(SEED)

# ---- 1. simulate ------------------------------------------------------------------------------
null_runs <- list()
for (s in NULL_SD_CHANGE) {
  nm <- paste0("sd", s)
  message("Simulating sd_change = ", s)
  null_runs[[nm]] <- simulate_and_classify_populations5(
    n_species              = NULL_N_SPECIES,
    initial_abundance_mean = NULL_INIT_MEAN,
    initial_abundance_sd   = NULL_INIT_SD,
    sd_change              = s,
    n_years                = NULL_N_YEARS,
    n_simulations          = NULL_N_SIM
  )
  saveRDS(null_runs[[nm]]$classifications,
          path_inter(paste0("null_classifications_", nm, ".rds")))
}

# models used (negative binomial / quasi-Poisson fall-back) and classifications
model_use <- do.call(rbind, lapply(names(null_runs), function(nm) {
  d <- null_runs[[nm]]$classifications
  data.frame(scenario = nm, as.data.frame(table(model_type = d$Model_Type, useNA = "ifany")))
}))
print(model_use)
write.csv(model_use, path_table("null_model_types.csv"), row.names = FALSE)
class_use <- do.call(rbind, lapply(names(null_runs), function(nm) {
  data.frame(scenario = nm, as.data.frame(null_runs[[nm]]$classification_summary))
}))
write.csv(class_use, path_table("null_model_classifications.csv"), row.names = FALSE)

# ---- 2. rate of change ~ initial relative abundance ------------------------------------------------
null_models <- list()
for (nm in names(null_runs)) {
  d <- null_runs[[nm]]$classifications
  m <- lm(Estimate ~ Relative_Abundance, data = d)
  print(summary(m))
  print(car::Anova(m))
  cf <- coef(summary(m))["Relative_Abundance", ]
  null_models[[nm]] <- data.frame(scenario = nm,
                                  n = nobs(m),
                                  estimate = cf[["Estimate"]],
                                  std_error = cf[["Std. Error"]],
                                  t_value = cf[["t value"]],
                                  p_value = cf[["Pr(>|t|)"]])
}
null_models <- do.call(rbind, null_models)
write.csv(null_models, path_table("null_models.csv"), row.names = FALSE)

# ---- 3. Fig. S1 ------------------------------------------------------------------------------------
axis_lab <- "Initial Relative Abundance (%)"
panels <- list()
for (nm in names(null_runs)) {
  d <- null_runs[[nm]]$classifications
  panels[[paste0(nm, "_bar")]] <- plot_bar(d) +
    theme(axis.text.x = element_text(angle = 0, size = 7, hjust = 0.5)) +
    guides(x = guide_axis(angle = 45)) + My_Theme + labs(x = axis_lab)
  panels[[paste0(nm, "_pt")]] <- plot_point(d) + ylim(-0.5, 0.5) + My_Theme +
    labs(x = axis_lab) +
    scale_x_continuous(labels = percent_format(scale = 100, suffix = ""))
}

fig_s1 <- wrap_plots(panels, ncol = 2, byrow = TRUE) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(path_fig("null_model_fig_S1.png"), fig_s1, width = 8, height = 10, dpi = 300)

save_session_info("05_qi_null_model")
