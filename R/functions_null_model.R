# =============================================================================
# functions_null_model.R
# Random-walk null model for abundance change.
# Methods: "Null model for abundance change"
# Used by scripts/05_qi_null_model.R
#
# Same model rule as the empirical analysis: a negative binomial GLM is fitted
# first; a quasi-Poisson GLM is used only if the negative binomial fit fails.
# (Unlike the empirical analysis there is no constant added to the counts and
# no effort offset; a population that reaches zero is classed as "Extinct".)
# =============================================================================

# -----------------------------------------------------------------------------
# Fit abundance ~ year; returns the year coefficient and the model used.
# -----------------------------------------------------------------------------
fit_year_trend <- function(species_data) {
  out <- list(model_type = NA_character_, estimate = NA_real_, std_error = NA_real_,
              z_value = NA_real_, p_value = NA_real_, dispersion = NA_real_)

  # 1. negative binomial (convergence warnings are suppressed; only errors trigger the fall-back)
  fit <- tryCatch(suppressWarnings(MASS::glm.nb(Abundance ~ Year, data = species_data)),
                  error = function(e) NULL)
  type <- "negative_binomial"

  # 2. quasi-Poisson if the negative binomial failed
  if (is.null(fit)) {
    fit  <- tryCatch(suppressWarnings(glm(Abundance ~ Year, family = quasipoisson(link = "log"),
                                          data = species_data)),
                     error = function(e) NULL)
    type <- "quasipoisson"
  }
  if (is.null(fit)) {
    out$model_type <- "all_failed"
    return(out)
  }

  sm <- summary(fit)
  out$model_type <- type
  out$estimate   <- unname(coef(sm)["Year", 1])
  out$std_error  <- unname(coef(sm)["Year", 2])
  out$z_value    <- unname(coef(sm)["Year", 3])
  out$p_value    <- unname(coef(sm)["Year", 4])
  out$dispersion <- if (type == "quasipoisson") sm$dispersion else NA_real_
  out
}

# -----------------------------------------------------------------------------
# Classify one simulated population time series
# -----------------------------------------------------------------------------
classify_species4 <- function(species_data) {
  # a population that reaches zero in any year is classed as extinct
  if (any(species_data$Abundance < 1)) {
    return(list(classification = "Extinct", model_type = NA_character_, estimate = NA_real_,
                std_error = NA_real_, z_value = NA_real_, p_value = NA_real_,
                dispersion = NA_real_))
  }

  res <- fit_year_trend(species_data)

  res$classification <-
    if (is.na(res$p_value)) {
      "Error in Classification"
    } else if (res$p_value < ALPHA_LEVEL) {
      if (res$estimate > 0) "Increased" else "Declined"
    } else {
      "No Significant Change"
    }
  res
}

# -----------------------------------------------------------------------------
# Simulate assemblages and classify every species
#
#   n_species               number of species per assemblage
#   initial_abundance_mean  mean of log10 initial abundance
#   initial_abundance_sd    sd of log10 initial abundance
#   initial_population      optional fixed vector of initial abundances
#   sd_change               sd of the annual (log) rate of change
#   mean_change             mean of the annual (log) rate of change (0 = null)
#   n_years                 number of time points (years) simulated
#   n_simulations           number of simulated assemblages
#
# Relative_Abundance = log10(initial abundance) / log10(largest initial abundance in the
# SAME simulated assemblage), as for the empirical data (relative to the site maximum).
# Returns a list; $classifications is a species-by-simulation data frame.
# -----------------------------------------------------------------------------
simulate_and_classify_populations5 <- function(n_species = 100,
                                               initial_abundance_mean = 1.44,
                                               initial_abundance_sd = 1,
                                               initial_population = NULL,
                                               sd_change = 0.1,
                                               n_years = 28,
                                               n_simulations = 1000,
                                               mean_change = 0) {

  if (!is.null(initial_population) && length(initial_population) != n_species) {
    stop("Length of initial_population must match n_species")
  }

  # years x species x simulations
  abundance_array <- array(0, dim = c(n_years, n_species, n_simulations))
  all_classifications <- vector("list", n_simulations)

  for (sim in seq_len(n_simulations)) {
    # fresh initial abundances for every simulated assemblage
    sim_initial_pops <- if (is.null(initial_population)) {
      ceiling(10^rnorm(n_species, mean = initial_abundance_mean, sd = initial_abundance_sd))
    } else {
      initial_population
    }
    abundance_array[1, , sim] <- sim_initial_pops

    # exponential growth with a new random rate for each year and species;
    # values < 1 become 0, others are rounded to integers
    for (year in 2:n_years) {
      rates <- rnorm(n_species, mean = mean_change, sd = sd_change)
      raw_abundance <- abundance_array[year - 1, , sim] * exp(rates)
      abundance_array[year, , sim] <- ifelse(raw_abundance < 1, 0, round(raw_abundance))
    }

    # relative abundance within this assemblage
    max_initial <- max(sim_initial_pops)
    rel_abund <- if (max_initial > 1) log10(sim_initial_pops) / log10(max_initial) else
      rep(0, n_species)

    sim_results <- lapply(seq_len(n_species), function(sp) {
      sp_data <- data.frame(Year = 1:n_years, Abundance = abundance_array[, sp, sim])
      r <- classify_species4(sp_data)
      data.frame(
        Species            = paste0("Species_", sp),
        Initial_Abundance  = abundance_array[1, sp, sim],
        Final_Abundance    = abundance_array[n_years, sp, sim],
        Relative_Abundance = rel_abund[sp],
        Classification     = r$classification,
        Model_Type         = r$model_type,
        Estimate           = r$estimate,
        Std_Error          = r$std_error,
        Z_value            = r$z_value,
        P_value            = r$p_value,
        Dispersion         = r$dispersion,
        Trial              = sim
      )
    })
    all_classifications[[sim]] <- do.call(rbind, sim_results)

    if (sim %% 100 == 0) cat("Completed simulation", sim, "of", n_simulations, "\n")
  }

  all_classifications_df <- do.call(rbind, all_classifications)

  list(
    mean_abundance = apply(abundance_array, c(1, 2), mean),
    lower_ci = apply(abundance_array, c(1, 2), function(x) quantile(x, 0.025)),
    upper_ci = apply(abundance_array, c(1, 2), function(x) quantile(x, 0.975)),
    classifications = all_classifications_df,
    classification_summary = table(all_classifications_df$Classification),
    model_summary = table(all_classifications_df$Model_Type, useNA = "ifany"),
    parameters = list(n_species = n_species, n_years = n_years, n_simulations = n_simulations,
                      sd_change = sd_change, mean_change = mean_change,
                      initial_abundance_mean = initial_abundance_mean,
                      initial_abundance_sd = initial_abundance_sd)
  )
}
