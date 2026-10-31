# =============================================================================
# functions_species_change.R
# Functions for Question i (site-level species change and initial abundance).
# Sourced by scripts/02, 03 and 04.
# =============================================================================

require(dplyr)

# -----------------------------------------------------------------------------
# Initial relative abundance classes (10% bins)
# Methods: "grouped all species observed at each site into 10% intervals of
# initial relative abundance from 0% to 100%"
# -----------------------------------------------------------------------------
bin_relative_abundance <- function(x) {
  dplyr::case_when(
    x < 0.1               ~ "p0_10",
    x >= 0.1 & x < 0.2    ~ "p10_20",
    x >= 0.2 & x < 0.3    ~ "p20_30",
    x >= 0.3 & x < 0.4    ~ "p30_40",
    x >= 0.4 & x < 0.5    ~ "p40_50",
    x >= 0.5 & x < 0.6    ~ "p50_60",
    x >= 0.6 & x < 0.7    ~ "p60_70",
    x >= 0.7 & x < 0.8    ~ "p70_80",
    x >= 0.8 & x < 0.9    ~ "p80_90",
    x >= 0.9              ~ "p90ab"
  )
}

PSAD_LABELS <- c("p0_10" = "<10",   "p10_20" = "10-20", "p20_30" = "20-30",
                 "p30_40" = "30-40", "p40_50" = "40-50", "p50_60" = "50-60",
                 "p60_70" = "60-70", "p70_80" = "70-80", "p80_90" = "80-90",
                 "p90ab"  = "90-100")

# -----------------------------------------------------------------------------
# Initial relative abundance of each species at each site
# Methods: "percentage of the log10-transformed abundance of the focal species
# relative to the log10-transformed abundance of the most abundant species in
# that assemblage", using the earliest count at each site.
# Returned as a proportion (0-1). Species with a zero count in the first count
# give -Inf and are removed later with is.finite().
# -----------------------------------------------------------------------------
calculate_relative_abundance_first_year <- function(data,
                                                    site_col = "SiteCode",
                                                    year_col = "Year",
                                                    count_col = "Count",
                                                    species_col = "Species") {
  require(dplyr)

  # Step 1: earliest year of each site
  data_with_first_year <- data %>%
    dplyr::group_by(.data[[site_col]]) %>%
    dplyr::mutate(site_min_year = min(.data[[year_col]])) %>%
    dplyr::ungroup()

  # Step 2: keep records from that year only
  first_year_data <- data_with_first_year %>%
    dplyr::filter(.data[[year_col]] == site_min_year)

  # Step 3: largest count at the site in that year
  first_year_data <- first_year_data %>%
    dplyr::group_by(.data[[site_col]]) %>%
    dplyr::mutate(site_max_count_first_year = max(.data[[count_col]])) %>%
    dplyr::ungroup()

  # Step 4: relative abundance
  result_data <- first_year_data %>%
    dplyr::mutate(
      relative_abundance = log10(.data[[count_col]]) / log10(site_max_count_first_year)
    ) %>%
    dplyr::select(-site_min_year, -site_max_count_first_year)

  return(result_data)
}

# Same calculation, but also reports WHY a relative abundance is NA
# (diagnostic only; used in scripts/02).
diagnose_relative_abundance_na <- function(data,
                                           site_col = "SiteCode",
                                           year_col = "Year",
                                           count_col = "Count") {
  require(dplyr)

  data %>%
    dplyr::group_by(.data[[site_col]]) %>%
    dplyr::mutate(site_min_year = min(.data[[year_col]], na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data[[year_col]] == site_min_year) %>%
    dplyr::group_by(.data[[site_col]]) %>%
    dplyr::mutate(site_max_count_first_year = max(.data[[count_col]], na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      relative_abundance = log10(.data[[count_col]]) / log10(site_max_count_first_year),
      na_reason = dplyr::case_when(
        is.na(.data[[count_col]])                 ~ "Missing count value",
        .data[[count_col]] <= 0                   ~ "Count is zero or negative",
        is.na(site_max_count_first_year)          ~ "Missing maximum count",
        site_max_count_first_year <= 0            ~ "Maximum count is zero or negative",
        site_max_count_first_year == 1            ~ "Maximum count is 1 (log10(1) = 0, division by zero)",
        is.na(relative_abundance) & !is.na(.data[[count_col]]) &
          !is.na(site_max_count_first_year)       ~ "Other calculation error",
        TRUE                                      ~ NA_character_
      )
    ) %>%
    dplyr::select(-site_min_year, -site_max_count_first_year)
}

# -----------------------------------------------------------------------------
# Helper: remove the earliest year of every site (regression-to-the-mean test)
# Methods: "versions left-censored by one, two, and three years"
# -----------------------------------------------------------------------------
drop_first_year <- function(data, site_col = "SiteCode", year_col = "Year") {
  data %>%
    dplyr::group_by(.data[[site_col]]) %>%
    dplyr::filter(.data[[year_col]] > min(.data[[year_col]])) %>%
    dplyr::ungroup()
}

# -----------------------------------------------------------------------------
# Site-level species change
# Methods: "Defining site-level species changes"
#
# For every site x species series:
#   Species gain : present in the first decade? no  & present in the last decade? yes
#   Species loss : present in the first decade? yes & present in the last decade? no
#   otherwise    : slope of abundance ~ year from a negative binomial GLM
#                  (quasi-Poisson if the negative binomial fails), with the
#                  Link & Sauer effort offset where Hours, muP and muB exist;
#                  p < 0.05 -> Increasing / Declining Species, else No
#                  significant change species.
#
# All species are classified, however rarely they were recorded.
# -----------------------------------------------------------------------------
analyze_species_trends_cbc <- function(data_complete,
                                       first_period = PERIOD_FIRST,
                                       last_period  = PERIOD_LAST) {
  results_list <- list()

  site_species <- unique(data_complete[c("SiteCode", "Species")])

  for (i in seq_len(nrow(site_species))) {
    site    <- site_species$SiteCode[i]
    species <- site_species$Species[i]

    # defaults
    category   <- NA_character_
    model_type <- NA_character_
    estimate   <- NA_real_
    std_error  <- NA_real_
    t_value    <- NA_real_
    p_value    <- NA_real_
    used_offset <- FALSE

    species_data <- data_complete %>%
      dplyr::filter(SiteCode == site, Species == species)

    # presence (non-zero counts) in the first and last decade
    present_first <- any(species_data$Count[species_data$Year >= first_period[1] &
                                              species_data$Year <= first_period[2]] > 0)
    present_last  <- any(species_data$Count[species_data$Year >= last_period[1] &
                                              species_data$Year <= last_period[2]] > 0)

    if (!present_first && present_last) {
      category <- "Species gain"
    } else if (present_first && !present_last) {
      category <- "Species loss"
    } else {
      if (sum(species_data$Count > 0) > 0) {
        # small constant, as in the original analysis
        species_data$Count_adj <- species_data$Count + 0.0001

        # effort offset only where Hours, muP and muB are all available
        has_offset <- all(c("muP", "muB", "Hours") %in% names(species_data)) &&
          !any(is.na(species_data$muP)) &&
          !any(is.na(species_data$muB)) &&
          !any(is.na(species_data$Hours))

        if (has_offset) {
          # Link & Sauer effort term: B * ((hours / mean hours)^p - 1) / p
          offset_term <- with(species_data,
                              muB * ((Hours / mean(Hours, na.rm = TRUE))^muP - 1) / muP)

          tryCatch({
            nb_model <- MASS::glm.nb(Count_adj ~ Year + offset(offset_term), data = species_data)
            model_summary <- summary(nb_model)
            model_type <- "negative_binomial_with_offset"
            estimate   <- coef(model_summary)[2]
            std_error  <- model_summary$coefficients[2, 2]
            t_value    <- model_summary$coefficients[2, 3]
            p_value    <- model_summary$coefficients[2, 4]
            used_offset <- TRUE
          }, error = function(e) {
            tryCatch({
              quasi_model <- glm(Count_adj ~ Year + offset(offset_term),
                                 family = quasipoisson(), data = species_data)
              model_summary <- summary(quasi_model)
              model_type <- "quasipoisson_with_offset"
              estimate   <- coef(model_summary)[2]
              std_error  <- model_summary$coefficients[2, 2]
              t_value    <- model_summary$coefficients[2, 3]
              p_value    <- model_summary$coefficients[2, 4]
              used_offset <- TRUE
            }, error = function(e) {
              model_type <- "model_failed"
            })
          })
        } else {
          tryCatch({
            nb_model <- MASS::glm.nb(Count_adj ~ Year, data = species_data)
            model_summary <- summary(nb_model)
            model_type <- "negative_binomial"
            estimate   <- coef(model_summary)[2]
            std_error  <- model_summary$coefficients[2, 2]
            t_value    <- model_summary$coefficients[2, 3]
            p_value    <- model_summary$coefficients[2, 4]
            used_offset <- FALSE
          }, error = function(e) {
            tryCatch({
              quasi_model <- glm(Count_adj ~ Year, family = quasipoisson(), data = species_data)
              model_summary <- summary(quasi_model)
              model_type <- "quasipoisson"
              estimate   <- coef(model_summary)[2]
              std_error  <- model_summary$coefficients[2, 2]
              t_value    <- model_summary$coefficients[2, 3]
              p_value    <- model_summary$coefficients[2, 4]
              used_offset <- FALSE
            }, error = function(e) {
              model_type <- "model_failed"
            })
          })
        }

        # category from statistical significance
        if (is.na(p_value) || model_type == "model_failed") {
          category <- "No significant change species"
        } else if (p_value < ALPHA_LEVEL) {
          if (estimate > 0) {
            category <- "Increasing Species"
          } else {
            category <- "Declining Species"
          }
        } else {
          category <- "No significant change species"
        }
      } else {
        category   <- "No significant change species"
        model_type <- "all_zeros"
      }
    }

    results_list[[i]] <- data.frame(
      SiteCode        = site,
      Species         = species,
      Category        = category,
      model_type      = model_type,
      estimate        = estimate,
      std_error       = std_error,
      t_value         = t_value,
      p_value         = p_value,
      used_offset     = used_offset,
      total_records   = nrow(species_data),
      nonzero_records = sum(species_data$Count > 0),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, results_list)
}
