# =============================================================================
# 06_qii_gambin_alpha.R
# Methods: Question ii -> "Estimating SAD shape using the gambin model"
#
# For every year, counts with >= 10 species and >= 1,000 individuals are
# subsampled (1,000 individuals), the gambin model is fitted, and this is
# repeated 1,000 times per count (gambin::mult_abundances).
#
# Input : data/w90.csv  (columns Year, SiteCode, Species, Count)
# Output: output/intermediate/gambin/<year>_all_alpha.csv   alpha of every subsample
#         output/intermediate/gambin/<year>_x2.csv          chi-squared statistics
#         output/intermediate/gambin/<year>_gambin.csv      mean alpha, median, SE per site
#
# Requires the gambin package:  install.packages("gambin")
# =============================================================================

source("R/config.R")
library(gambin)
library(dplyr)
library(tidyr)

require_files(path_data(FILES$qii_counts))
set.seed(SEED)

OUT_GAMBIN <- file.path(DIR_INTER, "gambin")
dir.create(OUT_GAMBIN, showWarnings = FALSE, recursive = TRUE)
OVERWRITE <- FALSE      # TRUE = refit years that already have output

sdata <- read.csv(path_data(FILES$qii_counts))

fit_gambin_year <- function(sdata, yr) {

  # counts of the selected year
  sdata2 <- dplyr::filter(sdata, Year == yr)
  ssites <- sdata2 %>% dplyr::select(Year, SiteCode, Species, Count)

  # add up the counts of each species at each site
  df2 <- ssites %>%
    group_by(SiteCode, Species) %>%
    summarise(Sum = sum(Count), .groups = "drop") %>%
    as.data.frame()

  # sites with >= 1,000 individuals and >= 10 species
  slist <- df2 %>%
    group_by(SiteCode) %>%
    summarise(n = sum(Sum), spp_no = length(Species), .groups = "drop")
  slist6 <- slist %>%
    filter(n >= GAMBIN_MIN_INDIVIDUALS, spp_no >= GAMBIN_MIN_SPECIES)

  if (nrow(slist6) == 0) {
    message("  no site qualifies in ", yr)
    return(invisible(NULL))
  }

  # species x site matrix of abundances
  df4 <- dplyr::filter(df2, SiteCode %in% slist6$SiteCode)
  df5 <- tidyr::spread(df4, SiteCode, Sum)
  df5 <- dplyr::select(df5, -1)          # drop the Species column
  df5[is.na(df5)] <- 0

  # gambin on every site (columns)
  gmod <- mult_abundances(df5, N = GAMBIN_N, subsample = GAMBIN_SUBSAMPLE, r = GAMBIN_R)

  # summaries of alpha across subsamples
  x     <- gmod$All.Alpha
  SE50  <- apply(x, 2, function(v) sd(v) / sqrt(length(v)))
  med   <- apply(x, 2, function(v) median(v))

  # site names are taken from the columns of the matrix that was fitted
  yrdata <- data.frame(yr, SiteCode = colnames(df5),
                       Mean.Stan.Alpha = gmod$Mean.Stan.Alpha, med, SE50)

  write.csv(x,         file.path(OUT_GAMBIN, paste0(yr, "_all_alpha.csv")))
  write.csv(gmod$X2,   file.path(OUT_GAMBIN, paste0(yr, "_x2.csv")))
  write.csv(yrdata,    file.path(OUT_GAMBIN, paste0(yr, "_gambin.csv")))
  invisible(yrdata)
}

for (yr in YEAR_MIN:YEAR_MAX) {
  done <- file.path(OUT_GAMBIN, paste0(yr, "_gambin.csv"))
  if (file.exists(done) && !OVERWRITE) {
    message("Skipping ", yr, " (already fitted)")
    next
  }
  message("Fitting gambin for ", yr, "  [", format(Sys.time()), "]")
  fit_gambin_year(sdata, yr)
}

save_session_info("06_qii_gambin_alpha")
