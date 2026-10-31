# =============================================================================
# 06b_qii_extract_alpha.R
# Methods: Question ii -> "Estimating SAD shape using the gambin model"
#   "...to obtain a mean alpha parameter ... and associated p-value for each count.
#    Counts with a mean p-value below 0.05 or a mean alpha parameter exceeding 15
#    were excluded."
#
# Collects the yearly gambin output of script 06 into one table and applies the
# two exclusion rules.
#
# Input : output/intermediate/gambin/<year>_gambin.csv   Mean.Stan.Alpha per site
#         output/intermediate/gambin/<year>_x2.csv       X2 p-values of the subsamples
# Output: output/intermediate/alpha_all_counts.csv   every count, with mean_p
#         output/intermediate/alpha_filtered.csv     counts kept (-> script 07)
#
# ASSUMPTION (check): gmod$X2 holds the goodness-of-fit p-value of each subsample,
# one column per site, so the mean p-value of a count is the column mean. The
# gambin documentation describes this element only as "X2 P-values for each site".
# The script stops if the values are not between 0 and 1.
# =============================================================================

source("R/config.R")
library(dplyr)

OUT_GAMBIN <- file.path(DIR_INTER, "gambin")
files <- list.files(OUT_GAMBIN, pattern = "_gambin\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No gambin output found; run script 06 first.")

read_year <- function(f) {
  yr <- as.integer(sub("_gambin\\.csv$", "", basename(f)))
  g  <- read.csv(f)
  x2 <- read.csv(file.path(OUT_GAMBIN, paste0(yr, "_x2.csv")), row.names = 1)

  if (ncol(x2) == nrow(g)) {            # subsamples x sites
    mean_p <- colMeans(x2)
  } else if (ncol(x2) == 1 && nrow(x2) == nrow(g)) {   # one value per site
    mean_p <- x2[[1]]
  } else {
    stop("Year ", yr, ": the X2 file does not match the number of sites.")
  }
  if (any(mean_p < 0 | mean_p > 1, na.rm = TRUE)) {
    stop("Year ", yr, ": X2 values are not between 0 and 1, so they may be statistics, not p-values.")
  }
  data.frame(yr = yr, SiteCode = g$SiteCode, mean_alpha = g$Mean.Stan.Alpha,
             mean_p = unname(mean_p), stringsAsFactors = FALSE)
}

alpha_all <- bind_rows(lapply(files, read_year))
write.csv(alpha_all, path_inter("alpha_all_counts.csv"), row.names = FALSE)

# exclusion rules
alpha_filtered <- alpha_all %>%
  filter(mean_p >= ALPHA_P_MIN, mean_alpha <= ALPHA_MAX) %>%
  mutate(logmean = log10(mean_alpha))

cat("Counts fitted:", nrow(alpha_all), "\n")
cat("Excluded (mean p < ", ALPHA_P_MIN, "): ", sum(alpha_all$mean_p < ALPHA_P_MIN), "\n", sep = "")
cat("Excluded (mean alpha > ", ALPHA_MAX, "): ", sum(alpha_all$mean_alpha > ALPHA_MAX), "\n", sep = "")
cat("Counts kept:", nrow(alpha_filtered), " at ", n_distinct(alpha_filtered$SiteCode), " sites\n", sep = "")
write.csv(alpha_filtered, path_inter("alpha_filtered.csv"), row.names = FALSE)

save_session_info("06b_qii_extract_alpha")
