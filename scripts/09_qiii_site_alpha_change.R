# =============================================================================
# 09_qiii_site_alpha_change.R
# Methods: Question iii -> "Quantifying assemblage and SAD shape changes"
#   "SAD shape change at each site was quantified as the estimated slope of a
#    linear model with log10-transformed alpha as the response variable and year
#    as the explanatory variable."
#
# Input : data/allPn2c_clean.csv   (logmean, yr, SiteCode, region)
# Output: output/intermediate/site_alpha_change.csv
#         columns SiteCode, estimate, std.error, statistic, p.value, change, region
#         (estimate = slope of log10 alpha per year)
# =============================================================================

source("R/config.R")
library(dplyr)
library(broom)

require_files(path_data(FILES$alpha_clean))
alpha_clean <- read.csv(path_data(FILES$alpha_clean))

MIN_COUNTS_FOR_SLOPE <- 3     # a slope with a standard error needs at least 3 counts

site_change <- alpha_clean %>%
  group_by(SiteCode) %>%
  filter(n() >= MIN_COUNTS_FOR_SLOPE) %>%
  group_modify(~ broom::tidy(lm(logmean ~ yr, data = .x)) %>%
                 filter(term == "yr") %>%
                 dplyr::select(-term)) %>%
  ungroup()

# direction of change: significant increase / decrease, otherwise fluctuating
site_change <- site_change %>%
  mutate(change = case_when(p.value < ALPHA_LEVEL & estimate > 0 ~ "increase",
                            p.value < ALPHA_LEVEL & estimate < 0 ~ "decrease",
                            TRUE                                  ~ "fluctuate"))

# region of each site
site_region <- alpha_clean %>% distinct(SiteCode, region)
if (anyDuplicated(site_region$SiteCode)) warning("Some sites have more than one region.")
site_change <- site_change %>% left_join(site_region, by = "SiteCode")

print(table(site_change$change))
cat("Sites:", nrow(site_change), "\n")
write.csv(site_change, path_inter("site_alpha_change.csv"), row.names = FALSE)

save_session_info("09_qiii_site_alpha_change")
