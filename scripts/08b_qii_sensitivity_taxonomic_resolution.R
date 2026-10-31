# =============================================================================
# 08b_qii_sensitivity_taxonomic_resolution.R
# Methods: Question ii -> "Sensitivity tests", taxonomic resolution  (Fig. S7)
#
# CBC counts: mean alpha calculated (a) from species-level records only
# (id_all.csv) and (b) after also treating each unidentified taxon group
# (e.g. "gull sp.") as one species (data/alpha_unid/).
#
# Input : data/id_all.csv      yr, SiteCode, mean_alpha, sy (= year + SiteCode)
#         data/alpha_unid/<n>_gambin.csv   gambin output (script 06 layout) for the data
#              WITH unidentified taxa; n = CBC count number (year - 1900), 91..121.
#              Column 2 = count number, 3 = SiteCode, 4 = mean alpha.
# Output: output/tables/sens_taxonomic_*.csv, output/intermediate/sy_taxonomic_test.csv
#         output/figures/fig_S7_taxonomic_resolution.png
# =============================================================================

source("R/config.R")
library(dplyr)
library(ggplot2)

UNID_DIR   <- file.path(DIR_DATA, "alpha_unid")
UNID_INDEX <- 91:121
require_files(path_data(FILES$cbc_alpha), file.path(UNID_DIR, paste0(UNID_INDEX, "_gambin.csv")))

# ---- alpha with unidentified taxa ------------------------------------------------------------
# (the 1991 file has one column fewer than the others; only columns 2-4 are used)
unid_all <- bind_rows(lapply(UNID_INDEX, function(x) {
  d <- read.csv(file.path(UNID_DIR, paste0(x, "_gambin.csv")))
  data.frame(count_no = d[[2]], sitecode = d[[3]], unid_alpha = d[[4]])
}))
unid_all$sy <- paste0(unid_all$count_no + 1900, unid_all$sitecode)

# ---- alpha from species-level records only -----------------------------------------------------
id_all <- read.csv(path_data(FILES$cbc_alpha)) %>%
  dplyr::select(yr, SiteCode, mean_alpha, sy)

all2 <- inner_join(id_all, unid_all, by = "sy") %>%
  filter(mean_alpha <= ALPHA_MAX & unid_alpha <= ALPHA_MAX)

cat("Sites:", n_distinct(all2$SiteCode), "\n")
cat("Counts:", nrow(all2), "\n")

# ---- agreement between the two versions ----------------------------------------------------------
model <- lm(mean_alpha ~ unid_alpha, data = all2)
print(summary(model))
write.csv(data.frame(term = rownames(coef(summary(model))), coef(summary(model))),
          path_table("sens_taxonomic_lm.csv"), row.names = FALSE)

cors <- data.frame(
  method   = c("pearson", "spearman", "kendall"),
  estimate = sapply(c("pearson", "spearman", "kendall"),
                    function(m) cor(all2$mean_alpha, all2$unid_alpha, method = m)))
print(cors)
write.csv(cors, path_table("sens_taxonomic_correlations.csv"), row.names = FALSE)

write.csv(data.frame(sy = all2$sy), path_inter("sy_taxonomic_test.csv"), row.names = FALSE)

# ---- Fig. S7 ---------------------------------------------------------------------------------------
p7 <- ggplot(all2, aes(x = mean_alpha, y = unid_alpha)) +
  geom_point(colour = "#5A5A5A", alpha = 0.3) +
  geom_smooth(method = lm, se = TRUE) +
  theme_classic() +
  labs(x = "Mean SAD alpha without non-species level records",
       y = "Mean SAD alpha with non-species level records")
ggsave(path_fig("fig_S7_taxonomic_resolution.png"), p7, width = 5, height = 5, dpi = 300)

save_session_info("08b_qii_sensitivity_taxonomic_resolution")
