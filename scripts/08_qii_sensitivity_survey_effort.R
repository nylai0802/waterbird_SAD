# =============================================================================
# 08_qii_sensitivity_survey_effort.R
# Methods: Question ii -> "Sensitivity tests", survey effort  (Fig. S6)
#
# Christmas Bird Count (CBC) counts only. Tests whether alpha depends on total
# survey hours and whether including survey hours changes the temporal trend.
#   (a) lm: mean alpha ~ total survey hours
#   (b) lmer: log10 alpha ~ year + survey hours + (1 | site); year-by-year estimates
#   (c) lmer: log10 alpha ~ year + (1 | site);                 year-by-year estimates
#
# Input : data/id_all.csv
#         data/20221108_Nga_Yee_Lai_CBC_Effort_Many_Types_Report.csv
# Output: output/tables/sens_effort_*.csv, output/intermediate/cwc_surv_eff_adj.csv
#         output/figures/fig_S6_survey_effort.png
# =============================================================================

source("R/config.R")
library(dplyr)
library(lme4)
library(effects)
library(ggplot2)
library(ggpubr)

require_files(path_data(FILES$cbc_alpha), path_data(FILES$cbc_effort))

alpha <- read.csv(path_data(FILES$cbc_alpha))
if (!"logmean" %in% names(alpha)) alpha$logmean <- log10(alpha$mean_alpha)

# ---- 1. survey effort per CBC count -------------------------------------------------------------
ef <- read.csv(path_data(FILES$cbc_effort))

# count_yr 1 = 1901 ... count_yr 123 = 2023
ef$yr_ad <- ef$count_yr + 1900
ef <- ef %>% filter(yr_ad >= YEAR_MIN)
ef$sy <- paste0(ef$yr_ad, ef$abbrev)          # site-year id, same format as alpha$sy

eff <- ef %>% filter(sy %in% unique(alpha$sy))
cat("Effort records by method of travel:\n"); print(table(eff$description))

car_hr  <- eff %>% filter(description == "car")  %>% dplyr::select(sy, car  = hours)
foot_hr <- eff %>% filter(description == "foot") %>% dplyr::select(sy, foot = hours)
all_hr  <- eff %>% group_by(sy) %>% summarise(sum_hr = sum(hours), .groups = "drop")

if (anyDuplicated(car_hr$sy) || anyDuplicated(foot_hr$sy)) {
  warning("More than one 'car' or 'foot' effort row for some site-years: the joins below ",
          "will duplicate those counts.")
}

all_hr <- all_hr %>% left_join(car_hr, by = "sy") %>% left_join(foot_hr, by = "sy")
alpha  <- alpha %>% left_join(all_hr, by = "sy")

# CBC counts are those with a recorded car effort (as in the original analysis);
# counts with alpha > 15 are removed
cwc <- alpha %>% filter(!is.na(car)) %>% filter(mean_alpha <= ALPHA_MAX)

cat("\nCBC counts with total survey hours:", sum(!is.na(cwc$sum_hr)), "of", nrow(cwc), "\n")
cat("Sites:", n_distinct(cwc$SiteCode), "\n")
cat("Counts:", nrow(cwc), "\n")

cwc$yrf <- as.factor(cwc$yr)
cwc$SiteCode <- as.factor(cwc$SiteCode)
write.csv(cwc, path_inter("cwc_surv_eff_adj.csv"), row.names = FALSE)

# ---- 2. (a) alpha ~ survey hours ------------------------------------------------------------------
m_a <- lm(mean_alpha ~ sum_hr, data = cwc)
print(summary(m_a))
m_a_log <- lm(logmean ~ sum_hr, data = cwc)      # the plotted quantity (log10 alpha)
print(summary(m_a_log))
write.csv(rbind(cbind(model = "mean_alpha ~ sum_hr", term = rownames(coef(summary(m_a))),
                      coef(summary(m_a))),
                cbind(model = "logmean ~ sum_hr",    term = rownames(coef(summary(m_a_log))),
                      coef(summary(m_a_log)))),
          path_table("sens_effort_alpha_vs_hours.csv"), row.names = FALSE)

# ---- 3. (b), (c) trend with and without survey hours -----------------------------------------------
m_with    <- lmer(logmean ~ yr + (1 | SiteCode) + sum_hr, data = cwc, REML = TRUE)
m_without <- lmer(logmean ~ yr + (1 | SiteCode),          data = cwc, REML = TRUE)
print(summary(m_with));    print(car::Anova(m_with))
print(summary(m_without)); print(car::Anova(m_without))

trend_tab <- rbind(
  data.frame(model = "with survey hours",    term = rownames(coef(summary(m_with))),    coef(summary(m_with))),
  data.frame(model = "without survey hours", term = rownames(coef(summary(m_without))), coef(summary(m_without))))
write.csv(trend_tab, path_table("sens_effort_trend_models.csv"), row.names = FALSE)

# year-by-year estimates for the figure
m_with_yrf    <- lmer(logmean ~ yrf + (1 | SiteCode) + sum_hr, data = cwc, REML = TRUE)
x_with <- as.data.frame(effects::effect(term = "yrf", mod = m_with_yrf))
x_with$yr <- as.numeric(as.character(x_with$yrf))

m_without_yrf <- lmer(logmean ~ yrf + (1 | SiteCode), data = cwc, REML = TRUE)
x_without <- as.data.frame(effects::effect(term = "yrf", mod = m_without_yrf))
x_without$yr <- as.numeric(as.character(x_without$yrf))

# ---- 4. Fig. S6 ---------------------------------------------------------------------------------------
ribbon <- function(x) {
  ggplot(data = x, aes(x = yr, y = fit)) +
    geom_line(linewidth = 0.8) +
    geom_ribbon(aes(x = yr, ymin = lower, ymax = upper), fill = "#4c73b5", alpha = 0.5) +
    labs(x = "Year", y = "Estimate of log10 mean \nalpha value ± 95%CI") +
    scale_x_continuous(breaks = c(1990, 2000, 2010, 2018)) +
    theme_classic() +
    theme(legend.position = "none", text = element_text(size = 12))
}

p6a <- ggplot(cwc, aes(x = sum_hr, y = logmean)) +
  geom_point(size = 1, shape = 19, colour = "#5A5A5A", alpha = 0.2) +
  geom_smooth(method = lm, se = TRUE) +
  xlim(0, max(cwc$sum_hr)) +
  labs(x = "Total number of survey hours of the count",
       y = "log10 mean alpha value of the count") +
  theme_classic()

fig_s6 <- ggarrange(p6a,
                    ggarrange(ribbon(x_with), ribbon(x_without), ncol = 2, labels = c("(b)", "(c)")),
                    labels = "(a)", nrow = 2)
ggsave(path_fig("fig_S6_survey_effort.png"), fig_s6, width = 9, height = 9, dpi = 300)

save_session_info("08_qii_sensitivity_survey_effort")
