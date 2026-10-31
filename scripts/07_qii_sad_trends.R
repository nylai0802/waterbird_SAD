# =============================================================================
# 07_qii_sad_trends.R
# Methods: Question ii -> "Global and regional SAD shape trends"
#          (Fig. 3, Extended Data / Fig. S5)
#
# Linear mixed models of log10(mean alpha) ~ year with site as a random factor,
# for the global dataset and for each region.
#
# Input : output/intermediate/alpha_filtered.csv (script 06b), data/sitelist_country2.csv,
#         data/allPn2c_clean.csv
# Output: output/tables/sad_trend_models_all.csv    global + regional model estimates
#         output/tables/sad_trend_sample_sizes.csv
#         output/tables/site_year_list.csv          site-years used (yr_SiteCode)
#         output/intermediate/allPn2c.csv, duplicate_site_years.csv
#         output/intermediate/sad_yearly_effects_<region>.csv
#         output/figures/fig3_*.png, fig_S5_*.png
# =============================================================================

source("R/config.R")
source("R/functions_plotting.R")
library(dplyr)
library(lme4)
library(effects)
library(broom)          # tidy() method for the car::Anova() table
library(broom.mixed)    # tidy() method for lmer models
library(ggplot2)
library(ggpubr)
library(scales)
library(stringr)

require_files(path_inter("alpha_filtered.csv"), path_data(FILES$site_country),
              path_data(FILES$alpha_clean))

# ---- 1. combine alpha estimates with site country / region ------------------------------------
# alpha estimates after excluding counts with mean p < 0.05 or mean alpha > 15 (script 06b)
alpha_raw <- read.csv(path_inter("alpha_filtered.csv"))
country   <- read.csv(path_data(FILES$site_country))
country   <- country[, 2:8]            # keep columns 2-8 of the site list (as in the original)

alpha_c <- alpha_raw %>% left_join(country, by = "SiteCode")
write.csv(alpha_c, path_inter("allPn2c.csv"), row.names = FALSE)

# duplicated site-years
alpha_c$sy <- paste(alpha_c$yr, "-", alpha_c$SiteCode)
dups <- as.data.frame(table(alpha_c$sy)) %>% filter(Freq > 1)
write.csv(dups, path_inter("duplicate_site_years.csv"), row.names = FALSE)
cat("Duplicated site-years in allPn2c:", nrow(dups), "\n")

# Duplicated site-years and Russian sites were removed before this script (outside R),
# giving data/allPn2c_clean.csv.
alpha_clean <- read.csv(path_data(FILES$alpha_clean))
alpha_clean$SiteCode <- as.factor(alpha_clean$SiteCode)
stopifnot(all(c("logmean", "yr", "SiteCode", "region") %in% names(alpha_clean)))
cat("Sites in the cleaned alpha table:", length(unique(alpha_clean$SiteCode)), "\n")

# site-year identifiers of the counts analysed
alpha_clean$scyr <- paste0(alpha_clean$yr, "_", alpha_clean$SiteCode)
write.csv(data.frame(scyr = alpha_clean$scyr), path_table("site_year_list.csv"),
          row.names = FALSE)

# sample sizes
n_global <- alpha_clean %>% summarise(region = "Global", n_sites = n_distinct(SiteCode),
                                      n_counts = n())
n_region <- alpha_clean %>% group_by(region) %>%
  summarise(n_sites = n_distinct(SiteCode), n_counts = n(), .groups = "drop")
write.csv(rbind(n_global, n_region), path_table("sad_trend_sample_sizes.csv"), row.names = FALSE)

# ---- 2. global model ---------------------------------------------------------------------------
model1 <- lmer(logmean ~ yr + (1 | SiteCode), data = alpha_clean, REML = TRUE)
print(summary(model1))
print(car::Anova(model1))
cat("\nGlobal change in log10 alpha per 100 years:", fixef(model1)[["yr"]] * 100, "\n")

# ---- 3. regional models ------------------------------------------------------------------------
regions <- unique(alpha_clean$region)
all_reg_glmm <- NULL
for (reg in regions) {
  ea  <- dplyr::filter(alpha_clean, region == reg & yr <= YEAR_MAX)
  eam <- lmer(logmean ~ yr + (1 | SiteCode), data = ea, REML = TRUE)
  reg_glmm <- data.frame(reg, broom.mixed::tidy(eam), broom.mixed::tidy(car::Anova(eam)))
  all_reg_glmm <- rbind(all_reg_glmm, reg_glmm)
}

# add the global model
reg <- "Global"
reg_glmm <- data.frame(reg, broom.mixed::tidy(model1), broom.mixed::tidy(car::Anova(model1)))
all_reg_glmm <- rbind(all_reg_glmm, reg_glmm)
write.csv(all_reg_glmm, path_table("sad_trend_models_all.csv"), row.names = FALSE)

reg_change <- all_reg_glmm %>% filter(term %in% "yr")
reg_change$u95 <- reg_change$estimate + reg_change$std.error * 1.96
reg_change$l95 <- reg_change$estimate - reg_change$std.error * 1.96
write.csv(reg_change, path_table("sad_trend_region_change.csv"), row.names = FALSE)

# ---- 4. Fig. 3: annual change per region ---------------------------------------------------------
p_forest <- ggplot(reg_change,
                   aes(x = factor(reg, levels = REGION_ORDER_FOREST), y = estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = estimate - std.error * 1.96, ymax = estimate + std.error * 1.96),
                width = 0.2, linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(x = "Region",
       y = stringr::str_wrap("Annual change in log-transformed alpha value with 95% confidence interval")) +
  theme(legend.position = "none", text = element_text(size = 14)) +
  coord_flip()
ggsave(path_fig("fig3_region_forest.png"), p_forest, width = 6, height = 5, dpi = 300)

# bar-chart versions (two sizes)
reg_change$sname <- gsub("South and Southeast Asia", "S & SE Asia", reg_change$reg)

bars_large <- plot_trend_bars(reg_change, text_size = 16, line_width = 1,   label_width = 10)
ggsave(path_fig("fig3_global_chart.png"),  bars_large, width = 24,   height = 5.64, units = "cm")

bars_small <- plot_trend_bars(reg_change, text_size = 8,  line_width = 0.5, label_width = 8)
ggsave(path_fig("fig3_global_chart2.png"), bars_small, width = 12.5, height = 5,    units = "cm")

# ---- 5. Fig. 3: year-by-year global estimate (ribbon) -----------------------------------------------
# NOTE: this year-by-year estimate uses alpha_raw (alpha_filtered.csv, i.e. before the
# removal of duplicated site-years and Russian sites).
alpha_raw$yrf <- as.factor(alpha_raw$yr)
model4 <- lmer(logmean ~ yrf + (1 | SiteCode), data = alpha_raw, REML = TRUE)
eff_fac <- effects::effect(term = "yrf", mod = model4)
x_urch <- as.data.frame(eff_fac)
x_urch$yr <- as.numeric(as.character(x_urch$yrf))

p_ribbon <- ggplot(data = x_urch, aes(x = yr, y = fit)) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(x = yr, ymin = lower, ymax = upper), fill = REGION_COLOURS[["Global"]],
              alpha = 0.5) +
  labs(x = "Year", y = "Estimate of log10 mean alpha value ± 95%CI") +
  ylim(0.46, 0.56) +
  scale_x_continuous(breaks = c(1990, 2000, 2010, 2018)) +
  theme_classic() +
  theme(legend.position = "none", text = element_text(size = 18))
ggsave(path_fig("fig3_global_ribbon.png"), p_ribbon, width = 7, height = 5, dpi = 300)

# ---- 6. Fig. S5: year-by-year estimate for each region ---------------------------------------------
alpha_clean$yrf <- as.factor(alpha_clean$yr)
plot_list <- list()

for (reg in intersect(REGION_ORDER_PANELS, regions)) {
  temp <- alpha_clean %>% filter(region %in% reg)
  temp$yrf <- droplevels(temp$yrf)

  tmod <- lmer(logmean ~ yrf + (1 | SiteCode), data = temp, REML = TRUE)
  eff_fac2 <- effects::effect(term = "yrf", mod = tmod)

  x_urch2 <- as.data.frame(eff_fac2)
  x_urch2$yr <- as.numeric(as.character(x_urch2$yrf))
  write.csv(x_urch2, path_inter(paste0("sad_yearly_effects_", reg, ".csv")), row.names = FALSE)

  plot_list[[reg]] <- ggplot(data = x_urch2, aes(x = yr, y = fit)) +
    geom_line(linewidth = 0.8) +
    geom_ribbon(aes(x = yr, ymin = lower, ymax = upper), fill = REGION_COLOURS[[reg]],
                alpha = 0.5) +
    labs(x = "Year", y = "Estimate of log10 mean \nalpha value ± 95%CI") +
    scale_x_continuous(breaks = c(1990, 2000, 2010, 2018)) +
    ggtitle(reg) +
    theme_classic() +
    theme(legend.position = "none", text = element_text(size = 12))
}
unplotted <- setdiff(regions, REGION_ORDER_PANELS)
if (length(unplotted) > 0) {
  warning("Regions in the data without a panel: ", paste(unplotted, collapse = ", "))
}

all_reg_plot <- ggarrange(plotlist = plot_list, ncol = 4, nrow = 2)
ggsave(path_fig("fig_S5_regions.png"), plot = all_reg_plot, width = 12, height = 6, dpi = 300)

all_reg_plot_same_y <- ggarrange(plotlist = lapply(plot_list, function(p) p + ylim(0.25, 1.1)),
                                 ncol = 4, nrow = 2)
ggsave(path_fig("fig_S5_regions_same_y.png"), plot = all_reg_plot_same_y,
       width = 12, height = 6, dpi = 300)

save_session_info("07_qii_sad_trends")
