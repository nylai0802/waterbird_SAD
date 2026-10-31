# =============================================================================
# 10_qiii_sad_vs_assemblage_change.R
# Methods: Question iii -> "Quantifying assemblage and SAD shape changes"
#          (Fig. 4, Extended Data / Fig. S8, S9)
#
# For each site: the proportion of species that were gained, lost,
# increasing or declining (script 02), and the per-site rate of change in
# log10 alpha. The latter is regressed on the four proportions, globally and
# for each region.
#
# Input : output/intermediate/qi_species_change.csv         (script 02)
#         output/intermediate/site_alpha_change.csv  per-site slope of log10 alpha ~ year,
#                                              columns SiteCode, estimate, ..., region (script 09)
# Output: output/tables/qiii_*.csv
#         output/figures/fig4_*.png, fig_S8_site_composition.png, fig_S9_correlations.png
# =============================================================================

source("R/config.R")
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
library(ggpubr)
library(scales)
library(PerformanceAnalytics)

require_files(path_inter("qi_species_change.csv"), path_inter("site_alpha_change.csv"))

change <- read.csv(path_inter("site_alpha_change.csv"))      # alpha change per site
spc    <- read.csv(path_inter("qi_species_change.csv"))     # species-level change

# ---- 1. proportion of species in each category at each site ----------------------------------------
# denominator = all species recorded at the site (including 'no significant change')
category_counts <- spc %>%
  group_by(SiteCode, Category) %>%
  summarise(CategoryCount = n(), .groups = "drop") %>%
  group_by(SiteCode) %>%
  mutate(total_count = sum(CategoryCount),
         prop = CategoryCount / total_count) %>%
  ungroup()

cat("Categories present:\n"); print(unique(category_counts$Category))

# keep the four categories analysed; 'no significant change' is not a model term
change_categories <- c("Species gain", "Species loss", "Declining Species", "Increasing Species")

cc2 <- category_counts %>%
  filter(Category %in% change_categories) %>%
  dplyr::select(SiteCode, Category, prop) %>%
  pivot_wider(names_from = Category, values_from = prop, values_fill = 0)

# As in the original analysis, a site in which NONE of the four categories occurs
# (all its species 'no significant change') has no row here and is not in the models.
cat("Sites without any gained, lost, increasing or declining species (not modelled):",
    n_distinct(spc$SiteCode) - nrow(cc2), "\n")

# ---- 2. collinearity (Fig. S9a) --------------------------------------------------------------------
# Fig. S9a: Pearson correlations between the proportions of gained, lost, increasing and
# declining species at each site (all sites; a category absent from a site has proportion 0)
pw2 <- category_counts %>%
  dplyr::select(SiteCode, Category, prop) %>%
  pivot_wider(names_from = Category, values_from = prop, values_fill = 0)

png(path_fig("fig_S9_correlations.png"), width = 1200, height = 800, units = "px", res = 150)
chart.Correlation(pw2[, change_categories])
dev.off()

# collinearity of the four predictors used in the Question iii model (gain, loss,
# declining, increasing; sites with at least one of them), for the Methods statement
# that no correlation exceeds 0.7
cor_model <- cor(cc2[, change_categories])
print(round(cor_model, 3))
write.csv(cor_model, path_table("qiii_correlations_model_predictors.csv"))
cat("Largest absolute correlation among model predictors:",
    max(abs(cor_model[upper.tri(cor_model)])), "(Methods: none above 0.7)\n")
print(cor.test(cc2$`Species gain`, cc2$`Species loss`))

# ---- 3. join the SAD change estimate ------------------------------------------------------------------
cc2 <- cc2 %>% left_join(change, by = "SiteCode")
cat("Sites with an alpha-change estimate:", sum(!is.na(cc2$estimate)),
    "of", nrow(cc2), "\n")

# ---- 4. global model ------------------------------------------------------------------------------------
new_mod <- glm(estimate ~ `Species gain` + `Species loss` + `Declining Species` +
                 `Increasing Species`, data = cc2)
print(summary(new_mod))

# Fig. S9b: the model without one of the two correlated variables
new_mod2 <- glm(estimate ~ `Species loss` + `Declining Species` + `Increasing Species`, data = cc2)
new_mod3 <- glm(estimate ~ `Species gain` + `Declining Species` + `Increasing Species`, data = cc2)
print(summary(new_mod2))
print(summary(new_mod3))
write.csv(rbind(cbind(model = "without Species gain", tidy(new_mod2)),
                cbind(model = "without Species loss", tidy(new_mod3))),
          path_table("qiii_models_without_gain_or_loss.csv"), row.names = FALSE)

# ---- 5. Fig. 4 scatter panels ---------------------------------------------------------------------------
panel_theme <- theme_classic() +
  theme(axis.text = element_text(size = 12), axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(), axis.title = element_text(size = 12))
y_lab <- "Rate of change in \nlog10 transformed SAD alpha"

scatter_panel <- function(var, colour, xlab, smooth = TRUE) {
  p <- ggplot(cc2, aes(x = .data[[var]], y = estimate)) +
    geom_point(color = colour, size = 0.9)
  if (smooth) p <- p + geom_smooth(method = "lm", col = "black")
  p + geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = xlab, y = y_lab) + panel_theme
}

p1 <- scatter_panel("Species gain",       CATEGORY_COLOURS[["Species gain"]],       "Prop. of species gain")
p2 <- scatter_panel("Species loss",       CATEGORY_COLOURS[["Species loss"]],       "Prop. of species loss")
p3 <- scatter_panel("Declining Species",  CATEGORY_COLOURS[["Declining Species"]],  "Prop. of declining species",
                    smooth = FALSE)       # no fitted line for declining species (as in the original)
p4 <- scatter_panel("Increasing Species", CATEGORY_COLOURS[["Increasing Species"]], "Prop. of increasing species")

fig4_scatter <- ggarrange(p1, p2, p4, p3, ncol = 2, nrow = 2, align = "h")
ggsave(path_fig("fig4_scatter.png"), fig4_scatter, width = 8, height = 7, dpi = 300)

# ---- 6. regional models -------------------------------------------------------------------------------------
cat("\nSites per region:\n"); print(table(cc2$region, useNA = "ifany"))
regions <- unique(cc2$region)

# regional scatter panel (Extended Data). In the original script the x-axis labels of
# the declining / increasing panels were swapped; here each label matches its variable.
regional_panel <- function(dat, var, colour, xlab, smooth) {
  p <- ggplot(dat, aes(x = .data[[var]], y = estimate)) + geom_point(color = colour)
  if (smooth) p <- p + geom_smooth(method = "lm", col = "black")
  p + geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = xlab, y = "Rate of change in SAD alpha") + panel_theme
}

all_models <- list()
for (sreg in regions) {
  rdec <- cc2 %>% filter(region %in% sreg)

  oaplot <- ggarrange(
    regional_panel(rdec, "Species gain",       CATEGORY_COLOURS[["Species gain"]],       "Prop. of species gain",       TRUE),
    regional_panel(rdec, "Species loss",       CATEGORY_COLOURS[["Species loss"]],       "Prop. of species loss",       FALSE),
    regional_panel(rdec, "Declining Species",  CATEGORY_COLOURS[["Declining Species"]],  "Prop. of declining species",  FALSE),
    regional_panel(rdec, "Increasing Species", CATEGORY_COLOURS[["Increasing Species"]], "Prop. of increasing species", FALSE) +
      rremove("x.text"),
    labels = c(paste0(sreg)), ncol = 2, nrow = 2, align = "h",
    font.label = list(size = 10), hjust = -3)
  ggsave(path_fig(paste0("fig4_regional_", gsub(" ", "_", sreg), ".png")), oaplot,
         width = 16, height = 9)

  rmod <- glm(estimate ~ `Species gain` + `Species loss` + `Declining Species` +
                `Increasing Species`, data = rdec)
  print(summary(rmod))
  print(car::Anova(rmod))
  all_models[[sreg]] <- cbind(sreg = sreg, tidy(rmod))
}

# add the global model
all_models[["global"]] <- cbind(sreg = "Global", tidy(new_mod))
atmod <- do.call(rbind, all_models)
atmod$term <- gsub("`", "", atmod$term)           # glm() keeps backticks around non-syntactic names
atmod$sign <- ifelse(atmod$p.value <= 0.05, "significant", "insignificant")
atmod <- atmod %>% filter(term != "(Intercept)")
write.csv(atmod, path_table("qiii_models_global_and_regional.csv"), row.names = FALSE)

# ---- 7. Fig. 4 forest plot (estimates by region) ----------------------------------------------------------------
# All four terms are displayed.
term_levels <- c("Declining Species", "Increasing Species", "Species loss", "Species gain")
aplot <- ggplot(atmod,
                aes(x = factor(sreg, levels = REGION_ORDER_FOREST),
                    y = estimate,
                    shape = factor(term, levels = term_levels),
                    color = factor(sign, levels = c("significant", "insignificant")))) +
  geom_pointrange(aes(ymin = estimate - (std.error * 1.96), ymax = estimate + (std.error * 1.96)),
                  position = position_dodge(0.6)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  scale_color_manual(values = c("#0072b2", "#ddcc77")) +
  scale_shape_manual(values = c(7, 17, 15, 16),
                     labels = c("Species declined", "Species increase", "Species loss", "Species gain")) +
  coord_flip() +
  theme(legend.position = "bottom") +
  labs(y = "Estimate", x = "Region", shape = "") +
  guides(color = guide_none(), shape = guide_legend(nrow = 2)) +
  scale_x_discrete(labels = label_wrap(10))

midp <- ggarrange(p1, p4, ncol = 1, nrow = 2, align = "h")
rigp <- ggarrange(p2, p3, ncol = 1, nrow = 2, align = "h")
fig4 <- ggarrange(aplot, midp, rigp, ncol = 3, nrow = 1, align = "h")
ggsave(path_fig("fig4_combined.png"), fig4, width = 14, height = 7, dpi = 300)

# ---- 8. Fig. S8: composition of change at every site ----------------------------------------------------------------
pw2 <- category_counts %>%
  dplyr::select(SiteCode, Category, prop) %>%
  pivot_wider(names_from = Category, values_from = prop, values_fill = 0)

# every site's proportions must add up to 1
site_total <- rowSums(pw2[, setdiff(names(pw2), "SiteCode")])
stopifnot(all(abs(site_total - 1) < 1e-8))

pw3 <- pw2 %>% pivot_longer(!SiteCode, names_to = "change", values_to = "value")

site_plot <- ggplot(data = pw3,
                    aes(fill = factor(change, levels = c("Species loss", "Declining Species",
                                                         "No significant change species",
                                                         "Increasing Species", "Species gain")),
                        x = SiteCode, y = value)) +
  geom_bar(position = "fill", stat = "identity") +
  theme(axis.text.x = element_blank(), axis.line.x = element_blank()) +
  labs(x = "Site", y = "Proportion of Species") +
  scale_fill_manual(values = CATEGORY_COLOURS,
                    name = "Change",
                    labels = c("species loss", "declining species",
                               "no significant change species", "increasing species",
                               "species gain")) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1))
ggsave(path_fig("fig_S8_site_composition.png"), site_plot, width = 1200, height = 500, units = "px")

save_session_info("10_qiii_sad_vs_assemblage_change")
