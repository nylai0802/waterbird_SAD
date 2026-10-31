# =============================================================================
# 03_qi_initial_abundance_models.R
# Methods: Question i -> "Relationship between abundance change and initial
#          abundance"  (Table S1; Fig. 2; Extended Data Fig. 3)
#
# Input : output/intermediate/qi_species_change_initial_ra.csv   (from script 02)
# Output: output/tables/qi_models_rate_of_change_*.csv          LMM: rate ~ initial abundance
#         output/tables/qi_models_change_category.csv           binomial GLMMs, one per category
#         output/tables/qi_n_by_bin.csv, qi_n_scatter.csv       sample sizes behind the figures
#         output/figures/fig2_initial_abundance.png, ed_fig3_range.png
#
# =============================================================================

source("R/config.R")
source("R/functions_species_change.R")
source("R/functions_plotting.R")
library(dplyr)
library(lme4)
library(effects)
library(ggplot2)
library(ggpubr)

require_files(path_inter("qi_species_change_initial_ra.csv"))
results_b <- read.csv(path_inter("qi_species_change_initial_ra.csv"))   # finite initial abundance
results_c <- results_b[!is.na(results_b$estimate), ]                    # + a fitted rate of change

# ---- 1. rate of change ~ initial relative abundance (LMM, species random) -----------
# Methods: "linear mixed models with the estimated rate of abundance change as the
# response variable, initial relative abundance as the explanatory variable, and
# species as a random factor"
results_c$Sp2 <- as.factor(results_c$Species)
m_rate <- lmer(estimate ~ relative_abundance + (1 | Sp2), data = results_c, REML = TRUE)
print(summary(m_rate))
print(car::Anova(m_rate))

write.csv(data.frame(term = rownames(coef(summary(m_rate))), coef(summary(m_rate)),
                     row.names = NULL),
          path_table("qi_models_rate_of_change_coefficients.csv"), row.names = FALSE)
write.csv(as.data.frame(car::Anova(m_rate)),
          path_table("qi_models_rate_of_change_anova.csv"))

# fitted line and 95% CI for Fig. 2
effects_rate <- effects::effect(term = "relative_abundance", mod = m_rate)
x_rate <- as.data.frame(effects_rate)

# ---- 2. probability of each type of change ~ initial relative abundance --------------
# One binomial GLMM per category (1 = species is in that category, 0 = not),
# initial relative abundance as the explanatory variable, species as a random factor.
bi_data <- results_b %>% dplyr::select(Species, Category, relative_abundance)
bi_data$Species <- as.factor(bi_data$Species)

categories <- c("Loss"       = "Species loss",
                "Declining"  = "Declining Species",
                "Increasing" = "Increasing Species",
                "No Change"  = "No significant change species")

extract_model_info <- function(model, model_name) {
  coefs <- fixef(model)
  se <- sqrt(diag(vcov(model)))
  z_values <- coefs / se
  p_values <- 2 * (1 - pnorm(abs(z_values)))     # Wald z test
  data.frame(Model = model_name, Term = names(coefs), Estimate = coefs,
             Std_Error = se, z_value = z_values, p_value = p_values, row.names = NULL)
}

binom_results <- list()
for (k in names(categories)) {
  bi_data$y <- as.integer(bi_data$Category == categories[[k]])
  m <- glmer(y ~ relative_abundance + (1 | Species),
             family = binomial(link = "logit"), data = bi_data)
  print(summary(m))
  binom_results[[k]] <- extract_model_info(m, k)
}
binom_results <- do.call(rbind, binom_results)
binom_results$significance <- ifelse(binom_results$p_value < 0.001, "***",
                              ifelse(binom_results$p_value < 0.01, "**",
                              ifelse(binom_results$p_value < 0.05, "*", "")))
binom_results <- binom_results[binom_results$Term == "relative_abundance", ]
print(binom_results)
write.csv(binom_results, path_table("qi_models_change_category.csv"), row.names = FALSE)

# ---- 3. Fig. 2 -----------------------------------------------------------------------
results_d <- results_b %>%
  group_by(pSAD2, Category) %>%
  summarise(sp_no = n(), .groups = "drop")

fig_bar     <- plot_category_bar(results_d)
fig_scatter <- plot_estimate_scatter(results_c, fit = x_rate) + ylim(-0.5, 0.5)

fig2 <- ggarrange(fig_bar, fig_scatter, ncol = 2, nrow = 1,
                  common.legend = TRUE, legend = "bottom")
ggsave(path_fig("fig2_initial_abundance.png"), fig2, width = 10, height = 5, dpi = 300)

# Extended Data Fig. 3: full range of estimates, CI drawn opaque
cat("\nRange of estimated rates of change:", range(results_c$estimate), "\n")
ed3 <- plot_estimate_scatter(results_c, fit = x_rate, ci_alpha = 1, ci_fill = "purple") +
  ylim(-3, 1) + theme(legend.position = "bottom")
ggsave(path_fig("ed_fig3_range.png"), ed3, width = 6, height = 5, dpi = 300)

# ---- 4. sample sizes behind the figures --------------------------------------------------
n_by_bin <- results_b %>%
  group_by(pSAD2) %>%
  summarise(unique_species = n_distinct(Species), total_rows = n(), .groups = "drop")
write.csv(n_by_bin, path_table("qi_n_by_bin.csv"), row.names = FALSE)

n_scatter <- results_c %>% summarise(unique_species = n_distinct(Species), total_rows = n())
write.csv(n_scatter, path_table("qi_n_scatter.csv"), row.names = FALSE)

# category x bin counts
write.csv(xtabs(sp_no ~ Category + pSAD2, data = results_d),
          path_table("qi_category_by_bin.csv"))

save_session_info("03_qi_initial_abundance_models")
