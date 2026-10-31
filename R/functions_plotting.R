# =============================================================================
# functions_plotting.R
# Plot helpers shared by several scripts.
# Requires ggplot2 (>= 3.4) and dplyr to be attached by the calling script.
# =============================================================================

# Theme used for the null-model figure.
My_Theme <- ggplot2::theme_classic()

# -----------------------------------------------------------------------------
# Stacked proportion bars: category of change by initial relative abundance bin
# df_d needs columns pSAD2, Category, sp_no.
# -----------------------------------------------------------------------------
plot_category_bar <- function(df_d) {
  lev <- c("Species loss", "Declining Species",
           "No significant change species", "Increasing Species")
  df_d <- df_d[df_d$Category %in% lev, ]

  ggplot2::ggplot(df_d, ggplot2::aes(fill = factor(Category, levels = lev),
                                     y = sp_no, x = pSAD2)) +
    ggplot2::geom_bar(position = "fill", stat = "identity") +
    ggplot2::scale_fill_manual(values = CATEGORY_COLOURS, name = "Change") +
    ggplot2::scale_x_discrete(labels = PSAD_LABELS) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Initial relative abundance (%)", y = "Proportion of species") +
    ggplot2::theme(legend.position = "bottom")
}

# -----------------------------------------------------------------------------
# Annual rate of change against initial relative abundance, optionally with the
# fitted mixed-model line and 95% CI (fit = data frame from effects::effect()).
# df_c needs columns Category, estimate, relative_abundance.
# -----------------------------------------------------------------------------
plot_estimate_scatter <- function(df_c, fit = NULL, ci_alpha = 0.3, ci_fill = "black") {
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = df_c,
                        ggplot2::aes(colour = Category, y = estimate, x = relative_abundance)) +
    ggplot2::scale_colour_manual(values = CATEGORY_COLOURS) +
    ggplot2::scale_x_continuous(limits = c(0, 1),
                                labels = scales::label_number(scale = 100)) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Initial relative abundance (%)", y = "Annual rate of change")

  if (!is.null(fit)) {
    p <- p +
      ggplot2::geom_line(data = fit, ggplot2::aes(x = relative_abundance, y = fit),
                         color = "black") +
      ggplot2::geom_ribbon(data = fit,
                           ggplot2::aes(x = relative_abundance, ymin = lower, ymax = upper),
                           alpha = ci_alpha, fill = ci_fill)
  }
  p
}

# -----------------------------------------------------------------------------
# Null-model figure panels (scripts/05_qi_null_model.R).
# data = the `classifications` table of simulate_and_classify_populations5()
# -----------------------------------------------------------------------------
# classes of the null model, in legend order (also the stacking order, top to bottom, in the bars)
NULL_CLASS_LABELS <- c("Extinct"               = "Species loss",
                       "Declined"              = "Declining Species",
                       "No Significant Change" = "No significant change species",
                       "Increased"             = "Increasing Species")

NULL_CLASS_COLOURS <- c("Species loss"                  = "#551017",
                        "Declining Species"             = "#a71111",
                        "No significant change species" = "#bac1ca",
                        "Increasing Species"            = "#759bc7")

NULL_BIN_LABELS <- c("p0_10" = "0-10%",   "p10_20" = "10-20%", "p20_30" = "20-30%",
                     "p30_40" = "30-40%", "p40_50" = "40-50%", "p50_60" = "50-60%",
                     "p60_70" = "60-70%", "p70_80" = "70-80%", "p80_90" = "80-90%",
                     "p90ab"  = "90-100%")

# relabel the classification to the category names used in the main figures
relabel_null_classes <- function(x) {
  factor(NULL_CLASS_LABELS[as.character(x)], levels = unname(NULL_CLASS_LABELS))
}

plot_bar <- function(data) {
  d <- data[is.finite(data$Relative_Abundance), ]
  d$pSAD2 <- bin_relative_abundance(d$Relative_Abundance)
  d$Category <- relabel_null_classes(d$Classification)
  d <- as.data.frame(table(pSAD2 = d$pSAD2, Category = d$Category), responseName = "sp_no")

  ggplot2::ggplot(d, ggplot2::aes(x = pSAD2, y = sp_no, fill = Category)) +
    ggplot2::geom_bar(position = "fill", stat = "identity") +
    ggplot2::scale_fill_manual(values = NULL_CLASS_COLOURS, name = "Change Category",
                               drop = FALSE) +
    ggplot2::scale_x_discrete(labels = NULL_BIN_LABELS) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Initial relative abundance (%)", y = "Proportion of species") +
    ggplot2::theme(legend.position = "bottom")
}

plot_point <- function(data) {
  data$Category <- relabel_null_classes(data$Classification)
  ggplot2::ggplot(data, ggplot2::aes(x = Relative_Abundance, y = Estimate, colour = Category)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.5, na.rm = TRUE) +
    ggplot2::geom_smooth(ggplot2::aes(x = Relative_Abundance, y = Estimate),
                         inherit.aes = FALSE, method = "lm", colour = "black",
                         fill = "grey60", linewidth = 0.6, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = NULL_CLASS_COLOURS, name = "Change Category",
                                 drop = FALSE) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Initial relative abundance (%)", y = "Annual rate of change") +
    ggplot2::theme(legend.position = "bottom")
}

# -----------------------------------------------------------------------------
# Bar chart of annual change in log10 alpha (global + regions) with 95% CI
# reg_change needs columns sname (short region label), reg, estimate, std.error.
# -----------------------------------------------------------------------------
plot_trend_bars <- function(reg_change, text_size, line_width, label_width) {
  lev <- c("Global", "North America", "South America", "Europe", "Africa",
           "Middle east", "East Asia", "S & SE Asia", "Oceania")

  ggplot2::ggplot(reg_change,
                  ggplot2::aes(x = factor(sname, levels = lev), y = estimate, fill = reg)) +
    ggplot2::geom_bar(stat = "identity", show.legend = FALSE, color = "black",
                      linewidth = line_width) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = estimate - std.error * 1.96,
                                        ymax = estimate + std.error * 1.96),
                           width = 0, linewidth = line_width,
                           position = ggplot2::position_dodge(.9)) +
    ggplot2::scale_fill_manual(values = REGION_COLOURS) +
    ggplot2::theme_classic() +
    ggplot2::labs(y = "Annual change in log-transformed alpha value with 95% CI") +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.text.x  = ggplot2::element_text(size = text_size),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.line.y  = ggplot2::element_line(linewidth = line_width),
                   axis.text.y  = ggplot2::element_text(size = text_size),
                   axis.title.y = ggplot2::element_blank(),
                   axis.line.x  = ggplot2::element_blank()) +
    ggplot2::geom_hline(yintercept = 0, linewidth = line_width) +
    ggplot2::ylim(-0.0075, 0.0075) +
    ggplot2::scale_x_discrete(labels = scales::label_wrap(label_width))
}
