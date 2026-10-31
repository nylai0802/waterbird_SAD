# =============================================================================
# config.R
# Shared settings for the waterbird SAD analysis. Sourced by every script.
#
# Run all scripts from the PROJECT ROOT (open waterbird-sad.Rproj in RStudio),
# e.g.   source("scripts/02_qi_species_change.R")
# =============================================================================

## ---- Paths ------------------------------------------------------------------
DIR_DATA   <- "data"                          # input data (not in the repository)
DIR_OUT    <- "output"                        # everything the scripts create
DIR_INTER  <- file.path(DIR_OUT, "intermediate")
DIR_TABLES <- file.path(DIR_OUT, "tables")
DIR_FIGS   <- file.path(DIR_OUT, "figures")
for (d in c(DIR_INTER, DIR_TABLES, DIR_FIGS)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

path_data  <- function(f) file.path(DIR_DATA, f)
path_inter <- function(f) file.path(DIR_INTER, f)
path_table <- function(f) file.path(DIR_TABLES, f)
path_fig   <- function(f) file.path(DIR_FIGS, f)

## ---- Input files (all expected in data/; see data/README.md) -----------------
FILES <- list(
  qi_raw_counts     = "ac_raw_data.csv",                    # Question i/iii: counts
  link_sauer        = "SupplementaryDataS4_390species.csv", # Link & Sauer effort parameters (muP, muB)
  rtm_sites         = "rtm_site.csv",                       # sites used for regression-to-the-mean
  qii_counts        = "w90.csv",                            # Question ii: counts for gambin fitting
  site_country      = "sitelist_country2.csv",              # site -> country / region lookup
  alpha_clean       = "allPn2c_clean.csv",                  # alpha table after clean-up (duplicates, Russia)
  cbc_alpha         = "id_all.csv",                         # alpha per count, with site-year id `sy` (sensitivity tests)
  cbc_effort        = "20221108_Nga_Yee_Lai_CBC_Effort_Many_Types_Report.csv"  # CBC survey effort report
)

## ---- Study settings (Materials and Methods) ---------------------------------
SEED       <- 42           # random seed for simulations and subsampling
YEAR_MIN   <- 1990
YEAR_MAX   <- 2018
SEASON_KEEP <- "DectoFeb"  # value of column Season.x that is kept
                           # "DectoFeb" is a season label; the counts are January-February

# Question i
MIN_INDIVIDUALS_QI  <- 1000   # minimum individuals per count
MIN_COUNTS_PER_SITE <- 5      # minimum number of such counts per site
PERIOD_FIRST <- c(1990, 1999) # "first decade"  (species loss / gain)
PERIOD_LAST  <- c(2010, 2018) # "last decade"   (species loss / gain)
ALPHA_LEVEL  <- 0.05

# Null model (Methods: "Null model for abundance change")
NULL_N_SPECIES   <- 100
NULL_INIT_MEAN   <- 1.44      # mean of log10 initial abundance
NULL_INIT_SD     <- 1.00
NULL_SD_CHANGE   <- c(0.1, 0.2, 0.3)
NULL_N_YEARS     <- 28
NULL_N_SIM       <- 1000

# Question ii (gambin)
GAMBIN_MIN_SPECIES     <- 10
GAMBIN_MIN_INDIVIDUALS <- 1000
GAMBIN_SUBSAMPLE       <- 1000   # individuals drawn per subsample
GAMBIN_N               <- 1000   # number of subsamples per count
GAMBIN_R               <- 3
ALPHA_MAX              <- 15     # counts with mean alpha > 15 are excluded
ALPHA_P_MIN            <- 0.05   # counts with mean p < 0.05 are excluded

## ---- Colours and orderings used in figures ----------------------------------
CATEGORY_COLOURS <- c(
  "Species gain"                  = "#475c92",
  "Species loss"                  = "#551017",
  "Declining Species"             = "#a71111",
  "Increasing Species"            = "#759bc7",
  "No significant change species" = "#bac1ca"
)

REGION_COLOURS <- c(
  "South and Southeast Asia" = "#f5d730",
  "North America"            = "#4c73b5",
  "Europe"                   = "#2f2d72",
  "Oceania"                  = "#3c653a",
  "Middle east"              = "#cc2529",
  "East Asia"                = "#92ccdd",
  "South America"            = "#86bc4f",
  "Africa"                   = "#d2792e",
  "Global"                   = "#5A5A5A"
)

REGION_ORDER_FOREST <- c("Oceania", "South and Southeast Asia", "East Asia",
                         "Middle east", "Africa", "Europe", "South America",
                         "North America", "Global")
REGION_ORDER_PANELS <- c("North America", "South America", "Europe", "Africa",
                         "Middle east", "East Asia", "South and Southeast Asia",
                         "Oceania")

## ---- Small utilities ---------------------------------------------------------
require_files <- function(...) {
  paths   <- c(...)
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop("Missing input file(s):\n  ", paste(missing, collapse = "\n  "),
         "\nSee data/README.md for what each file is.", call. = FALSE)
  }
  invisible(TRUE)
}

save_session_info <- function(script_name) {
  writeLines(capture.output(sessionInfo()),
             file.path(DIR_INTER, paste0("sessionInfo_", script_name, ".txt")))
}
