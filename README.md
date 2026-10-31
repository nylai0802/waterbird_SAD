# Code for: Three decades of directional change in global waterbird assemblage structure

Analysis code accompanying the manuscript "Three decades of directional change in global waterbird assemblage structure" by Nga-Yee Lai, Sam Nicol, Richard A. Fuller, Taej Mundkur, Moriaki Yasuhara and Tatsuya Amano (manuscript in preparation; the full reference will be added once published).

Archived release: https://doi.org/10.5281/zenodo.23830973

The study uses yearly waterbird counts from the International Waterbird Census (IWC) and the
Christmas Bird Count (CBC), 1990-2018, to ask (i) how the abundance of initially common and rare
species changed, (ii) how the shape of the species abundance distribution (SAD, parameter α of
the gambin model) changed over time, and (iii) how SAD shape change relates to species gain,
loss, increase and decline at each site.

## Repository layout

```
R/                      shared code, sourced by the scripts
  config.R              paths, thresholds, seeds, colours (edit settings here)
  functions_species_change.R   species-level change and initial relative abundance
  functions_null_model.R       random-walk null model
  functions_plotting.R         plot helpers
scripts/                numbered scripts, in the order of the Materials and Methods
data/                   input data (NOT in the repository, see data/README.md)
output/                 created by the scripts (not tracked by git)
```

## How to run

1. Put the input files in `data/` (see `data/README.md`).
2. Open `waterbird-sad.Rproj` in RStudio, or start R in the project root. Every script must be run
   from the project root.
3. Install the packages (below), then run the scripts in numerical order, e.g.
   `source("scripts/02_qi_species_change.R")`.

### Packages

`dplyr`, `tidyr`, `MASS`, `lme4`, `car`, `effects`, `broom`, `broom.mixed`, `ggplot2` (>= 3.5.0),
`ggpubr`, `patchwork`, `scales`, `stringr`, `PerformanceAnalytics`, `gambin`.

Every script saves `sessionInfo()` to `output/intermediate/` so the exact versions are recorded.
For a fully pinned environment, run `renv::init()` once and commit the resulting `renv.lock`.

## Materials and Methods → scripts

| Methods section | Script | Main outputs |
|---|---|---|
| Waterbird count data; data filters of Question i | `01_check_input_files.R` | checks that the input files are in `data/` |
| Question i: defining site-level species changes; initial commonness and rarity | `02_qi_species_change.R` | species-by-site change categories and initial relative abundance |
| Question i: relationship between abundance change and initial abundance | `03_qi_initial_abundance_models.R` | Table S1, Fig. 2, Extended Data Fig. 3 |
| Question i: testing for regression-to-the-mean | `04_qi_regression_to_mean.R` | Fig. S4 |
| Question i: null model for abundance change | `05_qi_null_model.R` | Fig. S1 |
| Question ii: estimating SAD shape using the gambin model | `06_qii_gambin_alpha.R` | α for every site and year |
| Question ii: gambin results, exclusion rules | `06b_qii_extract_alpha.R` | counts kept for the trend analysis |
| Question ii: global and regional SAD shape trends | `07_qii_sad_trends.R` | Fig. 3, Fig. S5 |
| Question ii: sensitivity tests (survey effort) | `08_qii_sensitivity_survey_effort.R` | Fig. S6 |
| Question ii: sensitivity tests (taxonomic resolution) | `08b_qii_sensitivity_taxonomic_resolution.R` | Fig. S7 |
| Question iii: SAD shape change at each site | `09_qiii_site_alpha_change.R` | per-site slope of log10 α |
| Question iii: relating SAD change to assemblage change | `10_qiii_sad_vs_assemblage_change.R` | Fig. 4, Figs. S8 and S9 |

Settings that appear in the Methods (year range, thresholds, subsample sizes, null-model
parameters) are set in `R/config.R`.

## Reproducibility notes

* Random seeds are set in `config.R` (`SEED`) for the null model and the gambin subsampling.
* The raw count data are not redistributed here; see the data availability statement below.

## Data availability
The waterbird count data used in this study are collated and managed by Wetlands International and the National Audubon Society.

## Licence and citation

Code is released under the MIT licence (`LICENSE`). If you use it, please cite the paper and the
archived release (`CITATION.cff`).
