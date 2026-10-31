# =============================================================================
# 01_check_input_files.R
# Checks that every input file listed in R/config.R is present in data/.
# The data are not part of this repository (see README.md, "Data availability").
# =============================================================================

source("R/config.R")

for (f in unlist(FILES)) {
  cat(sprintf("%-60s %s\n", f, if (file.exists(path_data(f))) "found" else "MISSING"))
}
cat(sprintf("%-60s %s\n", "alpha_unid/ (gambin output, script 08b)",
            if (dir.exists(file.path(DIR_DATA, "alpha_unid"))) "found" else "MISSING"))

save_session_info("01_check_input_files")
