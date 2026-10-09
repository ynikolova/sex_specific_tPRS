# Shared path setup. Start R in the repository root.
if (!file.exists("README.md") || !dir.exists("scripts")) {
  stop("Run this script from the model1-repository root directory.")
}
output_dir <- file.path("outputs", analysis_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# Preserve CSV contents while isolating each analysis variant.
write_results_csv <- function(x, file, ...) {
  utils::write.csv(x, file = file.path(output_dir, basename(file)), ...)
}

# Several original scripts call these without loading the defining package.
# Bind the original package functions explicitly, without changing calculations.
rbindlist <- data.table::rbindlist
Winsorize <- DescTools::Winsorize
