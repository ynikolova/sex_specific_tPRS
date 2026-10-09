# Install declared dependencies; original package versions were not supplied.
packages <- c("DescTools", "dplyr", "ggplot2", "lme4", "mediation", "nlme", "zoo")
missing <- setdiff(packages, rownames(installed.packages()))
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
