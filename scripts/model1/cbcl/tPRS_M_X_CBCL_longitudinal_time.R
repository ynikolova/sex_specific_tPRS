# tPRS M X CBCL longitudinal time 
#
# Inputs: input_data/symptoms_followup.csv, input_data/symptoms_baseline.csv, input_data/predictors_primary.csv
# Models below retain the original formulas, filtering, and FDR families.
# Run from the repository root; see README.md and docs/review_notes.md.
analysis_id <- "model1/tPRS_M_X_CBCL_longitudinal_time"
source("R/setup.R")

library(data.table)
library(dplyr)
library(lme4)
library(nlme)
library(zoo)
library(ggplot2)


### set variables
# Paths are resolved relative to the repository root.
scores <- 'MetaXcan_tPRS_M'
pheno <- 'CBCL'

### load and merge input files
### load input files
# Read a private input table; preserve the documented column order.
t1 <- read.csv("input_data/symptoms_baseline.csv") #5002
# Columns 1:76 retain the original CBCL/predictor block and exclude later transformed-score columns.
# Exact names of all 76 fields require the original CBCL table header.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
t1<- t1[,c(1:76)]
## new PCs
# Read a private input table; preserve the documented column order.
df <- read.csv("input_data/predictors_primary.csv",h=T)
# Predictor table: column 1 is expected to be SID; columns 97:106 are expected to be
# the ten updated ancestry PCs (eigenvec.1 through eigenvec.10). Verify names in the input header.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
df <- df[,c(1, 97:106)] # only need SID and new PCs
### merge
t1 <- merge(t1, df, by='SID')


# Read a private input table; preserve the documented column order.
t2 <- read.csv("input_data/symptoms_followup.csv") #4569
# Columns 1:76 retain the original CBCL/predictor block and exclude later transformed-score columns.
# Exact names of all 76 fields require the original CBCL table header.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
t2<- t2[,c(1:76)]
## new PCs
# Read a private input table; preserve the documented column order.
df <- read.csv("input_data/predictors_primary.csv",h=T)
# Predictor table: column 1 is expected to be SID; columns 97:106 are expected to be
# the ten updated ancestry PCs (eigenvec.1 through eigenvec.10). Verify names in the input header.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
df <- df[,c(1, 97:106)] # only need SID and new PCs
### merge
t2 <- merge(t2, df, by='SID')

t1$time <- 1
t2$time <- 2

## keeping t1 data only for those with t2 data
t1 <- t1[c(which(t1$SID %in% t2$SID)),]
nrow(t1) #4569

## binding both datasets (time 1 and time 2)
df_cbcl <- rbind(t1,t2)
nrow(df) # 9138 (4569 x2)

### remove one subject from each relative pair
df_cbcl <- df_cbcl[df_cbcl$has_rel==0,]

### subset by males
df_cbcl <- df_cbcl[df_cbcl$sex==0,]

### keep only longitudinal data
df_cbcl <- df_cbcl[duplicated(df_cbcl$SID) | duplicated(df_cbcl$SID, fromLast = TRUE), ]
### order according to time (0 then 1)
df_cbcl <- df_cbcl[order(df_cbcl$time),]

### standardize
# Replace this outcome with standardized residuals from the covariate model.
df_cbcl[df_cbcl$time==1,]$cbcl_scr_syn_withdep_r <- rstandard(lm(df_cbcl$cbcl_scr_syn_withdep_r ~ as.numeric(interview_age) + eigenvec.1 + eigenvec.2 + eigenvec.3 + eigenvec.4 + eigenvec.5 + eigenvec.6 + eigenvec.7 + eigenvec.8 + eigenvec.9 + eigenvec.10, data=df_cbcl, na.action=na.exclude, subset= time==1))
# Replace this outcome with standardized residuals from the covariate model.
df_cbcl[df_cbcl$time==2,]$cbcl_scr_syn_withdep_r <- rstandard(lm(df_cbcl$cbcl_scr_syn_withdep_r ~ as.numeric(interview_age) + eigenvec.1 + eigenvec.2 + eigenvec.3 + eigenvec.4 + eigenvec.5 + eigenvec.6 + eigenvec.7 + eigenvec.8 + eigenvec.9 + eigenvec.10, data=df_cbcl, na.action=na.exclude, subset= time==2))

### convert rows to the differences between the two times for each SID
# Merged CBCL column 76 is intended to be cbcl_scr_syn_withdep_r (withdrawn/depressed score).
# diff computes within-SID differences in current row order; verify the column name and time order.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
df_cbcl <- df_cbcl %>% group_by(SID) %>% mutate_at(76,diff)
### remove the duplicated row
df_cbcl <- df_cbcl[-which(duplicated(df_cbcl$SID)),]


### create output file
results <- list()


### set sample and print relevant column numbers
### test effects in males
df <- df_cbcl
sample <- 'male'



### test main effects
effect <- 'main'

# Fit the specified mixed model; inspect the random-effects structure in review notes.
  lme <- lme(cbcl_scr_syn_withdep_r ~ MetaXcan_tPRS_M,
             data=df, random=~1|site, na.action=na.omit)
  sum <- summary(lme)
# nlme summary tTable: row 2 is intended to select the first non-intercept effect
# (e.g., tPRS, time, sex, or the selected MRI predictor, as shown in the model formula);
# column 4 is the t statistic. Check rownames(sum$tTable) for the term.
# PORTABILITY: The numbered coefficient row may change when model terms change.
# Check matrix row names and adjust the index to select the intended effect.
# Matrix columns refer to the summary statistics described above, not input-data fields.
  tval <- sum$tTable[2,4]
# Extract the original coefficient by row index; verify its term when changing a formula.
# nlme summary tTable: row 2 is intended to select the first non-intercept effect
# (e.g., tPRS, time, sex, or the selected MRI predictor, as shown in the model formula);
# column 5 is the p-value. Check rownames(sum$tTable) for the term.
# PORTABILITY: The numbered coefficient row may change when model terms change.
# Check matrix row names and adjust the index to select the intended effect.
# Matrix columns refer to the summary statistics described above, not input-data fields.
  pval <- sum$tTable[2,5]
# Column 76 is treated as the CBCL withdrawn/depressed outcome for missingness counting;
# confirm names(df)[76]. This count is not necessarily the fitted-model sample size.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
  results <- append(results, list(data.frame(score="MetaXcan_tPRS_M", dep="cbcl_scr_syn_withdep_r", sample=sample, effect=effect, t=round(tval,3), p=round(pval,3), n=nrow(df)-sum(is.na(df[,76])))))


### compile outputs
results <- as.data.frame(rbindlist(results))

### apply FDR correction
results[results$sample=='male', "pFDR"] <- round(p.adjust(results[results$sample=='male',]$p, method='fdr'),3)


### save results
# Save under outputs/<script-name>/, using the generic result filename results.csv.
write_results_csv(results, "results.csv", row.names=FALSE)

