# tPRS F x time FS separate T1 T2 
#
# Inputs: input_data/predictors_primary.csv, input_data/structural_measures.txt
# Models below retain the original formulas, filtering, and FDR families.
# Run from the repository root; see README.md and docs/review_notes.md.
analysis_id <- "model1/tPRS_F_x_time_FS_separate_T1_T2"
source("R/setup.R")

library(data.table)
library(dplyr)
library(lme4)
library(nlme)
library(zoo)
library(DescTools)


### set variables
# Paths are resolved relative to the repository root.
scores <- 'MetaXcan_tPRS_F'
SCV_label <- 'smri_vol_scs_hpuslh|smri_vol_scs_hpusrh|smri_vol_scs_caudaterh|smri_vol_scs_aar'
pheno <- 'longitudinal_FS_winsorized_separate_t1_t2_effects'

### load and merge input files
# Read a private input table; preserve the documented column order.
t1 <- read.csv("input_data/predictors_primary.csv",h=T)
# Predictor table: retain original fields 1:33 and 44:50, excluding fields 34:43;
# append updated ancestry PCs from 97:106. Exact field names in the retained blocks
# cannot be verified without the original predictor-table header.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
t1 <- t1[,c(1:33, 44:50, 97:106)] # we don't need all the other scores for these analyses

# Read a private input table; preserve the documented column order.
t2 <- read.table('input_data/structural_measures.txt',h=T)
# INPUT FORMAT: -1 removes the first data row, expected to contain field descriptions.
# Its presence can vary by export; adjust/remove this step if row 1 is a participant record.
t2 <- t2[-1,] # remove row 1 as it is a description of the column
names(t2)[names(t2) == 'src_subject_id'] <- 'SID'

# FreeSurfer export: columns 11:484 are treated as numeric measurement fields.
# These refer to the original export order, before selecting columns or merging.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
t2 <- t2 %>% mutate_at(c(11:484), as.numeric) ## convert from chr to num
# FreeSurfer export positions: 5 = SID (renamed from src_subject_id);
# 7 = interview_age when included; 9 = eventname when included.
# Measurement selection: 11:78, 225:292, 444:447, 451, 452, 455, 462:467, 469.
# The original code describes these collectively as cortical thickness, surface area,
# and subcortical volume measures, excluding sulcal depth. Exact per-position names
# and the metric assigned to each block cannot be verified without the export header.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
t2 <- t2[,c(5,7,9,11:78,225:292,444:447,451,452,455,462:467,469)] ### selecting relevant columns only (SID, interview_age, eventname, thickness, surface area, and subcortical volume (dropped sulcal depth))

t3 <- merge(t1, t2, by='SID')
t3 <- t3[t3$has_rel==0,]

t3$time <- ifelse(t3$eventname=="baseline_year_1_arm_1",0,1)

### keep only longitudinal data
t3 <- t3[duplicated(t3$SID) | duplicated(t3$SID, fromLast = TRUE), ]

### order according to time (0 then 1)
df <- t3[order(t3$time),]

# Clamp selected columns using DescTools defaults; selected columns are unchanged.
# Merged-table columns 53:202 are the intended 150 MRI outcome fields.
# Positions refer to this merged table, not the original FreeSurfer export; offsets
# vary with predictor selection and inclusion of age/event/time columns.
# Verify with names(df)[c(53:202)] before interpreting the selected measures.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
df[df$time==0,53:202] <- lapply(df[df$time==0,53:202], Winsorize)
# Clamp selected columns using DescTools defaults; selected columns are unchanged.
# Merged-table columns 53:202 are the intended 150 MRI outcome fields.
# Positions refer to this merged table, not the original FreeSurfer export; offsets
# vary with predictor selection and inclusion of age/event/time columns.
# Verify with names(df)[c(53:202)] before interpreting the selected measures.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
df[df$time==1,53:202] <- lapply(df[df$time==1,53:202], Winsorize)

### set sample and print relevant column numbers
### test effects in females
df <- df[df$sex==1,]
sample <- 'female'
colno_scores <- grep(scores, names(df))
colno_SCV_regions <- grep(SCV_label, names(df))

for (i in colno_SCV_regions){
# Replace this outcome with standardized residuals from the covariate model.
  df[,i] <- rstandard(lm(df[,i] ~ as.numeric(interview_age) + eTIV + eigenvec.1 + eigenvec.2 + eigenvec.3 + eigenvec.4 + eigenvec.5 + eigenvec.6 + eigenvec.7 + eigenvec.8 + eigenvec.9 + eigenvec.10, data=df, na.action=na.exclude))}

### two Time Points
df_time_1 <- df[df$time=="0",]
df_time_2 <- df[df$time=="1",]

### create output file
results <- list()

### test effect AT T1
effect <- 'Baseline'

### for SCV
for (i in colno_scores){
  for (k in colno_SCV_regions){
    tprs <- df_time_1[,i]
    region <- df_time_1[,k]
# Fit the specified mixed model; inspect the random-effects structure in review notes.
    lme <- lme(region ~ tprs,
               data=df_time_1, random=~1|site, control = lmeControl(opt = "optim"), na.action=na.omit)
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
    results <- append(results, list(data.frame(score=names(df_time_1)[i], region=names(df_time_1)[k], sample=sample, effect=effect, t=round(tval,3), p=round(pval,3), n=nrow(df_time_1)-sum(is.na(df_time_1[,i])))))}}

### test effect AT T2
effect <- 'Follow-Up'

### for SCV
for (i in colno_scores){
  for (k in colno_SCV_regions){
    tprs <- df_time_2[,i]
    region <- df_time_2[,k]
# Fit the specified mixed model; inspect the random-effects structure in review notes.
    lme <- lme(region ~ tprs,
               data=df_time_2, random=~1|site, control = lmeControl(opt = "optim"), na.action=na.omit)
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
    results <- append(results, list(data.frame(score=names(df_time_2)[i], region=names(df_time_2)[k], sample=sample, effect=effect, t=round(tval,3), p=round(pval,3), n=nrow(df_time_2)-sum(is.na(df_time_2[,i])))))}}


### compile outputs
results <- as.data.frame(rbindlist(results))

### save results
# Save under outputs/<script-name>/, using the generic result filename results.csv.
write_results_csv(results, "results.csv", row.names=FALSE)

