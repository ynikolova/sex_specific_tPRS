# tPRS F X FS longitudinal time  posthoc
#
# Inputs: input_data/predictors_primary.csv, input_data/structural_measures.txt
# Models below retain the original formulas, filtering, and FDR families.
# Run from the repository root; see README.md and docs/review_notes.md.
analysis_id <- "model1/tPRS_F_X_FS_longitudinal_time_posthoc"
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
label <- 'smri_vol_scs_hpuslh'
pheno <- 'SCV'


### load and merge input files
### delta DS for LH SCV in females
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

### FS
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
t2$time <- ifelse(t2$eventname=="baseline_year_1_arm_1",0,1)

# FreeSurfer export positions: 5 = SID (renamed from src_subject_id);
# 7 = interview_age when included; 9 = eventname when included.
# Measurement selection: 11:78, 225:292, 444:447, 451, 452, 455, 462:467, 469.
# The original code describes these collectively as cortical thickness, surface area,
# and subcortical volume measures, excluding sulcal depth. Exact per-position names
# and the metric assigned to each block cannot be verified without the export header.
# Column 486 is expected to be the appended time field; verify after creating time.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
t2 <- t2[,c(5,7,11:78,225:292,444:447,451,452,455,462:467,469,486)] ### selecting relevant columns only (SID, interview_age, thickness, surface area, and subcortical volume (dropped sulcal depth), time)

### merge
df_FS <- merge(t1, t2, by='SID')
### remove one subject from each relative pair
df_FS <- df_FS[df_FS$has_rel==0,]

# winsorize
# Clamp selected columns using DescTools defaults; selected columns are unchanged.
# Merged-table columns 52:201 are the intended 150 MRI outcome fields.
# Positions refer to this merged table, not the original FreeSurfer export; offsets
# vary with predictor selection and inclusion of age/event/time columns.
# Verify with names(df_FS)[c(52:201)] before interpreting the selected measures.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
df_FS[df_FS$time==0,52:201] <- lapply(df_FS[df_FS$time==0,52:201], Winsorize)
# Clamp selected columns using DescTools defaults; selected columns are unchanged.
# Merged-table columns 52:201 are the intended 150 MRI outcome fields.
# Positions refer to this merged table, not the original FreeSurfer export; offsets
# vary with predictor selection and inclusion of age/event/time columns.
# Verify with names(df_FS)[c(52:201)] before interpreting the selected measures.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
df_FS[df_FS$time==1,52:201] <- lapply(df_FS[df_FS$time==1,52:201], Winsorize)

## subset by females
df_FS <- df_FS[df_FS$sex==1,]


### keep only longitudinal data
df_FS <- df_FS[duplicated(df_FS$SID) | duplicated(df_FS$SID, fromLast = TRUE), ]
### order according to time (0 then 1)
df_FS <- df_FS[order(df_FS$time),]

### standardize
# Replace this outcome with standardized residuals from the covariate model.
df_FS[df_FS$time==0,]$smri_vol_scs_hpuslh <- rstandard(lm(df_FS$smri_vol_scs_hpuslh ~ as.numeric(interview_age) + eTIV + eigenvec.1 + eigenvec.2 + eigenvec.3 + eigenvec.4 + eigenvec.5 + eigenvec.6 + eigenvec.7 + eigenvec.8 + eigenvec.9 + eigenvec.10, data=df_FS, na.action=na.exclude, subset= time==0))
# Replace this outcome with standardized residuals from the covariate model.
df_FS[df_FS$time==1,]$smri_vol_scs_hpuslh <- rstandard(lm(df_FS$smri_vol_scs_hpuslh ~ as.numeric(interview_age) + eTIV + eigenvec.1 + eigenvec.2 + eigenvec.3 + eigenvec.4 + eigenvec.5 + eigenvec.6 + eigenvec.7 + eigenvec.8 + eigenvec.9 + eigenvec.10, data=df_FS, na.action=na.exclude, subset= time==1))

### convert rows to the differences between the two times for each SID
# Merged columns 52:201 are intended to contain the 150 MRI outcomes.
# diff computes within-SID differences for each field; verify their names and time order.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
df_FS <- df_FS %>% group_by(SID) %>% mutate_at(c(52:201),diff)
### remove the duplicated row
df_FS <- df_FS[-which(duplicated(df_FS$SID)),]


### create output file
results <- list()

### set sample and print relevant column numbers
### test effects in females
df <- df_FS
sample <- 'female'
colno_scores <- grep(scores, names(df))
colno_regions <- grep(label, names(df))


### test main effects
effect <- 'main'
# Fit the specified mixed model; inspect the random-effects structure in review notes.
    lme <- lme(smri_vol_scs_hpuslh ~ MetaXcan_tPRS_F,
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
# Column 192 is used as a missingness-count field in the posthoc hippocampal analysis.
# Its exact name cannot be verified from the code alone: inspect names(df)[192].
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
    results <- append(results, list(data.frame(score="MetaXcan_tPRS_F", region="smri_vol_scs_hpuslh", sample=sample, effect=effect, t=round(tval,3), p=round(pval,3), n=nrow(df)-sum(is.na(df[,192])))))

### compile outputs
results <- as.data.frame(rbindlist(results))

### apply FDR correction
results[results$sample=='female', "pFDR"] <- round(p.adjust(results[results$sample=='female',]$p, method='fdr'),5)


### save results
# Save under outputs/<script-name>/, using the generic result filename results.csv.
write_results_csv(results, "results.csv", row.names=FALSE)




