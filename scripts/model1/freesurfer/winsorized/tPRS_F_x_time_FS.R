# tPRS F x time FS 
#
# Inputs: input_data/predictors_primary.csv, input_data/structural_measures.txt
# Models below retain the original formulas, filtering, and FDR families.
# Run from the repository root; see README.md and docs/review_notes.md.
analysis_id <- "model1/tPRS_F_x_time_FS"
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
CSA_SCV_label <- 'smri_area|smri_vol'
CT_label <- 'smri_thick'
pheno <- 'FS_intx_winsorized'

### load and merge input files
# Read a private input table; preserve the documented column order.
t1 <- read.csv("input_data/predictors_primary.csv",h=T)
# Predictor table: retain columns 1:50 (SID, original covariates and selected scores),
# plus 97:106 (the ten updated ancestry PCs). Individual names in 1:50 were not supplied.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
t1 <- t1[,c(1:50, 97:106)] # we don't need all the other scores for these analyses

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

t2$time <- ifelse(t2$eventname=="baseline_year_1_arm_1",0,1)

### keep only longitudinal data
t2 <- t2[duplicated(t2$SID) | duplicated(t2$SID, fromLast = TRUE), ]

### order according to time (0 then 1)
t2 <- t2[order(t2$time),]

df <- merge(t1, t2, by='SID')

### remove one subject from each relative pair
df <- df[df$has_rel==0,]

### create output file
results <- list()

# winsorize
# Clamp selected columns using DescTools defaults; selected columns are unchanged.
# Merged-table columns 63:212 are the intended 150 MRI outcome fields.
# Positions refer to this merged table, not the original FreeSurfer export; offsets
# vary with predictor selection and inclusion of age/event/time columns.
# Verify with names(df)[c(63:212)] before interpreting the selected measures.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
df[df$time==0,63:212] <- lapply(df[df$time==0,63:212], Winsorize)
# Clamp selected columns using DescTools defaults; selected columns are unchanged.
# Merged-table columns 63:212 are the intended 150 MRI outcome fields.
# Positions refer to this merged table, not the original FreeSurfer export; offsets
# vary with predictor selection and inclusion of age/event/time columns.
# Verify with names(df)[c(63:212)] before interpreting the selected measures.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
df[df$time==1,63:212] <- lapply(df[df$time==1,63:212], Winsorize)

### set sample and print relevant column numbers
### test effects in females
df <- df[df$sex==1,]
sample <- 'female'
colno_scores <- grep(scores, names(df))
colno_CSA_SCV_regions <- grep(CSA_SCV_label, names(df))
colno_CT_regions <- grep(CT_label, names(df))


for (i in colno_CSA_SCV_regions){
# Replace this outcome with standardized residuals from the covariate model.
  df[,i] <- rstandard(lm(df[,i] ~ as.numeric(interview_age) + eTIV + eigenvec.1 + eigenvec.2 + eigenvec.3 + eigenvec.4 + eigenvec.5 + eigenvec.6 + eigenvec.7 + eigenvec.8 + eigenvec.9 + eigenvec.10, data=df, na.action=na.exclude))}

for (i in colno_CT_regions){
# Replace this outcome with standardized residuals from the covariate model.
  df[,i] <- rstandard(lm(df[,i] ~ as.numeric(interview_age) + eigenvec.1 + eigenvec.2 + eigenvec.3 + eigenvec.4 + eigenvec.5 + eigenvec.6 + eigenvec.7 + eigenvec.8 + eigenvec.9 + eigenvec.10, data=df, na.action=na.exclude))}

### test interaction effects
effect <- 'intx'

### for SCV and CSA (eTIV as a covariate)
for (i in colno_scores){
  for (k in colno_CSA_SCV_regions){
    tprs <- df[,i]
    region <- df[,k]
# Fit the specified mixed model; inspect the random-effects structure in review notes.
    lme <- lme(region ~ tprs*time,
               data=df, random= list(SID=~time,site=~1), control = lmeControl(opt = "optim"), na.action=na.omit)
    sum <- summary(lme)
# nlme summary tTable: row 4 is intended to select the interaction effect
# (tPRS-by-time or sex-by-time, as shown in the preceding model formula);
# column 4 is the t statistic. Check rownames(sum$tTable) for the term.
# PORTABILITY: The numbered coefficient row may change when model terms change.
# Check matrix row names and adjust the index to select the intended effect.
# Matrix columns refer to the summary statistics described above, not input-data fields.
    tval <- sum$tTable[4,4]
# Extract the original coefficient by row index; verify its term when changing a formula.
# nlme summary tTable: row 4 is intended to select the interaction effect
# (tPRS-by-time or sex-by-time, as shown in the preceding model formula);
# column 5 is the p-value. Check rownames(sum$tTable) for the term.
# PORTABILITY: The numbered coefficient row may change when model terms change.
# Check matrix row names and adjust the index to select the intended effect.
# Matrix columns refer to the summary statistics described above, not input-data fields.
    pval <- sum$tTable[4,5]
    results <- append(results, list(data.frame(score=names(df)[i], region=names(df)[k], sample=sample, effect=effect, t=round(tval,3), p=round(pval,3), n=nrow(df)-sum(is.na(df[,i])))))}}

### for SCV and CSA (eTIV as a covariate)
for (i in colno_scores){
  for (k in colno_CT_regions){
    tprs <- df[,i]
    region <- df[,k]
# Fit the specified mixed model; inspect the random-effects structure in review notes.
    lme <- lme(region ~ tprs*time,
               data=df, random= list(SID=~time,site=~1), control = lmeControl(opt = "optim"), na.action=na.omit)
    sum <- summary(lme)
# nlme summary tTable: row 4 is intended to select the interaction effect
# (tPRS-by-time or sex-by-time, as shown in the preceding model formula);
# column 4 is the t statistic. Check rownames(sum$tTable) for the term.
# PORTABILITY: The numbered coefficient row may change when model terms change.
# Check matrix row names and adjust the index to select the intended effect.
# Matrix columns refer to the summary statistics described above, not input-data fields.
    tval <- sum$tTable[4,4]
# Extract the original coefficient by row index; verify its term when changing a formula.
# nlme summary tTable: row 4 is intended to select the interaction effect
# (tPRS-by-time or sex-by-time, as shown in the preceding model formula);
# column 5 is the p-value. Check rownames(sum$tTable) for the term.
# PORTABILITY: The numbered coefficient row may change when model terms change.
# Check matrix row names and adjust the index to select the intended effect.
# Matrix columns refer to the summary statistics described above, not input-data fields.
    pval <- sum$tTable[4,5]
    results <- append(results, list(data.frame(score=names(df)[i], region=names(df)[k], sample=sample, effect=effect, t=round(tval,3), p=round(pval,3), n=nrow(df)-sum(is.na(df[,i])))))}}


### compile outputs
results <- as.data.frame(rbindlist(results))

#### apply FDR correction across all separate families of analyses
results[grep("smri_area_cdk_",results$region), "pFDR"] <- round(p.adjust(results[grep("smri_area_cdk_",results$region),]$p, method='fdr'),3)
results[grep("smri_vol",results$region), "pFDR"] <- round(p.adjust(results[grep("smri_vol",results$region),]$p, method='fdr'),3)
results[grep("smri_thick",results$region), "pFDR"] <- round(p.adjust(results[grep("smri_thick",results$region),]$p, method='fdr'),3)

### save results
# Save under outputs/<script-name>/, using the generic result filename results.csv.
write_results_csv(results, "results.csv", row.names=FALSE)

