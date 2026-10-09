# Model 2: test tPRS F x FS time 2 main in females pFDR families 
# See the shared README and review notes for data requirements and model definitions.
analysis_id <- "model2/test_tPRS_F_x_FS_time_2_main_in_females_pFDR_families"
source("R/setup.R")

library(data.table)
library(dplyr)
library(lme4)
library(nlme)

### set variables
# Paths are relative to the repository root.
scores <- 'MetaXcan_tPRS_F'
CSA_SCV_label <- 'smri_area|smri_vol'
CT_label <- 'smri_thick'
pheno <- 'FS_time_2_main_effects_females_winsorized_MODEL_2'

### load and merge input files

t1 <- read.csv("input_data/predictors_primary.csv",h=T)
# Predictor table: retain columns 1:50 (SID, original covariates and selected scores),
# plus 97:106 (the ten updated ancestry PCs). Individual names in 1:50 were not supplied.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
t1 <- t1[,c(1:50, 97:106)] # we don't need all the other scores for these analyses 

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
t2 <- t2[t2$eventname=="2_year_follow_up_y_arm_1",]
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
t2 <- t2[,c(5,7,11:78,225:292,444:447,451,452,455,462:467,469)] ### selecting relevant columns only (SID, interview_age, thickness, surface area, and subcortical volume (dropped sulcal depth))

df <- merge(t1, t2, by='SID')

### remove one subject from each relative pair
df <- df[df$has_rel==0,]

### create output file
results <- list()

# Merged-table columns 62:211 are the intended 150 MRI outcome fields.
# Positions refer to this merged table, not the original FreeSurfer export; offsets
# vary with predictor selection and inclusion of age/event/time columns.
# Verify with names(df)[c(62:211)] before interpreting the selected measures.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
df[,62:211] <- lapply(df[,62:211], Winsorize)

### set sample and print relevant column numbers
sample <- 'females'
colno_scores <- grep(scores, names(df))
colno_CSA_SCV_regions <- grep(CSA_SCV_label, names(df))
colno_CT_regions <- grep(CT_label, names(df))

### setting sample to females only
df <- df[df$sex==1,]

### test main effects in females
effect <- 'main'

### for SCV and CSA (eTIV as a covariate)
for (i in colno_scores){
  for (k in colno_CSA_SCV_regions){
    tprs <- df[,i]
    region <- df[,k]
    ctrl <- lmeControl(opt='optim')
    lme <- lme(region ~ tprs + as.numeric(interview_age) + eTIV + eigenvec.1 + eigenvec.2 + eigenvec.3 + eigenvec.4 + eigenvec.5 + eigenvec.6 + eigenvec.7 + eigenvec.8 + eigenvec.9 + eigenvec.10 + MetaXcan_tPRS_M, 
               data=df, random=~1|site, na.action=na.omit, control=ctrl)
    sum <- summary(lme)
# nlme summary tTable: row 2 is intended to select the first non-intercept effect
# (e.g., tPRS, time, sex, or the selected MRI predictor, as shown in the model formula);
# column 4 is the t statistic. Check rownames(sum$tTable) for the term.
# PORTABILITY: The numbered coefficient row may change when model terms change.
# Check matrix row names and adjust the index to select the intended effect.
# Matrix columns refer to the summary statistics described above, not input-data fields.
    tval <- sum$tTable[2,4]
# nlme summary tTable: row 2 is intended to select the first non-intercept effect
# (e.g., tPRS, time, sex, or the selected MRI predictor, as shown in the model formula);
# column 5 is the p-value. Check rownames(sum$tTable) for the term.
# PORTABILITY: The numbered coefficient row may change when model terms change.
# Check matrix row names and adjust the index to select the intended effect.
# Matrix columns refer to the summary statistics described above, not input-data fields.
    pval <- sum$tTable[2,5]
    results <- append(results, list(data.frame(score=names(df)[i], region=names(df)[k], sample=sample, effect=effect, t=round(tval,3), p=round(pval,3), n=nrow(df)-sum(is.na(df[,i])))))}}

### for CT (eTIV not included as a covariate)
for (i in colno_scores){
  for (k in colno_CT_regions){
    tprs <- df[,i]
    region <- df[,k]
    ctrl <- lmeControl(opt='optim')
    lme <- lme(region ~ tprs + as.numeric(interview_age) + eigenvec.1 + eigenvec.2 + eigenvec.3 + eigenvec.4 + eigenvec.5 + eigenvec.6 + eigenvec.7 + eigenvec.8 + eigenvec.9 + eigenvec.10 + MetaXcan_tPRS_M, 
               data=df, random=~1|site, na.action=na.omit, control=ctrl)
    sum <- summary(lme)
# nlme summary tTable: row 2 is intended to select the first non-intercept effect
# (e.g., tPRS, time, sex, or the selected MRI predictor, as shown in the model formula);
# column 4 is the t statistic. Check rownames(sum$tTable) for the term.
# PORTABILITY: The numbered coefficient row may change when model terms change.
# Check matrix row names and adjust the index to select the intended effect.
# Matrix columns refer to the summary statistics described above, not input-data fields.
    tval <- sum$tTable[2,4]
# nlme summary tTable: row 2 is intended to select the first non-intercept effect
# (e.g., tPRS, time, sex, or the selected MRI predictor, as shown in the model formula);
# column 5 is the p-value. Check rownames(sum$tTable) for the term.
# PORTABILITY: The numbered coefficient row may change when model terms change.
# Check matrix row names and adjust the index to select the intended effect.
# Matrix columns refer to the summary statistics described above, not input-data fields.
    pval <- sum$tTable[2,5]
    results <- append(results, list(data.frame(score=names(df)[i], region=names(df)[k], sample=sample, effect=effect, t=round(tval,3), p=round(pval,3), n=nrow(df)-sum(is.na(df[,i])))))}}


### compile outputs
results <- as.data.frame(rbindlist(results))

### apply FDR correction across all separate families of analyses 
results[grep("smri_area_cdk_",results$region), "pFDR"] <- round(p.adjust(results[grep("smri_area_cdk_",results$region),]$p, method='fdr'),3)
results[grep("smri_vol",results$region), "pFDR"] <- round(p.adjust(results[grep("smri_vol",results$region),]$p, method='fdr'),3)
results[grep("smri_thick",results$region), "pFDR"] <- round(p.adjust(results[grep("smri_thick",results$region),]$p, method='fdr'),3)


### save results
write_results_csv(results, "results.csv", row.names=FALSE)

