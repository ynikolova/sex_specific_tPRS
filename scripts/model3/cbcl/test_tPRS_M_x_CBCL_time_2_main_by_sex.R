# Model 3: test tPRS M x CBCL time 2 main by sex 
# See the shared README and review notes for data requirements and model definitions.
analysis_id <- "model3/test_tPRS_M_x_CBCL_time_2_main_by_sex"
source("R/setup.R")

library(data.table)
library(dplyr)
library(lme4)
library(nlme)
library(DescTools)

### set variables
# Paths are relative to the repository root.
scores <- 'MetaXcan_tPRS_M'
label <- 'cbcl_scr_syn_anxdep_r|cbcl_scr_syn_withdep_r'
pheno <- 'CBCL_time_2_main_effects_by_sex_MODEL_3'


### load input files (already merged)
df <- read.csv('input_data/symptoms_followup.csv')

# Columns 1:76 retain the original CBCL/predictor block and exclude later transformed-score columns.
# Exact names of all 76 fields require the original CBCL table header.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
df <- df[,c(1:76)] ## don't need transformed scores

## new PCs
t1 <- read.csv("input_data/predictors_primary.csv",h=T)
# Predictor table: column 1 is expected to be SID; columns 97:106 are expected to be
# the ten updated ancestry PCs (eigenvec.1 through eigenvec.10). Verify names in the input header.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
t1 <- t1[,c(1, 97:106)] # only need SID and new PCs 


### merge 
df <- merge(df, t1, by='SID')

### create output file
results <- list()

### test main effects in females
df1 <- df[df$sex==1,]
sample <- 'female'
colno_scores <- grep(scores, names(df1))
colno_regions <- grep(label, names(df1))
effect <- 'main'
for (i in colno_scores){
  for (k in colno_regions){
    tprs <- df1[,i]
    region <- df1[,k]
    lme <- lme(region ~ tprs + as.numeric(interview_age) + eigenvec.1 + eigenvec.2 + eigenvec.3 + eigenvec.4 + eigenvec.5 + eigenvec.6 + eigenvec.7 + eigenvec.8 + eigenvec.9 + eigenvec.10 + MetaXcan_tPRS_F + TPRS_v1_unweighted + PGC_PRS, 
               data=df1, random=~1|site, na.action=na.omit)
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
    results <- append(results, list(data.frame(score=names(df1)[i], region=names(df1)[k], sample=sample, effect=effect, t=round(tval,3), p=round(pval,3), n=nrow(df1)-sum(is.na(df1[,i])))))}}

### test main effects in males
df1 <- df[df$sex==0,]
sample <- 'male'
colno_scores <- grep(scores, names(df1))
colno_regions <- grep(label, names(df1))
effect <- 'main'
for (i in colno_scores){
  for (k in colno_regions){
    tprs <- df1[,i]
    region <- df1[,k]
    lme <- lme(region ~ tprs + as.numeric(interview_age) + eigenvec.1 + eigenvec.2 + eigenvec.3 + eigenvec.4 + eigenvec.5 + eigenvec.6 + eigenvec.7 + eigenvec.8 + eigenvec.9 + eigenvec.10 + MetaXcan_tPRS_F + TPRS_v1_unweighted + PGC_PRS, 
               data=df1, random=~1|site, na.action=na.omit)
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
    results <- append(results, list(data.frame(score=names(df1)[i], region=names(df1)[k], sample=sample, effect=effect, t=round(tval,3), p=round(pval,3), n=nrow(df1)-sum(is.na(df1[,i])))))}}

### compile outputs
results <- as.data.frame(rbindlist(results))

### apply FDR correction
results[results$sample=='female', "pFDR"] <- round(p.adjust(results[results$sample=='female',]$p, method='fdr'),3)
results[results$sample=='male', "pFDR"] <- round(p.adjust(results[results$sample=='male',]$p, method='fdr'),3)

### save results
write_results_csv(results, "results.csv", row.names=FALSE)
