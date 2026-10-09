# test tPRS F x FS x CBCL mediation time 2 
#
# Inputs: input_data/symptoms_followup.csv, input_data/predictors_primary.csv, input_data/structural_measures.txt
# Models below retain the original formulas, filtering, and FDR families.
# Run from the repository root; see README.md and docs/review_notes.md.
analysis_id <- "model1/test_tPRS_F_x_FS_x_CBCL_mediation_time_2"
source("R/setup.R")

library(data.table)
library(dplyr)
library(mediation)
library(DescTools)


### set variables
# Paths are resolved relative to the repository root.
score <- 'MetaXcan_tPRS_F'
SCV_label <- 'smri_vol_scs_aar|smri_vol_scs_hpusrh'
CBCL_label <- 'cbcl_scr_syn_withdep_r'

### load data
### load input files
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
t3 <- t2[t2$eventname=="2_year_follow_up_y_arm_1",]
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
t3 <- t3[,c(5,11:78,225:292,444:447,451,452,455,462:467,469)] ### selecting relevant columns only (SID, thickness, surface area, and subcortical volume (dropped sulcal depth))

# Read a private input table; preserve the documented column order.
df1 <- read.csv("input_data/symptoms_followup.csv")
# Columns 1:76 retain the original CBCL/predictor block and exclude later transformed-score columns.
# Exact names of all 76 fields require the original CBCL table header.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
df1 <- df1[,c(1:76)]
## new PCs
# Read a private input table; preserve the documented column order.
df2 <- read.csv("input_data/predictors_primary.csv",h=T)
# Predictor table: column 1 is expected to be SID; columns 97:106 are expected to be
# the ten updated ancestry PCs (eigenvec.1 through eigenvec.10). Verify names in the input header.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
df2 <- df2[,c(1, 97:106)] # only need SID and new PCs
### merge
df1 <- merge(df1, df2, by='SID')

### merge files
t <- merge(df1, t3, by='SID')
### remove one subject from each relative pair
t <- t[t$has_rel==0,]


### create output file
results <- list()

# Clamp selected columns using DescTools defaults; selected columns are unchanged.
# Merged-table columns 87:236 are the intended 150 MRI outcome fields.
# Positions refer to this merged table, not the original FreeSurfer export; offsets
# vary with predictor selection and inclusion of age/event/time columns.
# Verify with names(t)[c(87:236)] before interpreting the selected measures.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
t[,87:236] <- lapply(t[,87:236], Winsorize)


### identify sample
df <- t[t$sex==1,]

### residualize mediator variables
colno_SCV <- grep(SCV_label, names(df))
for (i in colno_SCV){
# Replace this outcome with standardized residuals from the covariate model.
  df[,i] <- rstandard(lm(df[,i] ~ as.numeric(interview_age) + eTIV + eigenvec.1 + eigenvec.2 + eigenvec.3 + eigenvec.4 + eigenvec.5 + eigenvec.6 + eigenvec.7 + eigenvec.8 + eigenvec.9 + eigenvec.10 + site, data=df, na.action=na.exclude))}

colno_CBCL <- grep(CBCL_label, names(df))
for (i in colno_CBCL){
# Replace this outcome with standardized residuals from the covariate model.
  df[,i] <- rstandard(lm(df[,i] ~ as.numeric(interview_age) + eigenvec.1 + eigenvec.2 + eigenvec.3 + eigenvec.4 + eigenvec.5 + eigenvec.6 + eigenvec.7 + eigenvec.8 + eigenvec.9 + eigenvec.10 + site, data=df, na.action=na.exclude))}


# ------------------------------------------------------------
# DV = cbcl_scr_syn_withdep_r
# IV = MetaXcan_tPRS_F
# MED = smri_vol_scs_aar
# ------------------------------------------------------------

label <- 'smri_vol_scs_aar'
scale <- 'cbcl_scr_syn_withdep_r'

### test mediator variable: smri_vol_scs_aar with cbcl_scr_syn_withdep_r
region <- 'smri_vol_scs_aar'
var <- paste0(scale, '|', score, '|', region)
df1 <- df[,grep(var, names(df))]
df1 <- df1[complete.cases(df1),]
fit.med <- lm(smri_vol_scs_aar ~ MetaXcan_tPRS_F, data=df1)
fit.dv <- lm(cbcl_scr_syn_withdep_r ~ MetaXcan_tPRS_F + smri_vol_scs_aar, data=df1)
fit.b <- lm(cbcl_scr_syn_withdep_r ~ smri_vol_scs_aar, data=df1)
mediate <- mediate(fit.med, fit.dv, treat='MetaXcan_tPRS_F', mediator='smri_vol_scs_aar', boot=T)
results <- append(results, list(data.frame(region=region, n=nrow(df1),
# Confidence-interval vector: element 1 is the lower bound; element 2 is the upper bound.
                                           ACME_est=summary(mediate)$d0, ACME_CI_low=summary(mediate)$d0.ci[1], ACME_CI_high=summary(mediate)$d0.ci[2], ACME_p=summary(mediate)$d0.p,
# Confidence-interval vector: element 1 is the lower bound; element 2 is the upper bound.
                                           ADE_est=summary(mediate)$z0, ADE_CI_low=summary(mediate)$z0.ci[1], ADE_CI_high=summary(mediate)$z0.ci[2], ADE_p=summary(mediate)$z0.p,
# Confidence-interval vector: element 1 is the lower bound; element 2 is the upper bound.
                                           TE_est=summary(mediate)$tau.coef, TE_CI_low=summary(mediate)$tau.ci[1], TE_CI_high=summary(mediate)$tau.ci[2], TE_p=summary(mediate)$tau.p,
# lm coefficient matrix: row 2 is the first non-intercept term; column 1 = estimate,
# column 4 = p-value. Legacy a_t/b_t labels store estimates, not t statistics.
# PORTABILITY: The numbered coefficient row may change when model terms change.
# Check matrix row names and adjust the index to select the intended effect.
# Matrix columns refer to the summary statistics described above, not input-data fields.
                                           a_t=summary(fit.med)$coefficients[2,1], a_p=summary(fit.med)$coefficients[2,4],
# lm coefficient matrix: row 2 is the first non-intercept term; column 1 = estimate,
# column 4 = p-value. Legacy a_t/b_t labels store estimates, not t statistics.
# PORTABILITY: The numbered coefficient row may change when model terms change.
# Check matrix row names and adjust the index to select the intended effect.
# Matrix columns refer to the summary statistics described above, not input-data fields.
                                           b_t=summary(fit.b)$coefficients[2,1], b_p=summary(fit.b)$coefficients[2,4])))


# ------------------------------------------------------------
# DV = cbcl_scr_syn_withdep_r
# IV = MetaXcan_tPRS_F
# MED = smri_vol_scs_hpusrh
# ------------------------------------------------------------

label <- 'smri_vol_scs_aar'
scale <- 'cbcl_scr_syn_withdep_r'

### test mediator variable: smri_vol_scs_hpusrh with cbcl_scr_syn_withdep_r
region <- 'smri_vol_scs_hpusrh'
var <- paste0(scale, '|', score, '|', region)
df1 <- df[,grep(var, names(df))]
df1 <- df1[complete.cases(df1),]
fit.med <- lm(smri_vol_scs_hpusrh ~ MetaXcan_tPRS_F, data=df1)
fit.dv <- lm(cbcl_scr_syn_withdep_r ~ MetaXcan_tPRS_F + smri_vol_scs_hpusrh, data=df1)
fit.b <- lm(cbcl_scr_syn_withdep_r ~ smri_vol_scs_hpusrh, data=df1)
mediate <- mediate(fit.med, fit.dv, treat='MetaXcan_tPRS_F', mediator='smri_vol_scs_hpusrh', boot=T)
results <- append(results, list(data.frame(region=region, n=nrow(df1),
# Confidence-interval vector: element 1 is the lower bound; element 2 is the upper bound.
                                           ACME_est=summary(mediate)$d0, ACME_CI_low=summary(mediate)$d0.ci[1], ACME_CI_high=summary(mediate)$d0.ci[2], ACME_p=summary(mediate)$d0.p,
# Confidence-interval vector: element 1 is the lower bound; element 2 is the upper bound.
                                           ADE_est=summary(mediate)$z0, ADE_CI_low=summary(mediate)$z0.ci[1], ADE_CI_high=summary(mediate)$z0.ci[2], ADE_p=summary(mediate)$z0.p,
# Confidence-interval vector: element 1 is the lower bound; element 2 is the upper bound.
                                           TE_est=summary(mediate)$tau.coef, TE_CI_low=summary(mediate)$tau.ci[1], TE_CI_high=summary(mediate)$tau.ci[2], TE_p=summary(mediate)$tau.p,
# lm coefficient matrix: row 2 is the first non-intercept term; column 1 = estimate,
# column 4 = p-value. Legacy a_t/b_t labels store estimates, not t statistics.
# PORTABILITY: The numbered coefficient row may change when model terms change.
# Check matrix row names and adjust the index to select the intended effect.
# Matrix columns refer to the summary statistics described above, not input-data fields.
                                           a_t=summary(fit.med)$coefficients[2,1], a_p=summary(fit.med)$coefficients[2,4],
# lm coefficient matrix: row 2 is the first non-intercept term; column 1 = estimate,
# column 4 = p-value. Legacy a_t/b_t labels store estimates, not t statistics.
# PORTABILITY: The numbered coefficient row may change when model terms change.
# Check matrix row names and adjust the index to select the intended effect.
# Matrix columns refer to the summary statistics described above, not input-data fields.
                                           b_t=summary(fit.b)$coefficients[2,1], b_p=summary(fit.b)$coefficients[2,4])))



### compile outputs
results <- as.data.frame(rbindlist(results))

### apply FDR correction
results$ACME_pFDR <- p.adjust(results$ACME_p, method='fdr')

### save results
# Save under outputs/<script-name>/, using the generic result filename results.csv.
write_results_csv(results, "results.csv", row.names=FALSE)
