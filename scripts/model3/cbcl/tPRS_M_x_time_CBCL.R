# Model 3: tPRS M x time CBCL 
# See the shared README and review notes for data requirements and model definitions.
analysis_id <- "model3/tPRS_M_x_time_CBCL"
source("R/setup.R")

library(data.table)
library(dplyr)
library(lme4)
library(nlme)
library(zoo)
library(DescTools)


### set variables
# Paths are relative to the repository root.
scores <- 'MetaXcan_tPRS_M'
CBCL_label <- 'cbcl_scr_syn_anxdep_r|cbcl_scr_syn_withdep_r'
pheno <- 'longitudinal_CBCL_separate_t1_t2_effects_MODEL_3'

### load input files
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
df <- rbind(t1,t2)
nrow(df) # 9138 (4569 x2)

### set sample and print relevant column numbers
### test effects in males
df <- df[df$sex==0,]
sample <- 'male'
colno_scores <- grep(scores, names(df))
colno_CBCL_regions <- grep(CBCL_label, names(df))

for (i in colno_CBCL_regions){
  df[,i] <- rstandard(lm(df[,i] ~ as.numeric(interview_age) + eigenvec.1 + eigenvec.2 + eigenvec.3 + eigenvec.4 + eigenvec.5 + eigenvec.6 + eigenvec.7 + eigenvec.8 + eigenvec.9 + eigenvec.10 + MetaXcan_tPRS_F + TPRS_v1_unweighted + PGC_PRS, data=df, na.action=na.exclude))}

### two Time Points
df_time_1 <- df[df$time=="1",]
df_time_2 <- df[df$time=="2",]

### create output file
results <- list()

### test effect AT T1
effect <- 'Baseline'

### for CBCL 
for (i in colno_scores){
  for (k in colno_CBCL_regions){
    tprs <- df_time_1[,i]
    region <- df_time_1[,k]
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

### for CBCL 
for (i in colno_scores){
  for (k in colno_CBCL_regions){
    tprs <- df_time_2[,i]
    region <- df_time_2[,k]
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
write_results_csv(results, "results.csv", row.names=FALSE)

