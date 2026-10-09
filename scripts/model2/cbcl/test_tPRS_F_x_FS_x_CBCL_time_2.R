# Model 2: test tPRS F x FS x CBCL time 2 
# See the shared README and review notes for data requirements and model definitions.
analysis_id <- "model2/test_tPRS_F_x_FS_x_CBCL_time_2"
source("R/setup.R")

  library(data.table)
  library(dplyr)
  library(lme4)
  library(nlme)
  library(DescTools)
  
  
  ### set variables
# Paths are relative to the repository root.
  scores <- 'MetaXcan_tPRS_F'
  label_vol <- 'smri_vol_scs_hpuslh|smri_vol_scs_hpusrh|smri_vol_scs_caudaterh|smri_vol_scs_aar'   
  label_depr <- 'cbcl_scr_syn_anxdep_r|cbcl_scr_syn_withdep_r'
  
  ### load input files
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
  t3 <- t3[,c(5,7,11:78,225:292,444:447,451,452,455,462:467,469)] ### selecting relevant columns only (SID, interview_age, thickness, surface area, and subcortical volume (dropped sulcal depth))
  
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

# Merged-table columns 88:237 are the intended 150 MRI outcome fields.
# Positions refer to this merged table, not the original FreeSurfer export; offsets
# vary with predictor selection and inclusion of age/event/time columns.
# Verify with names(t)[c(88:237)] before interpreting the selected measures.
# PORTABILITY: These numbered data columns refer to the original table layout.
# Their positions may vary with data release, predictor version, column order,
# and preceding merges/selections, and will likely need adjustment for your inputs.
# Check the intended fields described above against names() on this exact table;
# verify the selected names before running, and update the indices as needed.
t[,88:237] <- lapply(t[,88:237], Winsorize)

### test main effect on volume in female subsample
df <- t[t$sex==1,]
sample <- 'female'
colno_vol <- grep(label_vol, names(df))
colno_depr <- grep(label_depr, names(df))


lme <- lme(cbcl_scr_syn_anxdep_r ~ smri_vol_scs_hpuslh + as.numeric(interview_age.x) + eTIV + eigenvec.1 + eigenvec.2 + eigenvec.3 + eigenvec.4 + eigenvec.5 + eigenvec.6 + eigenvec.7 + eigenvec.8 + eigenvec.9 + eigenvec.10 + MetaXcan_tPRS_M,
        data=df, random=~1|site, na.action=na.omit)
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
  results <- append(results, list(data.frame(sample=sample, dep="cbcl_scr_syn_anxdep_r",region="smri_vol_scs_hpuslh", t=round(tval,3), p=round(pval,3), n=nrow(df)-sum(is.na(df[,grep("cbcl_scr_syn_anxdep_r", names(df))])))))


  lme <- lme(cbcl_scr_syn_withdep_r ~ smri_vol_scs_hpuslh + as.numeric(interview_age.x) + eTIV + eigenvec.1 + eigenvec.2 + eigenvec.3 + eigenvec.4 + eigenvec.5 + eigenvec.6 + eigenvec.7 + eigenvec.8 + eigenvec.9 + eigenvec.10 + MetaXcan_tPRS_M,
             data=df, random=~1|site, na.action=na.omit)
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
  results <- append(results, list(data.frame(sample=sample, dep="cbcl_scr_syn_withdep_r",region="smri_vol_scs_hpuslh", t=round(tval,3), p=round(pval,3), n=nrow(df)-sum(is.na(df[,grep("cbcl_scr_syn_withdep_r", names(df))])))))


   lme <- lme(cbcl_scr_syn_anxdep_r ~ smri_vol_scs_caudaterh + as.numeric(interview_age.x) + eTIV + eigenvec.1 + eigenvec.2 + eigenvec.3 + eigenvec.4 + eigenvec.5 + eigenvec.6 + eigenvec.7 + eigenvec.8 + eigenvec.9 + eigenvec.10 + MetaXcan_tPRS_M,
             data=df, random=~1|site, na.action=na.omit)
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
  results <- append(results, list(data.frame(sample=sample, dep="cbcl_scr_syn_anxdep_r",region="smri_vol_scs_caudaterh", t=round(tval,3), p=round(pval,3), n=nrow(df)-sum(is.na(df[,grep("cbcl_scr_syn_anxdep_r", names(df))])))))



  lme <- lme(cbcl_scr_syn_withdep_r ~ smri_vol_scs_caudaterh + as.numeric(interview_age.x) + eTIV + eigenvec.1 + eigenvec.2 + eigenvec.3 + eigenvec.4 + eigenvec.5 + eigenvec.6 + eigenvec.7 + eigenvec.8 + eigenvec.9 + eigenvec.10 + MetaXcan_tPRS_M,
             data=df, random=~1|site, na.action=na.omit)
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
  results <- append(results, list(data.frame(sample=sample, dep="cbcl_scr_syn_withdep_r", region="smri_vol_scs_caudaterh", t=round(tval,3), p=round(pval,3), n=nrow(df)-sum(is.na(df[,grep("cbcl_scr_syn_withdep_r", names(df))])))))



  lme <- lme(cbcl_scr_syn_anxdep_r ~ smri_vol_scs_hpusrh + as.numeric(interview_age.x) + eTIV + eigenvec.1 + eigenvec.2 + eigenvec.3 + eigenvec.4 + eigenvec.5 + eigenvec.6 + eigenvec.7 + eigenvec.8 + eigenvec.9 + eigenvec.10 + MetaXcan_tPRS_M,
             data=df, random=~1|site, na.action=na.omit)
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
  results <- append(results, list(data.frame(sample=sample, dep="cbcl_scr_syn_anxdep_r",region="smri_vol_scs_hpusrh", t=round(tval,3), p=round(pval,3), n=nrow(df)-sum(is.na(df[,grep("cbcl_scr_syn_anxdep_r", names(df))])))))
  
  
  lme <- lme(cbcl_scr_syn_withdep_r ~ smri_vol_scs_hpusrh + as.numeric(interview_age.x) + eTIV + eigenvec.1 + eigenvec.2 + eigenvec.3 + eigenvec.4 + eigenvec.5 + eigenvec.6 + eigenvec.7 + eigenvec.8 + eigenvec.9 + eigenvec.10 + MetaXcan_tPRS_M,
             data=df, random=~1|site, na.action=na.omit)
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
  results <- append(results, list(data.frame(sample=sample, dep="cbcl_scr_syn_withdep_r",region="smri_vol_scs_hpusrh", t=round(tval,3), p=round(pval,3), n=nrow(df)-sum(is.na(df[,grep("cbcl_scr_syn_withdep_r", names(df))])))))
  
  
  lme <- lme(cbcl_scr_syn_anxdep_r ~ smri_vol_scs_aar + as.numeric(interview_age.x) + eTIV + eigenvec.1 + eigenvec.2 + eigenvec.3 + eigenvec.4 + eigenvec.5 + eigenvec.6 + eigenvec.7 + eigenvec.8 + eigenvec.9 + eigenvec.10 + MetaXcan_tPRS_M,
             data=df, random=~1|site, na.action=na.omit)
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
  results <- append(results, list(data.frame(sample=sample, dep="cbcl_scr_syn_anxdep_r",region="smri_vol_scs_aar", t=round(tval,3), p=round(pval,3), n=nrow(df)-sum(is.na(df[,grep("cbcl_scr_syn_anxdep_r", names(df))])))))
  
  
  lme <- lme(cbcl_scr_syn_withdep_r ~ smri_vol_scs_aar + as.numeric(interview_age.x) + eTIV + eigenvec.1 + eigenvec.2 + eigenvec.3 + eigenvec.4 + eigenvec.5 + eigenvec.6 + eigenvec.7 + eigenvec.8 + eigenvec.9 + eigenvec.10 + MetaXcan_tPRS_M,
             data=df, random=~1|site, na.action=na.omit)
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
  results <- append(results, list(data.frame(sample=sample, dep="cbcl_scr_syn_withdep_r",region="smri_vol_scs_aar", t=round(tval,3), p=round(pval,3), n=nrow(df)-sum(is.na(df[,grep("cbcl_scr_syn_withdep_r", names(df))])))))
  

### compile outputs
results <- as.data.frame(rbindlist(results))

### apply FDR correction
results[results$sample=='female', "pFDR"] <- p.adjust(results[results$sample=='female',]$p, method='fdr')

### save results
write_results_csv(results, "results.csv", row.names=FALSE)

