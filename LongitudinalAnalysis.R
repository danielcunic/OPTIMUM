library(dplyr)
library(tidyverse)
##stats packages
library(lme4)
library(buildmer)
library(perm)
library(coin)
library(nlme)

##read in files
rh_stats = read.csv("C:/Users/Daniel Cunic/OneDrive/Documents/OPTIMUM/data/rh_stats_F.csv")
lh_stats = read.csv("C:/Users/Daniel Cunic/OneDrive/Documents/OPTIMUM/data/lh_stats_F.csv")
aseg_stats = read.csv("C:/Users/Daniel Cunic/OneDrive/Documents/OPTIMUM/data/aseg_stats_F-DESKTOP-HUO6GG3.csv")
qc = read_tsv("C:/Users/Daniel Cunic/OneDrive/Documents/OPTIMUM/data/QC_ABCD_outputs.txt") 
npsych = read.csv("C:/Users/Daniel Cunic/OneDrive/Documents/OPTIMUM/data/dbo_DKEFS_RBANS_merged20240228_date051624.csv")

##clean qc table and rename lh 
qc_F <- qc |> rename(Subject.id = "Subject-id") |> select("Subject.id", "SHN", "Pass/Fail", "Re-run", "Doubtful")

##clean np table
npsych_F <- npsych

ix= npsych_F$codeCDATE.03vsCDATE.01<=0;
npsych_F$codeCDATE.03vsCDATE.01 [ix]=NA;
ix= npsych_F$dspanCDATE2_vsCDATE1<=0;
npsych_F$dspanCDATE2_vsCDATE1[ix]=NA;
##npsych_F$zeroDATE=zeros(length(npsych_F$dspanCDATE2_vsCDATE1),1);
npsych_F$dspanCDATE2_vsCDATE1==0;
npsych_F$codeCDATE_03vsCDATE_01;
npsych_F$dspanCDATE2_vsCDATE1;
npsych_F$dspanCDATE2_vsCDATE1<=0;
ix= npsych_F$ED>=80;
npsych_F$ED[ix]=NA;

ix= npsych_F$CDRSCORE_01>3
npsych_F$CDRSCORE_01[ix]=NA;
ix= npsych_F$CDRSCORE_01<0
npsych_F$CDRSCORE_01[ix]=NA;

ix= npsych_F$CDRSCORE_02>3;
npsych_F$CDRSCORE_02[ix]=NA;
ix= npsych_F$CDRSCORE_02<0;
npsych_F$CDRSCORE_02[ix]=NA;

ix= npsych_F$CDRSCORE_03>3;
npsych_F$CDRSCORE_03[ix]=NA;
ix= npsych_F$CDRSCORE_03<0;
npsych_F$CDRSCORE_03[ix]=NA;

##merging tables
lh_rh <- merge(rh_stats, lh_stats, by="Subject.id")
aseg_bothhemis <- merge(lh_rh, aseg_stats, by="Subject.id")
MRI_QC_Long <- merge(aseg_bothhemis, qc_F, by="Subject.id")
MRI_QC_Npsych_Long <- merge(MRI_QC_Long, npsych_F, by="ID_complete")

##Remove all failed qc participants
##identify the rows with a Pass/Fail value of 0, identify the associated ID_completes, turn char vector into data frame
IDof_Failed_QCs <- MRI_QC_Npsych_Long$ID_complete[MRI_QC_Npsych_Long$`Pass/Fail`== 0];
IDof_Failed_QCs <- data.frame(Name=IDof_Failed_QCs);
IDof_Failed_QCs <- IDof_Failed_QCs |> rename(ID_complete = "Name")
##delete rows in MRI_QC_Npsych_Long which are common in the two tables
MRI_QC_Npsych_Long <- anti_join(MRI_QC_Npsych_Long, IDof_Failed_QCs)

#########################
#####     CDR 1     #####
#########################

##clean table for analysis with CDR Score 1
IDof_Missing_CDR1 <- MRI_QC_Npsych_Long$ID_complete[is.na(MRI_QC_Npsych_Long$CDRSCORE_01)];
IDof_Missing_CDR1 <- data.frame(Name=IDof_Missing_CDR1);
IDof_Missing_CDR1 <- IDof_Missing_CDR1 |> rename(ID_complete = "Name")
Tbl4CDR1 <- anti_join(MRI_QC_Npsych_Long, IDof_Missing_CDR1, by = "ID_complete")

##finding all participants w/ a single session
P_singleSesh <- summarise(group_by(Tbl4CDR1, ID_complete), repeats = n())
P_singleSesh <- P_singleSesh |> filter(repeats == 1)

##remove participants with a single session
Tbl4CDR1 <- anti_join(Tbl4CDR1, P_singleSesh)

  ##index time_days column to change odd rows to 0
##odd_ix <- seq(1, 358, 2)
##MRI_QC_Npsych_Long$time_days[odd_ix] = 0

  ##make participants w/ no dates or neg dates values = NA
  #ix_neg_0 <- MRI_QC_Npsych_Long$time_days < 0;  
  #MRI_QC_Npsych_Long$time_days[ix_neg_0] = NA;

##create new column for time in years to make analysis more clear
#MRI_QC_Npsych_Long$time_years <- MRI_QC_Npsych_Long$time_days/365
even_time_ix <- seq(2, 322, 2)
Tbl4CDR1$time_days[even_time_ix] = 0.5
odd_time_ix <- seq(1, 322, 2)
Tbl4CDR1$time_days[odd_time_ix] = 0

##rename column
Tbl4CDR1 <- Tbl4CDR1 |> rename(time_years = "time_days")

##to clean dates and be consistent (not required)
#odd_ix <- seq(1, 364, 2)
#MRI_QC_Npsych_Long$imindex_CDATE_02vsCDATE_01[odd_ix] = 0
#MRI_QC_Npsych_Long$imindex_CDATE_03vsCDATE_01[odd_ix] = 0

###########################
## Stats Experimentation ##
###########################

#baseline <- lmer(rh_bankssts_thickness ~ AGE + GENDER + ED + time_days + (1+time_days|ID_complete), data=MRI_QC_Npsych_Long)

baseline<-lme(Left.VentralDC ~ AGE + GENDER + ED + time_years*CDRSCORE_01, random = ~ time_years | ID_complete, data=Tbl4CDR1)
summary(baseline)

baseline<-lme(lh_entorhinal_thickness ~ AGE + GENDER + ED + time_years, random = ~ time_years | ID_complete, data=CDR_0)
summary(baseline)

baseline<-lme(Right.Hippocampus ~ AGE + GENDER + ED + CDRSCORE_01*time_years, random = ~ time_years | ID_complete, data=Tbl4CDR1)
summary(baseline)

baseline<-lme(rh_precuneus_thickness ~ AGE + GENDER + ED + CDRSCORE_01*time_years, random = ~ time_years | ID_complete, data=Tbl4CDR1)
summary(baseline)

#########################
####  R Entorhinal   ####
#########################

##create table w/ baseline entorhinal, change in.., and CDR score
odd_ix <- seq(1, 322, 2)
even_ix <- seq(2, 322, 2)
delta_entorhinal <- Tbl4CDR1$rh_entorhinal_thickness[odd_ix] - Tbl4CDR1$rh_entorhinal_thickness[even_ix]
baseline_entorhinal <- Tbl4CDR1$rh_entorhinal_thickness[odd_ix]
cdr<-Tbl4CDR1$CDRSCORE_01[odd_ix]

##turning these tables into data frames
baseline_entorhinal <- data.frame(Name=baseline_entorhinal)
delta_entorhinal <- data.frame(Name=delta_entorhinal)
cdr <- data.frame(Name=cdr)

##adding ID_complete column
baseline_entorhinal$ID_complete <- Tbl4CDR1$ID_complete[odd_ix]
delta_entorhinal$ID_complete <- Tbl4CDR1$ID_complete[odd_ix]
cdr$ID_complete <- Tbl4CDR1$ID_complete[odd_ix]

##merge tables into final table for CDR1 analysis
F_Tbl4CDR1 <- merge(baseline_entorhinal, delta_entorhinal, by = "ID_complete")
F_Tbl4CDR1 <- merge(F_Tbl4CDR1, cdr, by = "ID_complete")

##rename columns
F_Tbl4CDR1 <- F_Tbl4CDR1 |> rename(b_r_entorhinal = "Name.x", d_r_entorhinal = "Name.y", CDR = "Name")

##add time_years column?
even_time_ix <- seq(2, 322, 2)
F_Tbl4CDR1$time_years <- Tbl4CDR1$time_years[even_time_ix]

##reorder columns
F_Tbl4CDR1 <- F_Tbl4CDR1 |> relocate(CDR:time_years, .after = ID_complete)

##Changing CDR column to char vectors in order to colour code graphs
ix = F_Tbl4CDR1$CDR == 0.0;
F_Tbl4CDR1$CDR[ix] = "No Impairment";

ix = F_Tbl4CDR1$CDR == 0.5;
F_Tbl4CDR1$CDR[ix] = "Questionable or Mild Impairment";


ix = F_Tbl4CDR1$CDR == 1.0;
F_Tbl4CDR1$CDR[ix] = "Questionable or Mild Impairment"

##visualization with box and violin 
ggplot(data = F_Tbl4CDR1, mapping = aes(x=CDR, y=d_r_entorhinal, fill=CDR)) + 
  geom_violin() + 
  geom_boxplot() + scale_fill_manual(values= c("#8cc5e3", "#ea801c")) +
  labs(title = "Predicting Change in R.H. Entorhinal Thickness Based on CDR Score", x = "CDR Score", y = "Change in L.H. Entorhinal Thickness", fill = "Cognitive Impairment")

##visualization; scatter plot of entorhinal thickness change over time according to CDR; needs work
ggplot(data=CDR_0.5_1, mapping = aes(y=lh_entorhinal_thickness, x=ED, fill=CDRSCORE_01)) +
  geom_point() +
  labs(title = "Predicting Change in R.H. Entorhinal Thickness Based on CDR Score", x = "Change in R.H. Entorhinal Thickness", fill = "Cognitive Impairment")

#########################
####  L Entorhinal   ####
#########################

##create table w/ baseline entorhinal, change in.., and CDR score
odd_ix <- seq(1, 322, 2)
even_ix <- seq(2, 322, 2)
d_l_entorhinal <- Tbl4CDR1$lh_entorhinal_thickness[odd_ix] - Tbl4CDR1$lh_entorhinal_thickness[even_ix]
b_l_entorhinal <- Tbl4CDR1$lh_entorhinal_thickness[odd_ix]

##turning these tables into data frames
b_l_entorhinal <- data.frame(Name=b_l_entorhinal)
d_l_entorhinal <- data.frame(Name=d_l_entorhinal)

##adding ID_complete column
b_l_entorhinal$ID_complete <- Tbl4CDR1$ID_complete[odd_ix]
d_l_entorhinal$ID_complete <- Tbl4CDR1$ID_complete[odd_ix]

##merge tables into final table for CDR1 analysis
F_Tbl4CDR1 <- merge(F_Tbl4CDR1, b_l_entorhinal, by = "ID_complete")
F_Tbl4CDR1 <- merge(F_Tbl4CDR1, d_l_entorhinal, by = "ID_complete")

##rename columns
F_Tbl4CDR1 <- F_Tbl4CDR1 |> rename(b_l_entorhinal = "Name.x", d_l_entorhinal = "Name.y")

##reorder
F_Tbl4CDR1 <- F_Tbl4CDR1 |> relocate(b_l_entorhinal:d_l_entorhinal, .after = d_r_entorhinal)

##visualization with box and violin 
ggplot(data = F_Tbl4CDR1, mapping = aes(x=CDR, y=d_l_entorhinal, fill=CDR)) + 
  geom_violin() + 
  geom_boxplot() + scale_fill_manual(values= c("#8cc5e3", "#ea801c")) +
  labs(title = "Predicting Change in L.H. Entorhinal Thickness Based on CDR Score", x = "CDR Score", y = "Change in L.H. Entorhinal Thickness", fill = "Cognitive Impairment")


###########################
## R Posterior Cingulate ##
###########################

##create table w/ baseline entorhinal, change in.., and CDR score
odd_ix <- seq(1, 322, 2)
even_ix <- seq(2, 322, 2)
d_pcingulate <- Tbl4CDR1$rh_posteriorcingulate_thickness[odd_ix] - Tbl4CDR1$rh_posteriorcingulate_thickness[even_ix]
b_pcingulate <- Tbl4CDR1$rh_posteriorcingulate_thickness[odd_ix]

##turning these tables into data frames
b_pcingulate <- data.frame(Name=b_pcingulate)
d_pcingulate <- data.frame(Name=d_pcingulate)

##adding ID_complete column
b_pcingulate$ID_complete <- Tbl4CDR1$ID_complete[odd_ix]
d_pcingulate$ID_complete <- Tbl4CDR1$ID_complete[odd_ix]

##merge tables into final table for CDR1 analysis
F_Tbl4CDR1 <- merge(F_Tbl4CDR1, b_pcingulate, by = "ID_complete")
F_Tbl4CDR1 <- merge(F_Tbl4CDR1, d_pcingulate, by = "ID_complete")

##rename columns
F_Tbl4CDR1 <- F_Tbl4CDR1 |> rename(b_r_pcingulate = "Name.x", d_r_pcingulate = "Name.y")

##Visualization with violin + box again for posterior cingulate
ggplot(data = F_Tbl4CDR1, mapping = aes(x=CDR, y=d_r_pcingulate, fill=CDR)) + 
  geom_violin() + 
  geom_boxplot() + scale_fill_manual(values= c("#8cc5e3", "#ea801c")) +
  labs(title = "Predicting Change in R.H. Posterior Cingulate Thickness Based on CDR Score", x = "CDR Score", y = "Change in RH. P. Cingulate Thickness", fill = "Cognitive Impairment")

###########################
## L Posterior Cingulate ##
###########################

##create table w/ baseline entorhinal, change in.., and CDR score
odd_ix <- seq(1, 322, 2)
even_ix <- seq(2, 322, 2)
d_l_pcingulate <- Tbl4CDR1$lh_posteriorcingulate_thickness[odd_ix] - Tbl4CDR1$lh_posteriorcingulate_thickness[even_ix]
b_l_pcingulate <- Tbl4CDR1$lh_posteriorcingulate_thickness[odd_ix]

##turning these tables into data frames
b_l_pcingulate <- data.frame(Name=b_l_pcingulate)
d_l_pcingulate <- data.frame(Name=d_l_pcingulate)

##adding ID_complete column
b_l_pcingulate$ID_complete <- Tbl4CDR1$ID_complete[odd_ix]
d_l_pcingulate$ID_complete <- Tbl4CDR1$ID_complete[odd_ix]

##merge tables into final table for CDR1 analysis
F_Tbl4CDR1 <- merge(F_Tbl4CDR1, b_l_pcingulate, by = "ID_complete")
F_Tbl4CDR1 <- merge(F_Tbl4CDR1, d_l_pcingulate, by = "ID_complete")

##rename columns
F_Tbl4CDR1 <- F_Tbl4CDR1 |> rename(b_l_pcingulate = "Name.x", d_l_pcingulate = "Name.y")

##Visualization with violin + box again for posterior cingulate
ggplot(data = F_Tbl4CDR1, mapping = aes(x=CDR, y=d_l_pcingulate, fill=CDR)) + 
  geom_violin() + 
  geom_boxplot() + scale_fill_manual(values= c("#8cc5e3", "#ea801c")) +
  labs(title = "Predicting Change in L.H. Posterior Cingulate Thickness Based on CDR Score", x = "CDR Score", y = "Change in L.H. P. Cingulate Thickness", fill = "Cognitive Impairment")


#########################
####   R Precuneus   ####
#########################

##create table w/ baseline precuneus, change in.., and CDR score
odd_ix <- seq(1, 322, 2)
even_ix <- seq(2, 322, 2)
d_r_precuneus <- Tbl4CDR1$rh_precuneus_thickness[odd_ix] - Tbl4CDR1$rh_precuneus_thickness[even_ix]
b_r_precuneus <- Tbl4CDR1$rh_precuneus_thickness[odd_ix]

##turning these tables into data frames
b_r_precuneus <- data.frame(Name=b_r_precuneus)
d_r_precuneus <- data.frame(Name=d_r_precuneus)

##adding ID_complete column
b_r_precuneus$ID_complete <- Tbl4CDR1$ID_complete[odd_ix]
d_r_precuneus$ID_complete <- Tbl4CDR1$ID_complete[odd_ix]

##merge tables into final table for CDR1 analysis
F_Tbl4CDR1 <- merge(F_Tbl4CDR1, b_r_precuneus, by = "ID_complete")
F_Tbl4CDR1 <- merge(F_Tbl4CDR1, d_r_precuneus, by = "ID_complete")

##rename columns
F_Tbl4CDR1 <- F_Tbl4CDR1 |> rename(b_r_precuneus = "Name.x", d_r_precuneus = "Name.y")

##Visualization with violin + box again for posterior cingulate
ggplot(data = F_Tbl4CDR1, mapping = aes(x=CDR, y=d_r_precuneus, fill=CDR)) + 
  geom_violin() + 
  geom_boxplot() + scale_fill_manual(values= c("#8cc5e3", "#ea801c")) +
  labs(title = "Predicting Change in R.H. Precuneus Thickness Based on CDR Score", x = "CDR Score", y = "Change in R.H. Precuneus Thickness", fill = "Cognitive Impairment")

#########################
####   L Precuneus   ####
#########################

##create table w/ baseline precuneus, change in.., and CDR score
odd_ix <- seq(1, 322, 2)
even_ix <- seq(2, 322, 2)
d_l_precuneus <- Tbl4CDR1$lh_precuneus_thickness[odd_ix] - Tbl4CDR1$lh_precuneus_thickness[even_ix]
b_l_precuneus <- Tbl4CDR1$lh_precuneus_thickness[odd_ix]

##turning these tables into data frames
b_l_precuneus <- data.frame(Name=b_l_precuneus)
d_l_precuneus <- data.frame(Name=d_l_precuneus)

##adding ID_complete column
b_l_precuneus$ID_complete <- Tbl4CDR1$ID_complete[odd_ix]
d_l_precuneus$ID_complete <- Tbl4CDR1$ID_complete[odd_ix]

##merge tables into final table for CDR1 analysis
F_Tbl4CDR1 <- merge(F_Tbl4CDR1, b_l_precuneus, by = "ID_complete")
F_Tbl4CDR1 <- merge(F_Tbl4CDR1, d_l_precuneus, by = "ID_complete")

##rename columns
F_Tbl4CDR1 <- F_Tbl4CDR1 |> rename(b_l_precuneus = "Name.x", d_l_precuneus = "Name.y")

##Visualization with violin + box again for posterior cingulate
ggplot(data = F_Tbl4CDR1, mapping = aes(x=CDR, y=d_l_precuneus, fill=CDR)) + 
  geom_violin() + 
  geom_boxplot() + scale_fill_manual(values= c("#8cc5e3", "#ea801c")) +
  labs(title = "Predicting Change in L.H. Precuneus Thickness Based on CDR Score", x = "CDR Score", y = "Change in L.H. Precuneus Thickness", fill = "Cognitive Impairment")

########################
####   R Isthmus    ####
########################

##create table w/ baseline entorhinal, change in.., and CDR score
odd_ix <- seq(1, 322, 2)
even_ix <- seq(2, 322, 2)
d_r_isthmus <- Tbl4CDR1$rh_isthmuscingulate_thickness[odd_ix] - Tbl4CDR1$rh_isthmuscingulate_thickness[even_ix]
b_r_isthmus <- Tbl4CDR1$rh_isthmuscingulate_thickness[odd_ix]

##turning these tables into data frames
b_r_isthmus <- data.frame(Name=b_r_isthmus)
d_r_isthmus <- data.frame(Name=d_r_isthmus)

##adding ID_complete column
b_r_isthmus$ID_complete <- Tbl4CDR1$ID_complete[odd_ix]
d_r_isthmus$ID_complete <- Tbl4CDR1$ID_complete[odd_ix]

##merge tables into final table for CDR1 analysis
F_Tbl4CDR1 <- merge(F_Tbl4CDR1, b_r_isthmus, by = "ID_complete")
F_Tbl4CDR1 <- merge(F_Tbl4CDR1, d_r_isthmus, by = "ID_complete")

##rename columns
F_Tbl4CDR1 <- F_Tbl4CDR1 |> rename(b_r_isthmus = "Name.x", d_r_isthmus = "Name.y")

##reorder
F_Tbl4CDR1 <- F_Tbl4CDR1 |> relocate(b_l_entorhinal:d_l_entorhinal, .after = d_r_entorhinal)

##visualization with box and violin 
ggplot(data = F_Tbl4CDR1, mapping = aes(x=CDR, y=d_r_isthmus, fill=CDR)) + 
  geom_violin() + 
  geom_boxplot() + scale_fill_manual(values= c("#8cc5e3", "#ea801c")) 
  #labs(title = "Predicting Change in L.H. Entorhinal Thickness Based on CDR Score", x = "CDR Score", y = "Change in L.H. Entorhinal Thickness", fill = "Cognitive Impairment")

########################
####   L Isthmus    ####
########################

##create table w/ baseline entorhinal, change in.., and CDR score
odd_ix <- seq(1, 322, 2)
even_ix <- seq(2, 322, 2)
d_l_isthmus <- Tbl4CDR1$lh_isthmuscingulate_thickness[odd_ix] - Tbl4CDR1$lh_isthmuscingulate_thickness[even_ix]
b_l_isthmus <- Tbl4CDR1$lh_isthmuscingulate_thickness[odd_ix]

##turning these tables into data frames
b_l_isthmus <- data.frame(Name=b_l_isthmus)
d_l_isthmus <- data.frame(Name=d_l_isthmus)

##adding ID_complete column
b_l_isthmus$ID_complete <- Tbl4CDR1$ID_complete[odd_ix]
d_l_isthmus$ID_complete <- Tbl4CDR1$ID_complete[odd_ix]

##merge tables into final table for CDR1 analysis
F_Tbl4CDR1 <- merge(F_Tbl4CDR1, b_l_isthmus, by = "ID_complete")
F_Tbl4CDR1 <- merge(F_Tbl4CDR1, d_l_isthmus, by = "ID_complete")

##rename columns
F_Tbl4CDR1 <- F_Tbl4CDR1 |> rename(b_l_isthmus = "Name.x", d_l_isthmus = "Name.y")

##reorder
F_Tbl4CDR1 <- F_Tbl4CDR1 |> relocate(b_l_isthmus:d_l_entorhinal, .after = d_r_entorhinal)

##visualization with box and violin 
ggplot(data = F_Tbl4CDR1, mapping = aes(x=CDR, y=d_l_isthmus, fill=CDR)) + 
  geom_violin() + 
  geom_boxplot() + scale_fill_manual(values= c("#8cc5e3", "#ea801c")) 
  #labs(title = "Predicting Change in L.H. Isthmus Cingulate Thickness Based on CDR Score", x = "CDR Score", y = "Change in L.H. Isthmus Cingulate Thickness", fill = "Cognitive Impairment")

#############################
## Significant Data Tables ##
#############################

##Table of t-stats and p-values for the 6 regions
p_value <- c(0.0149, 0.0645, 0.0341, 0.3284, 0.4851, 0.6188)
t_stat <- c(2.4621, 1.8616, 2.1368, 0.9804, 0.6998, 0.4986)
region <- c("entorhinal", "entorhinal", "posterior cingulate", "posterior cingulate", "precuneus", "precuneus")
hemisphere <- c("right", "left", "right", "left", "right", "left")

CDR_Time_Interaction <- data.frame(hemisphere, region, t_stat, p_value)
print(CDR_Time_Interaction)

##Main effects of CDR on non-significant CDR*time interactions
p_value <- c(0.9008, 0.2025, 0.2803)
t_stat <- c(0.1249, -1.2798, -1.0834)
region <- c("posterior cingulate", "precuneus", "precuneus")
hemisphere <- c("left", "right", "left")

CDR_maineffects <- data.frame(hemisphere, region, t_stat, p_value)

########################
## Follow-Up CDR*Time ##
########################

## split sample into groups of Cdr 0 and 0,5/1, then run models without the 
##cdr score variable on the two groups to test main effects of time

CDR_0.5_1 <- Tbl4CDR1[!(Tbl4CDR1$CDRSCORE_01 == 0.0),]

CDR_0 <- Tbl4CDR1[Tbl4CDR1$CDRSCORE_01 == 0.0,]


#######################
### Stats for CDR 0 ###
#######################

baseline <- lme(rh_entorhinal_thickness ~ AGE + GENDER + ED + time_years, random = ~ time_years | ID_complete, data=CDR_0)
summary(baseline)

baseline <- lme(lh_entorhinal_thickness ~ AGE + GENDER + ED + time_years, random = ~ time_years | ID_complete, data=CDR_0)
summary(baseline)

baseline <- lme(rh_posteriorcingulate_thickness ~ AGE + GENDER + ED + time_years, random = ~ time_years | ID_complete, data=CDR_0)
summary(baseline)

#########################
## Stats for CDR 0.5/1 ##
#########################

baseline <- lme(rh_entorhinal_thickness ~ AGE + GENDER + ED + time_years, random = ~ time_years | ID_complete, data=CDR_0.5_1)
summary(baseline)

baseline <- lme(lh_entorhinal_thickness ~ AGE + GENDER + ED + time_years, random = ~ time_years | ID_complete, data=CDR_0.5_1)
summary(baseline)

baseline <- lme(rh_posteriorcingulate_thickness ~ AGE + GENDER + ED + time_years, random = ~ time_years | ID_complete, data=CDR_0.5_1)
summary(baseline)

############################
## Table 4 Strat Analysis ##
############################

hemisphere <- c("right", "left", "right")
region <- c("entorhinal", "entorhinal", "posterior cingulate")
p_val_0 <- c(0.4783, 0.2645, 0.2532)
t_stat_0 <- c(-0.7118, -1.1223, -1.1496)
p_val_1 <- c(0.0418, 0.1474, 0.1795)
t_stat_1 <- c(2.0779, 1.4669, 1.3574)
Strat_Analysis_CDR_0 <- data.frame(hemisphere, region, p_val_0, t_stat_0)
Strat_Analysis_CDR_05_1 <- data.frame(hemisphere, region, p_val_1, t_stat_1)

##############################
## GGSEG 4 Important Regions##
##############################
someData = data.frame(
  region = c("entorhinal", "entorhinal", "posterior cingulate", "posterior cingulate", "precuneus", "precuneus", "isthmus cingulate", "isthmus cingulate"), 
  p = c(0.0149, 0.0645, 0.0341, 0.3284, 0.4851, 0.6188, 0.0574, 0.5488),
  t = c(2.4621, 1.8616, 2.1368, 0.9804, 0.6998, 0.4986, 1.9142, 0.6008),
  hemi = c("right", "left", "right", "left", "right", "left", "right", "left")
)

someData <- someData |> mutate(effectsize = t*2/sqrt(161))

ggseg(.data=someData, mapping=aes(fill=p), atlas = ("aseg"), view = "medial", colour="black") +
  theme_void() +
  scale_fill_gradient(low = "#8cc5e3", high = "#ea801c")

ggseg(.data = someData,
      colour = "black",
      mapping = aes(fill = effectsize), position = "dispersed", na.rm=FALSE, view = "medial") +
  scale_fill_gradient(low = "#8cc5e3", high = "#ea801c") +
  labs(fill="Effect Size") +
  theme_void()

effectsize=tstat*2/sqrt(161)

sum(F_Tbl4CDR1$CDR == "Questionable or Mild Impairment")  


#################################
## Loop 4 Whole Brain Analysis ##
#################################

dep_var_1 <- colnames(MRI_QC_Npsych_Long)[3:36]
dep_var_2 <- colnames(MRI_QC_Npsych_Long)[40:73]
dep_var <- c(dep_var_1, dep_var_2)
print(dep_var)
models <- c()

for (i in seq_along(dep_var)){
  outcome <- dep_var[i]
  formula <- as.formula(paste(outcome, "~ AGE + GENDER + ED + time_years"))
  model <- lme(formula, random = ~ 1 + time_years | ID_complete, data = Tbl4CDR1)
  models[[outcome]] <- summary(model)
}
print(models)

for (outcome in dep_var){
  cat("Summary for model with dependent variable:", outcome, "\n")
  print(models[[outcome]])
  cat("\n\n")
}

significant_pvalues <- list()

# Extract p-values <= 0.05 from each model summary
for (outcome in names(models)) {
  # Extract the tTable, which contains the p-values
  tTable <- models[[outcome]]$tTable
  
  # Identify p-values <= 0.05
  sig_pvalues <- tTable[tTable[, "p-value"] <= 0.05, "p-value"]
  
  # Store in the list with the outcome as the key
  significant_pvalues[[outcome]] <- sig_pvalues
}

# Print the significant p-values for each outcome
print(significant_pvalues)