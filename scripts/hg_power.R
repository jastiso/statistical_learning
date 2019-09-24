library(ggplot2)
library(dplyr)
library(coin)
library(lmPerm)
library(car)
library(aplpack)
library(lmerTest)
library(RColorBrewer)
library(wesanderson)
library(R.matlab)
library(ez)
library(plyr)
library(lm.beta)
setwd("/Users/stiso/Documents/Python/graphLearning/ECoG data/")

s = 4
load('behavior_preprocessed/clean.RData')
data = readMat(paste('/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_analysis/subj',s,'/hg_power.mat',sep=""))
theta = readMat(paste('/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_analysis/subj',s,'/theta_peaks.mat',sep=""))

# format 
df_subj = dplyr::filter(df_correct, subj == s)
elecs = unlist(data$elec.labels)
regions = unlist(data$region)

# get only elecs where there are theta oscillations
elecs_osc = elecs[as.logical(theta$ps)]
elecs_other = elecs[!as.logical(theta$ps)]
regions_osc = regions[as.logical(theta$ps)]
regions_other = regions[!as.logical(theta$ps)]
# get only elecs in grey matter
elecs_osc = elecs_osc[unlist(lapply(regions_osc, function(x) x != "No_label"))]
elecs_other = elecs_other[unlist(lapply(regions_other, function(x) x != "No_label"))]


df_pow <- data.frame(matrix(ncol = length(elecs) + 2, nrow = length(t(data$good.trial.idx))))
nam <- c("order", "mod", elecs)
colnames(df_pow) <- nam
df_pow$order = t(data$good.trial.idx)
df_pow$mod = data$module.idx

for (e in 1:length(elecs)){
  elec <- elecs[e]
  curr = data$hg.power[e,]
  df_pow[elec] = curr
  
}
# combine
df_pow = merge(df_subj,df_pow, by.x = 'order')
drops = c('X','Unnamed..0','graph','walk','typing_raw','subj','correct_raw','resp','sess')
df_fit = df_pow[,!(names(df_pow) %in% drops)]




##################
# H1: Decreased activity at transitions
##################
elec_group = elecs_osc
region_group = regions_osc
ext = '_hg_osc'

##################
# H2: Differential cortical activity
##################
elec_group = elecs_other
region_group = regions_other
ext = ''


##### Stats #################################
#save all the models
models = list()
ps = list()
betas = list()
for (e in elec_group){
  formula = paste(e, '~ transition + order + finger + hand + hand_transition')
  fit = lm(data=df_fit, formula)
  anova(fit)
  ps[[e]] = (anova(fit)$`Pr(>F)`[1])
  betas[[e]] = summary(lm.beta(fit))$coefficients[2,2] # standard beta
  models[[e]] = fit
}
p.adjust(ps, n=length(elec_group), method="fdr")
sum(p.adjust(ps, n=length(elec_group), method="fdr") < 0.05)
betas[p.adjust(ps, n=length(elec_group), method="fdr") < 0.05]
print(paste("Electrode with statistically significant module contrast: ", unlist(elec_group[p.adjust(ps, n=length(elec_group), method="fdr") < 0.05]), sep=""))
# save file
mod_stats = data.frame(elecs = elec_group, p = p.adjust(ps, n=length(elec_group), method="fdr"), betas = unlist(betas), 
                       region = region_group[unlist(lapply(region_group, function(x) x != "No_label"))])
write.csv(mod_stats, paste('/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_analysis/subj',s,'/mod_stats', ext, '.csv', sep=""))

# repeat for ramping contrast
models_ramp = list()
ps_ramp = list()
betas_ramp = list()
for (e in elec_group){
  formula = paste(e, '~ mod + order + finger + hand + hand_transition')
  fit = lm(data=df_fit, formula)
  anova(fit)
  ps_ramp[[e]] = (anova(fit)$`Pr(>F)`[1])
  betas_ramp[[e]] = summary(lm.beta(fit))$coefficients[2,2] # standard beta
  models_ramp[[e]] = fit
}
p.adjust(ps_ramp, n=length(elec_group), method="fdr")
sum(p.adjust(ps_ramp, n=length(elec_group), method="fdr") < 0.05)
betas_ramp[p.adjust(ps_ramp, n=length(elec_group), method="fdr") < 0.05]
print(paste("Electrode with statistically significant ramping contrast: ", unlist(elec_group[p.adjust(ps_ramp, n=length(elec_group), method="fdr") < 0.05]), sep=""))
# save file
ramp_stats = data.frame(elecs = elec_group, p = p.adjust(ps_ramp, n=length(elec_group), method="fdr"), betas = unlist(betas_ramp), 
                        region = region_group[unlist(lapply(region_group, function(x) x != "No_label"))])
write.csv(ramp_stats, paste('/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_analysis/subj',s,'/ramp_stats', ext, '.csv', sep=""))


