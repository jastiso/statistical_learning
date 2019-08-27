# instantaneous amplitude
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
ia = readMat(paste('/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_analysis/subj',s,'/inst_amp.mat',sep=""))

# format 
df_subj = dplyr::filter(df_correct, subj == s)
elecs = unlist(data$elec.labels)
regions = unlist(data$region)

# get only elecs in grey matter
elecs = elecs[unlist(lapply(regions, function(x) x != "No_label"))]

df_ia <- data.frame(matrix(ncol = length(elecs) + 2, nrow = length(t(ia$good.trial.idx))))
nam <- c("order", "mod", elecs)
colnames(df_ia) <- nam
df_ia$order = t(data$good.trial.idx)
df_ia$mod = data$module.idx

for (e in 1:length(elecs)){
  elec <- elecs[e]
  
  curr_ia = ia$ia.mat[e,]
  df_ia[elec] = curr_ia
}

df_ia = merge(df_subj,df_ia, by.x = 'order')
drops = c('X','Unnamed..0','graph','walk','typing_raw','subj','correct_raw','resp','sess')
df_ia_fit = df_ia[,!(names(df_ia) %in% drops)]

##### Stats for Inst Amp #################################
#save all the models
models_ia = list()
ps_ia = list()
betas_ia = list()
for (e in elecs){
  formula = paste(e, '~ transition*order + finger + hand + hand_transition')
  fit = lm(data=df_ia_fit, formula)
  anova(fit)
  ps_ia[[e]] = (anova(fit)$`Pr(>F)`[1])
  betas_ia[[e]] = summary(lm.beta(fit))$coefficients[2,2] # standard beta
  models_ia[[e]] = fit
}
p.adjust(ps_ia, n=length(elecs), method="fdr")
sum(p.adjust(ps_ia, n=length(elecs), method="fdr") < 0.05)
betas_ia[p.adjust(ps_ia, n=length(elecs), method="fdr") < 0.05]
print(paste("Electrode with statistically significant module contrast: ", unlist(elecs[p.adjust(ps_ia, n=length(elecs), method="fdr") < 0.05]), sep=""))
# save file
mod_stats_ia = data.frame(elecs = elecs, p = p.adjust(ps_ia, n=length(elecs), method="fdr"), betas = unlist(betas_ia), 
                          region = regions[unlist(lapply(regions, function(x) x != "No_label"))])
write.csv(mod_stats_ia, paste('/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_analysis/subj',s,'/mod_stats_ia.csv', sep=""))

# repeat for ramping contrast
models_ramp_ia = list()
ps_ramp_ia = list()
betas_ramp_ia = list()
for (e in elecs){
  formula = paste(e, '~ mod*order + finger + hand + hand_transition')
  fit = lm(data=df_ia_fit, formula)
  anova(fit)
  ps_ramp_ia[[e]] = (anova(fit)$`Pr(>F)`[1])
  betas_ramp_ia[[e]] = summary(lm.beta(fit))$coefficients[2,2] # standard beta
  models_ramp_ia[[e]] = fit
}
p.adjust(ps_ramp_ia, n=length(elecs), method="fdr")
sum(p.adjust(ps_ramp_ia, n=length(elecs), method="fdr") < 0.05)
betas_ramp_ia[p.adjust(ps_ramp_ia, n=length(elecs), method="fdr") < 0.05]
print(paste("Electrode with statistically significant ramping contrast: ", unlist(elecs[p.adjust(ps_ramp_ia, n=length(elecs), method="fdr") < 0.05]), sep=""))
# save file
ramp_stats_ia = data.frame(elecs = elecs, p = p.adjust(ps_ramp_ia, n=length(elecs), method="fdr"), betas = unlist(betas_ramp_ia), 
                           region = regions[unlist(lapply(regions, function(x) x != "No_label"))])
write.csv(ramp_stats_ia, paste('/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_analysis/subj',s,'/ramp_stats_ia.csv', sep=""))

