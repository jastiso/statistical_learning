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
setwd("/Users/stiso/Documents/Code/graph_learning/ECoG_data/")

s = 10
load('behavior_preprocessed/clean.RData')
data = readMat(paste('/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_analysis/subj',s,'/peak_power.mat',sep=""))
max_ent = readMat(paste('/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_analysis/subj',s,'/max_ent_contrast.mat',sep=""))

# format 
df_subj = dplyr::filter(df_correct, subj == s)
elecs = unlist(data$elec.labels)
regions = unlist(data$region)

# get only elecs in grey matter
elecs = elecs[!grepl("White Matter", regions)]
regions = regions[!grepl("White Matter", regions)]

df_pow <- data.frame(matrix(ncol = length(elecs) + 3, nrow = length(t(data$good.trial.idx))))
nam <- c("order", "mod", "max_ent", elecs)
colnames(df_pow) <- nam
df_pow$order = t(data$good.trial.idx)
df_pow$mod = data$module.idx

for (e in 1:length(elecs)){
  elec <- elecs[e]
  curr = data$peak.power[e,]
  df_pow[elec] = curr

}
# combine
df_pow = merge(df_subj,df_pow, by.x = 'order')
drops = c('X','Unnamed..0','graph','walk','typing_raw','subj','correct_raw','resp','sess')
df_fit = df_pow[,!(names(df_pow) %in% drops)]


##### Stats #################################
#save all the models
models = list()
ps = list()
betas = list()
for (e in elecs){
  formula = paste(e, '~ transition + order + finger + hand + hand_transition + log(recency_fact) + block + rt')
  fit = lm(data=df_fit, formula)
  anova(fit)
  ps[[e]] = (anova(fit)$`Pr(>F)`[1])
  betas[[e]] = summary(lm.beta(fit))$coefficients[2,2] # standard beta
  models[[e]] = fit
}
p.adjust(ps, n=length(elecs), method="fdr")
sum(p.adjust(ps, n=length(elecs), method="fdr") < 0.05)
betas[p.adjust(ps, n=length(elecs), method="fdr") < 0.05]
print(paste("Electrode with statistically significant module contrast: ", unlist(elecs[p.adjust(ps, n=length(elecs), method="fdr") < 0.05]), sep=""))
# save file
mod_stats = data.frame(elecs = elecs, p = p.adjust(ps, n=length(elecs), method="fdr"), betas = unlist(betas), 
                       region = regions[unlist(lapply(regions, function(x) x != "No_label"))])
write.csv(mod_stats, paste('/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_analysis/subj',s,'/mod_stats.csv', sep=""))

# repeat for ramping contrast
models_ramp = list()
ps_ramp = list()
betas_ramp = list()
for (e in elecs){
  formula = paste(e, '~ mod + order + finger + hand + hand_transition + log(recency_fact) + block + rt')
  fit = lm(data=df_fit, formula)
  anova(fit)
  ps_ramp[[e]] = (anova(fit)$`Pr(>F)`[1])
  betas_ramp[[e]] = summary(lm.beta(fit))$coefficients[2,2] # standard beta
  models_ramp[[e]] = fit
}
p.adjust(ps_ramp, n=length(elecs), method="fdr")
sum(p.adjust(ps_ramp, n=length(elecs), method="fdr") < 0.05)
betas_ramp[p.adjust(ps_ramp, n=length(elecs), method="fdr") < 0.05]
print(paste("Electrode with statistically significant ramping contrast: ", unlist(elecs[p.adjust(ps_ramp, n=length(elecs), method="fdr") < 0.05]), sep=""))
# save file
ramp_stats = data.frame(elecs = elecs, p = p.adjust(ps_ramp, n=length(elecs), method="fdr"), betas = unlist(betas_ramp), 
                        region = regions[unlist(lapply(regions, function(x) x != "No_label"))])
write.csv(ramp_stats, paste('/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_analysis/subj',s,'/ramp_stats.csv', sep=""))



# repeat for max-ent contrast
models_max_ent = list()
ps_max_ent = list()
betas_max_ent = list()
for (e in elecs){
  formula = paste(e, '~ mod + order + finger + hand + hand_transition')
  fit = lm(data=df_fit, formula)
  anova(fit)
  ps_max_ent[[e]] = (anova(fit)$`Pr(>F)`[1])
  betas_max_ent[[e]] = summary(lm.beta(fit))$coefficients[2,2] # standard beta
  models_max_ent[[e]] = fit
}
p.adjust(ps_max_ent, n=length(elecs), method="fdr")
sum(p.adjust(ps_max_ent, n=length(elecs), method="fdr") < 0.05)
betas_max_ent[p.adjust(ps_max_ent, n=length(elecs), method="fdr") < 0.05]
print(paste("Electrode with statistically significant max_ent contrast: ", unlist(elecs[p.adjust(ps_max_ent, n=length(elecs), method="fdr") < 0.05]), sep=""))
# save file
max_ent_stats = data.frame(elecs = elecs, p = p.adjust(ps_max_ent, n=length(elecs), method="fdr"), betas = unlist(betas_max_ent), 
                        region = regions[unlist(lapply(regions, function(x) x != "No_label"))])
write.csv(max_ent_stats, paste('/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_analysis/subj',s,'/max_ent_stats.csv', sep=""))

