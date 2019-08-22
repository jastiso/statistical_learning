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
setwd("/Users/stiso/Documents/Python/graphLearning/ECoG data/")

s = 4
load('behavior_preprocessed/clean.RData')
data = readMat(paste('/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_analysis/subj',s,'/peak_power.mat',sep=""))

# format 
df_subj = dplyr::filter(df_correct, subj == s)
elecs = unlist(data$elec.labels)
df_pow <- data.frame(matrix(ncol = length(elecs) + 2, nrow = length(t(data$good.trial.idx))))
nam <- c("order", "mod", elecs)
colnames(df_pow) <- nam
df_pow$order = t(data$good.trial.idx)
df_pow$mod = data$module.idx
for (e in 1:length(elecs)){
  elec <- elecs[e]
  curr = data$peak.power[e,]
  df_pow[elec] = curr
}
# combine
df = merge(df_subj,df_pow, by.x = 'order')
drops = c('X','Unnamed..0','graph','walk','typing_raw','subj','correct_raw','resp','sess')
df_fit = df[,!(names(df) %in% drops)]


##### Stats #################################
#save all the models
models = list()
ps = list()
for (e in elecs){
  formula = paste(e, '~ transition*order + finger + hand + hand_transition')
  fit = lm(data=df_fit, formula)
  anova(fit)
  ps[[e]] = (anova(fit)$`Pr(>F)`[1])
  models[[e]] = fit
}
p.adjust(ps, n=length(elecs), method="fdr")
sum(p.adjust(ps, n=length(elecs), method="fdr") < 0.05)

# repeat for ramping contrast
models_mod = list()
ps_mod = list()
for (e in elecs){
  formula = paste(e, '~ mod*order + finger + hand + hand_transition')
  fit = lm(data=df_fit, formula)
  anova(fit)
  ps_mod[[e]] = (anova(fit)$`Pr(>F)`[1])
  models_mod[[e]] = fit
}
p.adjust(ps_mod, n=length(elecs), method="fdr")
sum(p.adjust(ps_mod, n=length(elecs), method="fdr") < 0.05)

