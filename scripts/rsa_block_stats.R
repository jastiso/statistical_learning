library(ggplot2)
library(dplyr)
library(coin)
library(lmPerm)
library(car)
library(aplpack)
library(lmerTest)
library(RColorBrewer)
library(wesanderson)
library(ez)
library(plyr)
library(lm.beta)
setwd("/Users/stiso/Documents/Code/graph_learning/ECoG_data/")

df = read.csv('ephys_analysis/block_searchlight.csv')
demo = read.csv('behavioral_data_raw/demo.csv')
demo$subj = as.factor(demo$subj)
df = merge(df, demo, by='subj')
df$subj = as.factor(df$subj)
summary(df)

stat = lmer(data=df,corr~block*space + sex + yob + (1|subj))
anova(stat)
summary(stat)

stat_latent = lmer(data=filter(df, space=='latent'),corr~block + sex + yob + (1|subj))
anova(stat_latent)
summary(stat_latent)

df$uniqueid = paste(df$subj, df$space, sep = "_")
df_avg = dplyr::summarise(group_by(df, subj,space,block, uniqueid), mean_corr = mean(corr), sd_corr = sd(corr)/sqrt(length(corr)))
pd = position_dodge(0.05)
ggplot(data=df_avg, aes(x=block, y=mean_corr, color=space, group=uniqueid)) + 
  geom_errorbar(aes(ymin=mean_corr-sd_corr, ymax = mean_corr+sd_corr), color='black', width=0.01, position=pd) +
  geom_line(aes(color=space),position=pd) + geom_point(aes(color=space),size=3, position=pd) + 
  theme_minimal() + scale_color_manual(values=c('grey','gold'))


b_df = read.csv('ephys_analysis/ahat_block.csv')
b_df = merge(b_df, demo, by='subj')
b_df$subj = as.factor(b_df$subj)
summary(b_df)

b_df$beta_diff = abs(b_df$beta_block - b_df$beta)
stat = lmer(data=b_df, scale(dist) ~ block + sex + yob + (1|subj))
anova(stat)
summary(stat)

stat = t.test(filter(b_df, block == 1)$dist, filter(b_df, block == 2)$dist, paired=TRUE)
stat


pd = position_dodge(0.05)
ggplot(data=b_df, aes(x=block, y=log10(dist), group=subj)) + 
  geom_line(position=pd) + geom_point(size=3, position=pd) + 
  theme_minimal() 

