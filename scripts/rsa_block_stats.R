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

############################################################### Ephys
df = read.csv('ephys_analysis/block_searchlight.csv')
demo = read.csv('behavioral_data_raw/demo.csv')
demo$subj = as.factor(demo$subj)
df = merge(df, demo, by='subj')
df = df[df$subj != 6,]
df$subj = as.factor(df$subj)
summary(df)

stat = lmer(data=df,corr~block*space + sex + yob + (1|subj))
anova(stat)
summary(stat)

stat_latent = lmer(data=filter(df, space=='latent'),corr~block + sex + yob + (1|subj))
anova(stat_latent)
summary(stat_latent)

stat = t.test(filter(df,space=='latent',block==1)$corr,filter(df,space=='latent', block==2)$corr, paired=TRUE)
stat

df$uniqueid = paste(df$subj, df$space, sep = "_")
df_avg = dplyr::summarise(group_by(df, subj,space,block, uniqueid), mean_corr = mean(corr), sd_corr = sd(corr)/sqrt(length(corr)))
pd = position_dodge(0.05)
p = ggplot(data=df_avg, aes(x=block, y=mean_corr, color=space, group=uniqueid)) + #geom_boxplot(aes(x=block, y=mean_corr, group=block), position=pd) +
  geom_errorbar(aes(ymin=mean_corr-sd_corr, ymax = mean_corr+sd_corr), color='black', width=0.01, position=pd) +
  geom_line(aes(color=space),position=pd) + geom_point(aes(color=space),size=3, position=pd) + 
  theme_minimal() + scale_color_manual(values=c('grey','gold'))
p
ggsave("ephys_img/block_rsa.pdf",p)


############################################# Ahat
b_df = read.csv('ephys_analysis/ahat_seq.csv')
b_df = merge(b_df, demo, by='subj')
b_df$subj = as.factor(b_df$subj)
b_df$beta_rank = as.factor(rank(b_df$beta))
summary(b_df)

p = ggplot(data=b_df, aes(x=trial, y=(dist), group = subj, color=as.factor(beta))) + 
  geom_line() + scale_color_manual(values=brewer.pal(9,'OrRd')) +
  theme_minimal() 
p
ggsave("ephys_img/ahat_seq.pdf",p)

b_df$block = floor(b_df$trial/501)
b_df_avg = dplyr::summarise(group_by(b_df, subj, beta, sex, yob, block), mean_dist = mean(dist), sd_dist = sd(dist))

p = ggplot(data=dplyr::filter(b_df_avg, block < 2), aes(x=block, y=mean_dist, group = subj, color=as.factor(beta))) + 
  geom_line() + scale_color_manual(values=brewer.pal(9,'OrRd')) +
  theme_minimal() 
p

################################################## Combined

data_all = merge(b_df,df, by=c('subj','sex','yob'))
summary(data_all)

df_avg = dplyr::summarise(group_by(data_all, subj,space,block, uniqueid,beta_rank,beta), mean_corr = mean(corr), sd_corr = sd(corr)/sqrt(length(corr)))

pd = position_dodge(0.05)
p = ggplot(data=filter(df_avg, space=='latent'), aes(x=block, y=mean_corr, color=beta_rank, group=uniqueid)) + 
  geom_errorbar(aes(ymin=mean_corr-sd_corr, ymax = mean_corr+sd_corr), color='black', width=0.01, position=pd) +
  geom_line(position=pd) + geom_point(size=3, position=pd) + 
  theme_minimal() + scale_color_manual(values=brewer.pal(9,'OrRd'))
p

pd = position_dodge(0.05)
p = ggplot(data=filter(data_all,block==1, space=='latent'), aes(x=corr, y=log10(beta), color=as.factor((beta)), group=subj)) + 
  #geom_errorbar(aes(ymin=mean_corr-sd_corr, ymax = mean_corr+sd_corr), color='black', width=0.01, position=pd) +
  geom_line(position=pd) + geom_point(size=3, position=pd) + 
  theme_minimal() + scale_color_manual(values=brewer.pal(9,'OrRd'))
p

b_df_avg$block = b_df_avg$block + 1
df_avg = merge(df_avg, b_df_avg, by = c('subj','block'))
stat = lmer(data=filter(df_avg, space == 'latent'), mean_corr~mean_dist + sex + yob + (1|subj))
anova(stat)
summary(stat)


p = ggplot(data=dplyr::filter(df_avg, space == 'latent'), aes(x=mean_corr, y=(mean_dist), color=as.factor(block), group=subj)) + 
  #geom_errorbar(aes(ymin=mean_corr-sd_corr, ymax = mean_corr+sd_corr), color='black', width=0.01, position=pd) +
  geom_point(size=3) + 
  theme_minimal() #+ scale_color_manual(values=brewer.pal(2,'OrRd'))
p

##########################################################
df = read.csv('ephys_analysis/mod_dist.csv')
demo = read.csv('behavioral_data_raw/demo.csv')
demo$subj = as.factor(demo$subj)
df = merge(df, demo, by='subj')
summary(df)

stat = lm(data=df, module_dist~log10(beta))
summary(stat)