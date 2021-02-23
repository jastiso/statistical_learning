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
df$subj = as.factor(df$subj)
df$is_lat = as.factor(strtoi(df$subj) %% 2)
summary(df)

stat = lmer(data=df,corr~block*space + sex + yob + (1|subj))
anova(stat)
summary(stat)

stat_latent = lmer(data=filter(df, space=='latent'),corr~block*is_lat + sex + yob + (1|subj))
anova(stat_latent)
summary(stat_latent)

stat = t.test(filter(df,space=='latent',block==1)$corr,filter(df,space=='latent', block==2)$corr, paired=TRUE)
stat

df$uniqueid = paste(df$subj, df$space, sep = "_")
df_avg = dplyr::summarise(group_by(df, subj,space,block, uniqueid), is_lat = is_lat[1], mean_corr = mean(corr), sd_corr = sd(corr)/sqrt(length(corr)))
pd = position_dodge(0.05)
p = ggplot(data=df_avg, aes(x=block, y=mean_corr, color=space, group=uniqueid)) + #geom_boxplot(aes(x=block, y=mean_corr, group=block), position=pd) +
  geom_errorbar(aes(ymin=mean_corr-sd_corr, ymax = mean_corr+sd_corr), color='black', width=0.01, position=pd) +
  geom_line(aes(color=space),position=pd) + geom_point(aes(color=space),size=3, position=pd) + 
  theme_minimal() + scale_color_manual(values=c('grey','gold'))
p
ggsave("ephys_img/block_rsa.eps",p)

pd = position_dodge(0.05)
p = ggplot(data=filter(df_avg,space=='latent'), aes(x=block, y=mean_corr, color=is_lat, group=uniqueid)) + #geom_boxplot(aes(x=block, y=mean_corr, group=block), position=pd) +
  geom_errorbar(aes(ymin=mean_corr-sd_corr, ymax = mean_corr+sd_corr), color='black', width=0.01, position=pd) +
  geom_line(aes(color=is_lat),position=pd) + geom_point(aes(color=is_lat),size=3, position=pd) + 
  theme_minimal() + scale_color_manual(values=c(rgb(101/255,111/255,147/255), rgb(125/255,138/255,95/255)))
p
ggsave("ephys_img/block_rsa.eps",p)


############################################# Ahat static beta
b_df = read.csv('ephys_analysis/ahat_seq.csv')
b_df = merge(b_df, demo, by='subj')
b_df$subj = as.factor(b_df$subj)
b_df$beta_rank = as.factor(rank(b_df$beta))
summary(b_df)

p = ggplot(data=b_df, aes(x=trial, y=(dist), group = subj, color=as.factor(beta))) + 
  geom_line() + scale_color_manual(values=colorRampPalette(brewer.pal(9,'OrRd'))(10)) +
  theme_minimal() 
p
ggsave("ephys_img/ahat_seq.pdf",p)

b_df$block = floor(b_df$trial/501)
b_df_avg = dplyr::summarise(group_by(b_df, subj, beta, sex, yob, block), mean_dist = mean(dist), sd_dist = sd(dist))

p = ggplot(data=dplyr::filter(b_df_avg, block < 2), aes(x=block, y=mean_dist, group = subj, color=as.factor(beta))) + 
  geom_line() +scale_color_manual(values=colorRampPalette(brewer.pal(9,'OrRd'))(10)) +
  theme_minimal() 
p

############################################# Ahat veriable beta
bv_df = read.csv('ephys_analysis/ahat_block.csv')
bv_df = merge(bv_df, demo, by='subj')
bv_df$subj = as.factor(bv_df$subj)
bv_df$beta_rank = as.factor(rank(bv_df$beta))
bv_df$beta_diff = log10(bv_df$beta) - log10(bv_df$beta_block)
#bv_df = bv_df[bv_df$beta < 1000 & bv_df$beta > 0,]
summary(bv_df)

pd = position_dodge(0.05)
p = ggplot(data=bv_df, aes(x=block, y=log10(beta_block), group = subj, color=as.factor(beta_rank))) + 
  geom_line(position=pd) + geom_point(size=3, position=pd) + 
  scale_color_manual(values=colorRampPalette(brewer.pal(9,'OrRd'))(10)) +
  theme_minimal()
p
ggsave("ephys_img/ahat_block.svg",p)
stat = lmp(data=bv_df, (beta_block)~block)
anova(stat)
stat = lmp(data=bv_df[bv_df$beta_diff != Inf,], abs(beta_diff)~block)
anova(stat)

pd = position_dodge(0.05)
p = ggplot(data=bv_df, aes(x=block, y=(dist), group = subj, color=as.factor(beta_rank))) + 
  geom_line(position=pd) + geom_point(size=3, position=pd) + 
  scale_color_manual(values=colorRampPalette(brewer.pal(9,'OrRd'))(10)) +
  theme_minimal() 
p
stat = lmp(data=bv_df, dist~block)
anova(stat)


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
df_null = read.csv('ephys_analysis/null_mod_dist.csv')
comp = read.csv('ephys_analysis/comp_data.csv')
demo = read.csv('behavioral_data_raw/demo.csv')
demo$subj = as.factor(demo$subj)
demo$is_lat = as.factor(strtoi(demo$subj) %% 2)
df = merge(df, demo, by='subj')
df_null = merge(df_null, demo, by='subj')
comp = merge(comp, demo, by='subj')
summary(df)

stat = lm(data=df, module_dist~log10(beta)+sex+yob)
summary(stat)

p = ggplot(data = df, aes(x=log10(beta), y=module_dist)) + 
  geom_smooth(method='lm', color='black') + geom_point(size=5, color=rgb(125/255,138/255,95/255)) + 
  theme_minimal() + 
  geom_line(data = df_null, aes(x=log10(beta), y=module_dist, group=set), stat='smooth', method='lm', 
              color=rgb(101/255,111/255,147/255), alpha=0.1)
p
ggsave("ephys_img/mod_dist_neur.svg",p)

p = ggplot(data = df_null, aes(x=log10(beta), y=module_dist, group=set)) + 
  geom_smooth(data = df_null, aes(x=log10(beta), y=module_dist, group=set), method='lm', color=rgb(101/255,111/255,147/255), se=FALSE) + 
  theme_minimal()
p


stat = lmp(data=comp, log10(compress)~is_lat)
summary(stat)
stat = perm(comp[comp$is_lat == 1,]$comp, comp[comp$is_lat == 0,]$comp)
stat

p = ggplot(data = comp, aes(x=log10(beta), y=(compress), color=is_lat, group = is_lat)) + 
  geom_smooth(method='lm', color='black') + geom_point(size=5) + scale_color_manual(values=c(rgb(101/255,111/255,147/255), rgb(125/255,138/255,95/255))) + 
  theme_minimal()
p





