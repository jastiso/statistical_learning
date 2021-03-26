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
  theme_minimal() + scale_color_manual(values=c(rgb(101/255,111/255,147/255), rgb(174/255,116/255,133/255)))
p


############################################# Ahat static beta
b_df = read.csv('ephys_analysis/ahat_seq.csv')
b_df = merge(b_df, demo, by='subj')
b_df$subj = as.factor(b_df$subj)
b_df$beta_rank = as.factor(rank(b_df$beta))
b_df$is_lat = as.factor(strtoi(b_df$subj) %% 2)

summary(b_df)

p = ggplot(data=b_df, aes(x=trial, y=(dist), group = subj, color=as.factor(beta))) + 
  geom_line() + scale_color_manual(values=colorRampPalette(brewer.pal(9,'OrRd'))(10)) +
  theme_minimal() 
p
ggsave("ephys_img/ahat_seq.pdf",p)

p = ggplot(data=b_df, aes(x=trial, y=(dist), group = subj, color=is_lat)) + 
  geom_line() + scale_color_manual(values=c(rgb(174/255,116/255,133/255), rgb(101/255,111/255,147/255)))+
  theme_minimal() 
p
ggsave("ephys_img/ahat_seq.pdf",p)

b_df_avg = dplyr::data_frame(subj = character(0), yob = numeric(), sex = character(0), beta = numeric(),
                      block = numeric(), mean_dist = numeric(), is_lat = numeric())
shift = 100
win = 500
nBlock = (1000 - win)/shift + 1
st = 1
block_inds = list()
for (i in seq(1,nBlock)){
  block_inds = c(block_inds, list(c(st, st+win-1)))
  st = st + shift
}

for (s in unique(b_df$subj)){
  blocked = list()
  l = nrow(b_df_avg) + 1
  b_df_avg[l:(l+nBlock-1),'subj'] = rep(s, ntimes=nBlock)
  b_df_avg$yob[l:(l+nBlock-1)] = b_df[b_df$subj==s,'yob'][1]
  b_df_avg$sex[l:(l+nBlock-1)] =  as.character(b_df[b_df$subj==s,'sex'][1])
  b_df_avg$beta[l:(l+nBlock-1)] =  b_df[b_df$subj==s,'beta'][1]
  b_df_avg$is_lat[l:(l+nBlock-1)] =  as.numeric(b_df[b_df$subj==s,'is_lat'][1]) - 1
  b_df_avg$block[l:(l+nBlock-1)] = seq(1,nBlock)
  curr_dist = b_df[b_df$subj == s, 'dist']
  for (b in seq(1,nBlock)){
    blocked = c(blocked, mean(curr_dist[seq(block_inds[[b]][1],block_inds[[b]][2])]))
  }
  b_df_avg$mean_dist[l:(l+nBlock-1)] = unlist(blocked)
}
b_df_avg$subj = as.factor(b_df_avg$subj)
b_df_avg$is_lat = as.factor(b_df_avg$is_lat)
p = ggplot(data=b_df_avg, aes(x=block, y=mean_dist, group = subj, color=is_lat)) + 
  geom_line() +scale_color_manual(values=c(rgb(174/255,116/255,133/255), rgb(101/255,111/255,147/255))) +
  theme_minimal() 
p

############################################# Ahat veriable beta
bv_df = read.csv('ephys_analysis/ahat_block.csv')
bv_df = merge(bv_df, demo, by='subj')
bv_df$subj = as.factor(bv_df$subj)
bv_df$beta_rank = as.factor(rank(bv_df$beta))
bv_df$is_lat = as.factor(strtoi(bv_df$subj) %% 2)
bv_df_ext = bv_df[bv_df$beta_block >= 1000 | bv_df$beta_block <= 0,]
bv_df = bv_df[bv_df$beta_block < 1000 & bv_df$beta_block > 0,]
bv_df$beta_diff = log10(bv_df$beta) - log10(bv_df$beta_block)
bv_df = bv_df[bv_df$beta < 1000 & bv_df$beta > 0,]
bv_df$beta_diff_log = log10(abs(bv_df$beta_diff))
summary(bv_df)

pd = position_dodge(0.05)
p = ggplot(data=bv_df, aes(x=block, y=log10(beta_block), group = subj, color=as.factor(beta_rank))) + 
  geom_line(position=pd) + geom_point(size=3, position=pd) + 
  scale_color_manual(values=colorRampPalette(brewer.pal(9,'OrRd'))(10)) +
  theme_minimal()
p
ggsave("ephys_img/ahat_block.eps",p)
stat = lmer(data=bv_df, log10(beta_block)~scale(block) + (1|subj))
anova(stat)
summary(stat)

p = ggplot(data=bv_df_ext, aes(x=block)) + 
  geom_bar(aes(fill=as.factor(beta_block)))  + #scale_fill_manual(values=c(rgb(125/255,138/255,95/255), rgb(101/255,111/255,147/255))) +
  theme_minimal() 
p
ggsave("ephys_img/ahat_block_ext.pdf",p)

################################################## Repeat for mturk
bv_df_mturk= read.csv('ephys_analysis/ahat_block_mturk.csv')
bv_df_mturk$subj = as.factor(bv_df_mturk$subj)
bv_df_mturk$beta_rank = as.factor(rank(bv_df_mturk$beta))
bv_df_mturk_ext = bv_df_mturk[bv_df_mturk$beta_block >= 1000 | bv_df_mturk$beta_block <= 0,]
bv_df_mturk = bv_df_mturk[bv_df_mturk$beta_block < 1000 & bv_df_mturk$beta_block > 0,]
bv_df_mturk$beta_diff = log10(bv_df_mturk$beta) - log10(bv_df_mturk$beta_block)
bv_df_mturk = bv_df_mturk[bv_df_mturk$beta < 1000 & bv_df_mturk$beta > 0,]
bv_df_mturk$beta_diff_log = log10(abs(bv_df_mturk$beta_diff))
summary(bv_df_mturk)

pd = position_dodge(0.05)
p = ggplot(data=bv_df_mturk, aes(x=block, y=log10(beta_block), group = subj, color=as.factor(beta_rank))) + 
  geom_line(position=pd) + geom_point(size=3, position=pd) + 
  scale_color_manual(values=colorRampPalette(brewer.pal(9,'OrRd'))(25)) +
  theme_minimal()
p
ggsave("ephys_img/ahat_block_mturk.eps",p)
stat = lmer(data=bv_df_mturk,beta_diff_log~(block)+ (1|subj))
anova(stat)
summary(stat)

p = ggplot(data=bv_df_mturk_ext, aes(x=block)) + 
  geom_bar(aes(fill=as.factor(beta_block)))  + #scale_fill_manual(values=c(rgb(125/255,138/255,95/255), rgb(101/255,111/255,147/255))) +
  theme_minimal() 
p
ggsave("ephys_img/ahat_block_mturk_ext.pdf",p)


################################################## Combined

b_df_avg$sex = as.factor(b_df_avg$sex)
b_df_avg = b_df_avg[b_df_avg$subj != '18',]
latent_df = dplyr::filter(df, space=='latent')
data_all = merge(b_df_avg, latent_df, by=c('subj','sex','yob','block','is_lat'))
summary(data_all)
data_all$beta_rank = as.factor(rank(data_all$beta))
df_avg = dplyr::summarise(group_by(data_all, subj,block,beta_rank,beta), mean_corr = mean(corr),
                          full_corr = mean(full_corr), sd_corr = sd(corr)/sqrt(length(corr)), sex = sex[1], yob = yob[1])


pd = position_dodge(0.05)
p = ggplot(data=filter(df_avg, block == nBlock), aes(x=log10(beta), y=mean_corr, color=beta_rank, group=subj)) + 
   geom_point(size=3) + 
  theme_minimal() + scale_color_manual(values=brewer.pal(9,'OrRd'))
p

pd = position_dodge(0.05)
tmp = filter(df_avg, block == 1)
stat = lm(data=tmp, mean_corr ~ full_corr + sex + yob)
tmp$resid = resid(stat)
p = ggplot(data=filter(tmp, beta != 1000), aes(x=log10(beta), y=resid)) + 
  geom_point(size=3, position=pd) + geom_smooth( method='lm') + 
  theme_minimal() + scale_color_manual(values=brewer.pal(9,'OrRd'))
p

stat = lm(data=filter(tmp, beta != 1000), mean_corr~log10(beta) + full_corr + sex + yob)
anova(stat)
summary(stat)


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
rs = list()
for (i in unique(df_null$set)){
  tmp = dplyr::filter(df_null, is_lat == 1 & set == i)
  rs = c(rs, cor(tmp$module_dist, log10(tmp$beta), method='spearman'))
}
max(unlist(rs))
stat = cor.test(df$module_dist, log10(df$beta), method='spearman')

p = ggplot(data = df, aes(x=log10(beta), y=(module_dist))) + 
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





