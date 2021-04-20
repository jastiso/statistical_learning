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
library(reshape)
setwd("/Users/stiso/Documents/Code/graph_learning/ECoG_data/")

############################################################### Ephys
# load data on time varying RSA
df = read.csv('ephys_analysis/block_searchlight.csv')
# demographics (sex and age)
demo = read.csv('behavioral_data_raw/demo.csv')
demo$subj = as.factor(demo$subj)
df = merge(df, demo, by='subj')
# load behavioral data to get betas
b_df = read.csv('ephys_analysis/ahat_seq.csv')
b_df = b_df %>%
  group_by(subj) %>%
  dplyr::summarise(beta = mean(beta))
df = merge(df, b_df, by='subj')
# formatting stuff
df$beta_rank = as.factor(rank(df$beta))
df$subj = as.factor(df$subj)
df$is_lat = as.factor(strtoi(df$subj) %% 2)
df$norm_corr = df$corr/df$full_corr
#There shouldn't be empy fields, and sex, beta, and yob should be here'
summary(df)

## plots! 
df$uniqueid = paste(df$subj, df$space, sep = "_")
# average over contacts
df_avg = dplyr::summarise(group_by(df, subj,space,block, uniqueid), is_lat = is_lat[1], mean_corr = mean(corr), sd_corr = sd(corr)/sqrt(length(corr)),
                          mean_norm_corr = mean(norm_corr), sd_norm_corr = sd(norm_corr)/sqrt(length(norm_corr)),
                          mean_full_corr = mean(full_corr), sex = sex[1], yob = yob[1], beta = beta[1], beta_rank = beta_rank[1])
pd = position_dodge(0.05)
# color by space
p = ggplot(data=df_avg, aes(x=block, y=mean_norm_corr, color=space, group=uniqueid)) + #geom_boxplot(aes(x=block, y=mean_corr, group=block), position=pd) +
  geom_errorbar(aes(ymin=mean_norm_corr-sd_norm_corr, ymax = mean_norm_corr+sd_norm_corr), color='black', width=0.01, position=pd) +
  geom_line(aes(color=space),position=pd) + geom_point(aes(color=space),size=3, position=pd) + 
  theme_minimal() + scale_color_manual(values=c('grey','gold'))
p
ggsave("ephys_img/block_rsa.eps",p)

# the same, but colored by graph
pd = position_dodge(0.05)
p = ggplot(data=filter(df_avg,space=='latent'), aes(x=block, y=mean_norm_corr, color=is_lat, group=uniqueid)) + #geom_boxplot(aes(x=block, y=mean_corr, group=block), position=pd) +
  geom_errorbar(aes(ymin=mean_norm_corr-sd_norm_corr, ymax = mean_norm_corr+sd_norm_corr), color='black', width=0.01, position=pd) +
  geom_line(aes(color=is_lat),position=pd) + geom_point(aes(color=is_lat),size=3, position=pd) + 
  theme_minimal() + scale_color_manual(values=c(rgb(174/255,116/255,133/255),rgb(101/255,111/255,147/255)))
p
# colored by graph values with only lines
pd = position_dodge(0.05)
p = ggplot(data=filter(dplyr::filter(df_avg, space=='latent')), aes(x=block, y=mean_norm_corr, color=is_lat, group=is_lat)) + 
  geom_smooth(method='lm', alpha=0.1) + 
  theme_minimal() + scale_color_manual(values=c(rgb(174/255,116/255,133/255), rgb(101/255,111/255,147/255))) 
p
ggsave("ephys_img/block_rsa_graph.svg",p)

# now, colored by beta
pd = position_dodge(0.05)
p = ggplot(data=filter(df_avg, space=='latent'), aes(x=block, y=mean_norm_corr, color=beta_rank, group=subj)) + #geom_boxplot(aes(x=block, y=mean_corr, group=block), position=pd) +
  geom_errorbar(aes(ymin=mean_norm_corr-sd_norm_corr, ymax = mean_norm_corr+sd_norm_corr), color='black', width=0.01, position=pd) +
  geom_line(position=pd) + geom_point(size=3, position=pd) + 
  theme_minimal() + scale_color_manual(values=colorRampPalette(brewer.pal(9,'OrRd'))(10))
p
# colored by beta values with only lines
pd = position_dodge(0.05)
p = ggplot(data=filter(df_avg, space=='latent', beta < 1000), aes(x=block, y=mean_norm_corr, color=beta_rank, group=subj, fill=beta_rank)) + 
  geom_smooth(method='lm', alpha=0.1) + 
  theme_minimal() + scale_color_manual(values=colorRampPalette(brewer.pal(9,'OrRd'))(9)) +
  scale_fill_manual(values=colorRampPalette(brewer.pal(9,'OrRd'))(9))
p
ggsave("ephys_img/block_rsa_beta.svg",p)

## Stats
# mixed effect model - are corrs different between vis and latent spaces? do cors change over time?
stat = lmer(data=df_avg,mean_norm_corr~block*space + sex + yob + (1|subj))
anova(stat)
summary(stat)
# ttest for the first block
stat = t.test(filter(df_avg, space=='latent', block==1)$mean_corr,
              filter(df_avg, space=='euclid', block==1)$mean_corr, paired=TRUE)
stat
# what about just for the latent space? Does the graph type matter?
stat_latent = lmer(data=dplyr::filter(df_avg, space=='latent'),mean_norm_corr~block*is_lat + sex + scale(yob) + (1|subj))
anova(stat_latent)
summary(stat_latent)
# what about just for the latent space? Does the graph type matter?
stat_latent = lmer(data=dplyr::filter(df_avg, space=='latent'),mean_norm_corr~block*is_lat + sex + scale(yob) + (1|subj))
anova(stat_latent)
summary(stat_latent)
# do slopes change with betas?
stat_latent = lmer(data=dplyr::filter(df_avg, space=='latent',beta <= 1000),mean_norm_corr~scale(block)*scale(log10(beta)) + sex + yob + (1|subj))
anova(stat_latent)
summary(stat_latent)

############################################# Ahat static beta
# load simulated A_hat(t)
b_df = read.csv('ephys_analysis/ahat_seq.csv')
# since we already loaded demos, merge it here
b_df = merge(b_df, demo, by='subj')
b_df$subj = as.factor(b_df$subj)
b_df$beta_rank = as.factor(rank(b_df$beta))
b_df$is_lat = as.factor(strtoi(b_df$subj) %% 2)
# shouldnt be empty variables, and should have demographics
summary(b_df)

# plot, colored by beta
p = ggplot(data=b_df, aes(x=trial, y=(dist), group = subj, color=as.factor(beta))) + 
  geom_line() + scale_color_manual(values=colorRampPalette(brewer.pal(9,'OrRd'))(10)) +
  theme_minimal() 
p
ggsave("ephys_img/ahat_seq.pdf",p)

# same plot, colored by graph
p = ggplot(data=b_df, aes(x=trial, y=(dist), group = subj, color=is_lat)) + 
  geom_line() + scale_color_manual(values=c(rgb(174/255,116/255,133/255), rgb(101/255,111/255,147/255)))+
  theme_minimal() 
p
ggsave("ephys_img/ahat_seq_graph.pdf",p)

# to make this comparable to ephys data, we want to group the similarity into blocks
b_df_avg = dplyr::data_frame(subj = character(0), yob = numeric(), sex = character(0), beta = numeric(),
                      block = numeric(), mean_dist = numeric(), is_lat = numeric())
# these should match the ephys script
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

# repeat the graph plot, but in blocks
p = ggplot(data=b_df_avg, aes(x=block, y=mean_dist, group = subj, color=is_lat)) + 
  geom_line() +scale_color_manual(values=c(rgb(174/255,116/255,133/255), rgb(101/255,111/255,147/255))) +
  theme_minimal() 
p

############################################# Ahat veriable beta
# load data from Ahat being calculated in blocks
bv_df = read.csv('ephys_analysis/ahat_block.csv')
bv_df = merge(bv_df, demo, by='subj')
bv_df$subj = as.factor(bv_df$subj)
bv_df$beta_rank = as.factor(rank(bv_df$beta))
bv_df$is_lat = as.factor(strtoi(bv_df$subj) %% 2)
# this needs extra processing to separate "extreme" beta values
bv_df_ext = bv_df[bv_df$beta_block >= 1000 | bv_df$beta_block <= 0,]
bv_df = bv_df[bv_df$beta_block < 1000 & bv_df$beta_block > 0,]
bv_df$beta_diff = log10(bv_df$beta) - log10(bv_df$beta_block)
bv_df = bv_df[bv_df$beta < 1000 & bv_df$beta > 0,]
bv_df$beta_diff_log = log10(abs(bv_df$beta_diff))
# beta values should have max < 1000, and min > 0
summary(bv_df)

# plot - color by original beta value, and only show blocked betas between 0 and 1000
pd = position_dodge(0.05)
p = ggplot(data=bv_df, aes(x=block, y=log10(beta_block), group = subj, color=as.factor(beta_rank))) + 
  geom_line(position=pd) + geom_point(size=3, position=pd) + 
  scale_color_manual(values=colorRampPalette(brewer.pal(9,'OrRd'))(10)) +
  theme_minimal()
p
ggsave("ephys_img/ahat_block.eps",p)

# same plot, colored by graph type
p = ggplot(data=bv_df, aes(x=block, y=log10(beta_block), group = subj, color=as.factor(is_lat))) + 
  geom_line(position=pd) + geom_point(size=3, position=pd) + 
  scale_color_manual(values=c(rgb(174/255,116/255,133/255), rgb(101/255,111/255,147/255))) +
  theme_minimal()
p
ggsave("ephys_img/ahat_block.eps",p)


# mixed effects model - any group differences over time?
stat = lmer(data=bv_df, log10(beta_block)~scale(block) + (1|subj))
anova(stat)
summary(stat)

# plot the extreme values
p = ggplot(data=bv_df_ext, aes(x=block)) + 
  geom_bar(aes(fill=as.factor(beta_block)))  + #scale_fill_manual(values=c(rgb(125/255,138/255,95/255), rgb(101/255,111/255,147/255))) +
  theme_minimal() 
p
ggsave("ephys_img/ahat_block_ext.pdf",p)

################################################## Repeat for mturk
# same data (blocked A_hat) but for the mturk data
bv_df_mturk= read.csv('ephys_analysis/ahat_block_mturk.csv')
bv_df_mturk$subj = as.factor(bv_df_mturk$subj)
bv_df_mturk$beta_rank = as.factor(rank(bv_df_mturk$beta))
# extra processing for separating extreme values
bv_df_mturk_ext = bv_df_mturk[bv_df_mturk$beta_block >= 1000 | bv_df_mturk$beta_block <= 0,]
bv_df_mturk = bv_df_mturk[bv_df_mturk$beta_block < 1000 & bv_df_mturk$beta_block > 0,]
bv_df_mturk$beta_diff = log10(bv_df_mturk$beta) - log10(bv_df_mturk$beta_block)
bv_df_mturk = bv_df_mturk[bv_df_mturk$beta < 1000 & bv_df_mturk$beta > 0,]
bv_df_mturk$beta_diff_log = log10(abs(bv_df_mturk$beta_diff))
# beta values should all be between 0 and 1000
summary(bv_df_mturk)

# plot - colored by original beta, with extreme values removed
pd = position_dodge(0.05)
p = ggplot(data=bv_df_mturk, aes(x=block, y=log10(beta_block), group = subj, color=as.factor(beta_rank))) + 
  geom_line(position=pd) + geom_point(size=3, position=pd) + 
  scale_color_manual(values=colorRampPalette(brewer.pal(9,'OrRd'))(25)) +
  theme_minimal()
p
ggsave("ephys_img/ahat_block_mturk.eps",p)

# mixed effects model - any group level changes?
stat = lmer(data=bv_df_mturk,beta_diff_log~(block)+ (1|subj))
anova(stat)
summary(stat)

# extreme values
p = ggplot(data=bv_df_mturk_ext, aes(x=block)) + 
  geom_bar(aes(fill=as.factor(beta_block)))  + #scale_fill_manual(values=c(rgb(125/255,138/255,95/255), rgb(101/255,111/255,147/255))) +
  theme_minimal() 
p
ggsave("ephys_img/ahat_block_mturk_ext.pdf",p)


########################################################## Linear Disc
# test linear discriminability both behavioral and neural data
df = read.csv('ephys_analysis/mod_dist.csv')
demo = read.csv('behavioral_data_raw/demo.csv')
loss = read.csv('ephys_analysis/loss_data.csv')
demo$subj = as.factor(demo$subj)
demo$is_lat = as.factor(strtoi(demo$subj) %% 2)
df = merge(df, demo, by='subj')
loss = merge(loss,demo, by='subj') # the proportion of misclassified nodes
loss = loss[loss$is_lat==0,]
loss = melt(loss, id=c('subj','yob','sex','beta','is_lat'))
summary(loss)

# plot loss between neural and behavioral
pd = position_dodge(0.7)
p = ggplot(data = dplyr::filter(loss, variable != 'losses', variable != 'loss_beh'), aes(x=variable, y=value, color=variable, group = subj)) + 
  theme_minimal() + geom_line(position=pd, color='black') +
   geom_point(size=5, position=pd) + scale_color_manual(values=c('darkgrey','darkgrey'))
p
ggsave("ephys_img/loss_neur.eps",p)

# paired ttest between neural and behavioral
stat = t.test(loss[loss$variable=='loss_pca','value'],loss[loss$variable=='losses_bpca','value'],paried=TRUE)
stat






