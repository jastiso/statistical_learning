library(ggplot2)
library(dplyr)
library(coin)
library(lmPerm)
library(car)
library(aplpack)
library(lmerTest)
library(RColorBrewer)
library(wesanderson)
library(nationalparkcolors)
library(ez)
library(plyr)
setwd("/Users/stiso/Documents/Python/graphLearning/old_tasks/mTurk-10-node-breaks")
ext = '1'

df = read.csv('experiment/data/preprocessed/taskdata.csv.gz')
summary(df)


###########
# Check accuracy across subjects
acc = dplyr::filter(df, nTries == 1, stage != 'demo') %>%
  group_by(workerid) %>%
  dplyr::summarise(total_acc = table(correct)["True"]/10)
p<-ggplot(acc, aes(x=total_acc)) + 
  geom_histogram(fill='lightblue', color='white')
p
summary(acc$total_acc)

# only run this is some people have acc less than 80
subj_rm = df$workerid == acc$workerid[acc$total_acc < 80]
df = df[!subj_rm,]


######
# Clean up - remove columns you dont need, only get correct trials, etc

# might need to change this depending on the task
df_clean = subset(df, select = -c(target,query,phase,node,event))
df_clean = na.omit(df_clean)
summary(df_clean)

# remove demo
df_clean = dplyr::filter(df_clean, stage != 'demo')

#make factors
cat_vars = c('keyCode', 'workerid','stage', 'hand','hand_transition','is_crosscluster', 'correct')
for (var in cat_vars){
  df_clean[var] = 
    as.factor(unlist(df_clean[var]))
  
}

# remove RT > 2s
df_clean = filter(df_clean, rt < 2000)
df_clean = filter(df_clean, rt > 50)

# log order, and get continuous trial
cum_trial = df_clean$trial
df_clean$stage_num = rep(0, times = length(cum_trial))
for (t in 1:length(cum_trial)){
  curr_stage = df_clean$stage[t]
  if (curr_stage == "walk_two"){
    cum_trial[t] = df_clean$trial[t] + 250
    df_clean$stage_num[t] = 2
  } else if (curr_stage == "walk_three"){
    cum_trial[t] = df_clean$trial[t] + 250*2
    df_clean$stage_num[t] = 3
  } else if (curr_stage == "walk_four"){
    cum_trial[t] = df_clean$trial[t] + 250*3
    df_clean$stage_num[t] = 4
  } else {
    df_clean$stage_num[t] = 1
  }
}

df_clean$cum_trial = cum_trial  
df_clean$log_cum_trial = log10(df_clean$cum_trial + 1)
df_clean$trial = df_clean$trial + 1

# add finger
#{81:'left pinky',
#  87:'left ring',
#  69:'left middle',
#  82:'left index',
#  86:'left thumb',
#  66:'right thumb',
#  85:'right index',
#  73:'right middle',
#  79:'right ring',
#  80:'right pinky',}

# log order, and get continuous trial
finger = character(length = length(df_clean$keyCode))
for (t in 1:length(finger)){
  curr_key = df_clean$keyCode[t]
  if (curr_key == 81 | curr_key == 80){
    finger[t] = 'pinky'
  } else if (curr_key == 87 | curr_key == 79){
    finger[t] = 'ring'
  } else if (curr_key == 69 | curr_key == 73){
    finger[t] = 'middle'
  } else if (curr_key == 82 | curr_key == 85){
    finger[t] = 'index'
  } else if (curr_key == 86 | curr_key == 66){
    finger[t] = 'thumb'
  } else {
    finger[t] = NA
  }
}
df_clean$finger = as.factor(finger)


p<-ggplot(df_clean, aes(x=rt)) + 
  geom_histogram(fill='pink', color='white')
p
ggsave(paste( '/data/preprocessed/images/rt', ext,'.png', sep = ''))
# make rts log
df_clean$rt_raw = df_clean$rt
df_clean$rt = log10(df_clean$rt)
p<-ggplot(df_clean, aes(x=rt)) + 
  geom_histogram(fill='lightblue', color='white')
p
ggsave(paste( 'experiment/data/preprocessed/images/rt_log', ext,'.png', sep = ''))

summary(df_clean)

# remove incorrect trials
df_correct = dplyr::filter(df_clean, correct == "True" & nTries == 1)
df_correct = subset(df_correct, select = -c(correct, nTries))

#df_correct = filter(df_correct, stage != "walk_one")

df_acc = filter(df_clean, nTries == 1)
df_acc$correct = as.factor(as.integer(df_acc$correct))

df_modular_acc = filter(df_acc, is_lattice == 0)

df_left = filter(df_modular, hand == 'left')
df_right = filter(df_modular, hand == 'right')
summary(df_correct)



################################################
# effect of lag 10
################################################
# helper function for reformatting recency
stop_num = 10
f = function(x) {
  if (x > stop_num) {
    y = stop_num
  } else {
    y = x
  }
  return(y)
}
nuissance_reg = lmer(data=df_correct, rt~scale((cum_trial)) + scale(log10(trial)) + finger + hand + hand_transition + stage_num + 
                       (1 + scale((cum_trial))|workerid))
recency_fact = lapply(df_correct$recency, f)
df_correct$recency_fact = unlist(recency_fact)
df_modular = filter(df_correct, is_lattice == 0)

recency_data = data.frame( rt = resid(nuissance_reg), lag10 = df_correct$lag10, recency = df_correct$recency, graph = as.factor(df_correct$is_lattice), 
                           subj = df_correct$workerid, recency_fact = (unlist(recency_fact)))
p = ggplot(data=recency_data, aes(y=rt, x=recency_fact, color = graph)) + geom_jitter(alpha=0.3) + theme_minimal()
p
ggsave(paste( 'experiment/data/preprocessed/images/lag10_mTurk', ext,'.png', sep = ''))


recency_avg <- recency_data %>%
  group_by(subj, recency_fact,graph) %>%
  dplyr::summarise(mean_rt = mean(rt)) 
p = ggplot(data=recency_avg, aes(y=mean_rt, x=recency_fact, color = graph)) + geom_jitter(width=.3, alpha=0.5, size=3)  + theme_minimal()
p #+ geom_line(aes(group=subj)) #+ theme(legend.position = "none")
ggsave(paste( 'experiment/data/preprocessed/images/recency_avg.png', sep = ''))




#################
# LMER

### Learning
# add *is lattice to trial if you have multiple graphs
stat_learn = lmer(data=df_correct, rt~scale((cum_trial))*is_lattice + stage_num*is_lattice + finger + hand + hand_transition + scale(log(recency_fact)) + 
                                        (1 + scale((cum_trial)) + scale(log(recency_fact))|workerid))
anova(stat_learn)

# save residuals
# for later plotting
df_modular$resid = resid(lmer(data=df_modular, rt~scale((cum_trial)) + scale(log(recency_fact)) + finger + hand + hand_transition + stage_num + log10(recency) + 
                                (1 + scale((cum_trial)) + scale(log(recency_fact))|workerid)))

# for max_ent analyses
max_ent_data = select(df_correct, c('walk_id', 'workerid', 'trial', 'is_lattice'))
max_ent_data$resid = resid(stat_learn)
write.csv(max_ent_data, file = 'data/preprocessed/residuals.csv')


### Modular
# change to df modular if multiple graphs
# most basic
stat_surprisal1 = lmer(data=df_modular, rt~scale((cum_trial))*is_crosscluster + stage_num + finger + hand + hand_transition + 
                        (1 + scale((cum_trial))*is_crosscluster|workerid))
anova(stat_surprisal1)
summary(stat_surprisal1)

# adding lag 10 - is it useful to include this with a random slope?
stat_surprisal2 = lmer(data=df_modular, rt~scale((cum_trial))*is_crosscluster +  stage_num + finger + hand + hand_transition + scale(log10(recency_fact)) + 
                         (1 + scale((cum_trial))*is_crosscluster + scale((recency_fact))|workerid))
anova(stat_surprisal2)
summary(stat_surprisal2)

# chi sq
anova(stat_surprisal1, stat_surprisal2, test="Chisq")



### Accuracy
contrasts(df_modular_acc$hand_transition) <- contr.helmert(2)/2
contrasts(df_modular_acc$hand) <- contr.helmert(2)/2
contrasts(df_modular_acc$correct) <- contr.helmert(2)/2
contrasts(df_modular_acc$finger) <- contr.helmert(5)
stat_acc = glmer(data=df_modular_acc, correct~scale(cum_trial)*is_crosscluster + scale(log10(trial)) + stage_num + finger + hand + 
                   hand_transition + (1|workerid), 
                 family = binomial(link = "logit"))
anova(stat_acc)



##################
# Try it again with down sampling number of subjects
#################

subjs = unique(df_modular$workerid)

for (i in seq(1,10)){
  subset_subj = sample(subjs, 9, replace=FALSE)
  idx = unlist(lapply(df_modular$workerid, function(x) is.element(x, subset_subj)))
  tmp_df = df_correct[idx,]
  subset_stat = lmer(data=df_modular, rt~scale((cum_trial)) + scale(log(recency_fact)) + finger + hand + hand_transition + stage_num + log10(recency) + 
                       (1 + scale((cum_trial)) + scale(log(recency_fact))|workerid))
  print(anova(subset_stat))
}

# same but for no pooling
n.sims=100
ps = rep (NA, n.sims)
for (i in seq(1,n.sims)){
  subset_subj = sample(subjs, 7, replace=FALSE)
  idx = unlist(lapply(df_modular$workerid, function(x) is.element(x, subset_subj)))
  tmp_df = df_modular[idx,]
  subset_stat = lm(data=tmp_df, rt~scale((cum_trial))*is_crosscluster + scale(log10(trial)) + stage_num + finger + hand + hand_transition + scale(lag10) + 
                     workerid)
  ps[i] = anova(subset_stat)$`Pr(>F)`[2] < 0.05
}
mean(ps)


###################
# Some plots that are easier in R

avg_data = df_correct %>%
  group_by(cum_trial, keyCode) %>%
  dplyr::summarise(mean_rt = mean(rt_raw), sd_rt = sd(rt_raw))

plot = ggplot(data=avg_data, aes(x=cum_trial, y=mean_rt, color=keyCode))
plot + geom_line(size=1) + ggtitle('RT over time, by Finger') +
  theme_minimal() + scale_color_manual(values = (brewer.pal(11, "Set3"))) + labs(x = 'Trial', y = 'RT (ms)')
ggsave(paste( 'data/preprocessed/images/rt_mTurk_finger.png', sep = ''))


avg_cluster = df_modular %>%
  group_by(cum_trial, is_crosscluster) %>%
  dplyr::summarise(mean_rt = mean(rt_raw), sd_rt = sd(rt_raw), mean_resid = mean(resid))
avg_cluster$is_crosscluster = factor(avg_cluster$is_crosscluster,levels(avg_cluster$is_crosscluster)[c(2,1)])

plot = ggplot(data=avg_cluster, aes(x=cum_trial, y=mean_rt, color=is_crosscluster))
plot + geom_line(size=1) + ggtitle('RT over time, by Transition') +
  theme_minimal() + labs(x = 'Trial', y = 'RT (ms)')
ggsave(paste( 'experiment/data/preprocessed/images/rt_mTurk_cc', ext,'.png', sep = ''))


avg_graph = df_correct %>%
  group_by(cum_trial, is_lattice) %>%
  dplyr::summarise(mean_rt = mean(rt_raw), sd_rt = sd(rt_raw))

plot = ggplot(data=avg_graph, aes(x=cum_trial, y=mean_rt, color=is_lattice))
plot + geom_line(size=1) + ggtitle('RT over time, by Graph') +
  theme_minimal() + scale_color_manual(values = (wes_palette("Royal2"))) + labs(x = 'Trial', y = 'RT (ms)')
ggsave(paste( 'experiment/data/preprocessed/images/rt_mTurk_graph', ext,'.pdf', sep = ''))

avg_rt = df_correct %>%
  group_by(cum_trial) %>%
  dplyr::summarise(mean_rt = mean(rt_raw), sd_rt = sd(rt_raw))

plot = ggplot(data=avg_rt, aes(x=cum_trial, y=mean_rt))
plot + geom_line(size=1) + ggtitle('RT over time, by Graph') +
  theme_minimal() + labs(x = 'Trial', y = 'RT (ms)') + 
  geom_ribbon(aes(x=cum_trial, y=mean_rt, ymax=mean_rt+sd_rt, ymin = mean_rt-sd_rt), color = 'grey', alpha = 0.2)
ggsave(paste( 'experiment/data/preprocessed/images/rt_mTurk', ext,'.png', sep = ''))

avg_prob = df_correct %>%
  group_by(cum_trial, probability) %>%
  dplyr::summarise(mean_rt = mean(rt_raw), sd_rt = sd(rt_raw))

plot = ggplot(data=avg_prob, aes(x=cum_trial, y=mean_rt, color=as.factor(probability)))
plot + geom_line(size=1) + ggtitle('RT over time, by Graph') +
  theme_minimal() + scale_color_manual(values = (wes_palette("Royal2"))) + labs(x = 'Trial', y = 'RT (ms)')
ggsave(paste( 'experiment/data/preprocessed/images/rt_mTurk_prob', ext,'.png', sep = ''))














#bin
nbin = 40
bin_data = data_frame(trial = tapply(avg_rt$cum_trial, cut(avg_rt$cum_trial, nbin), mean), mean_rt = tapply(avg_rt$mean_rt, cut(avg_rt$cum_trial, nbin), mean))

plot = ggplot(data=bin_data, aes(x=trial, y=mean_rt))
plot + geom_line(size=1) + ggtitle('RT over time') +
  theme_minimal() + labs(x = 'Trial', y = 'RT (ms)')  
ggsave(paste( 'experiment/data/preprocessed/images/rt_mTurk_bin', ext,'.png', sep = ''))


plot = ggplot(data=avg_rt, aes(x=cum_trial, y=mean_rt))
plot + geom_point() + geom_smooth(method='lm') + ggtitle('RT over time') +
  theme_minimal() + labs(x = 'Trial', y = 'RT (ms)')
ggsave(paste( 'data/preprocessed/images/rt_mTurk_scatter', ext,'.png', sep = ''))


## Graph
nbin = 40
bin_data_graph= data_frame(trial = c(tapply(avg_graph$cum_trial, cut(avg_graph$cum_trial, nbin), mean), 
                                     tapply(avg_graph$cum_trial, cut(avg_graph$cum_trial, nbin), mean)),
                           mean_rt = c(tapply(filter(avg_graph, is_lattice == 0)$mean_rt, cut(filter(avg_graph, is_lattice == 0)$cum_trial, nbin), mean), 
                                       tapply(filter(avg_graph, is_lattice == 1)$mean_rt, cut(filter(avg_graph, is_lattice == 1)$cum_trial, nbin), mean)),
                           mean_std = c(tapply(filter(avg_graph, is_lattice == 0)$sd_rt, cut(filter(avg_graph, is_lattice == 0)$cum_trial, nbin), mean), 
                                       tapply(filter(avg_graph, is_lattice == 1)$sd_rt, cut(filter(avg_graph, is_lattice == 1)$cum_trial, nbin), mean)),
                           transition = c(rep("modular", times = length(tapply(avg_graph$cum_trial, cut(avg_graph$cum_trial, nbin), mean))), 
                                          rep("lattice", times = length(tapply(avg_graph$cum_trial, cut(avg_graph$cum_trial, nbin), mean)))))


plot = ggplot(data=bin_data_graph, aes(x=trial, y=mean_rt/1000, color = transition))
plot + geom_line(size=1) + ggtitle('RT over time, by Graph') +
  geom_ribbon(aes(x=trial, y=mean_rt/1000, ymax=mean_rt/1000+mean_std/1000, ymin = mean_rt/1000-mean_std/1000, fill=transition), alpha = 0.2) +
  theme_minimal() + labs(x = 'Trial', y = 'RT (ms)') + scale_color_manual(values = c(rgb(101/255,111/255,147/255), rgb(125/255,138/255,95/255))) +
  scale_fill_manual(values = c(rgb(101/255,111/255,147/255), rgb(125/255,138/255,95/255)))
ggsave(paste( 'experiment/data/preprocessed/images/rt_bin_mturk_graph.pdf', sep = ''))


# cross_cluster
bin_data_cluster = data_frame(trial = c(tapply(avg_cluster$cum_trial, cut(avg_cluster$cum_trial, nbin), mean), 
                                        tapply(avg_cluster$cum_trial, cut(avg_cluster$cum_trial, nbin), mean)),
                              mean_rt = c(tapply(filter(avg_cluster, is_crosscluster == "False")$mean_rt, cut(filter(avg_cluster, is_crosscluster == "False")$cum_trial, nbin), mean), 
                                          tapply(filter(avg_cluster, is_crosscluster == "True")$mean_rt, cut(filter(avg_cluster, is_crosscluster == "True")$cum_trial, nbin), mean)),
                              is_crosscluster = c(rep("False", times = length(tapply(avg_cluster$cum_trial, cut(avg_cluster$cum_trial, nbin), mean))), 
                                                  rep("True", times = length(tapply(avg_cluster$cum_trial, cut(avg_cluster$cum_trial, nbin), mean)))),
                              mean_resid = c(tapply(filter(avg_cluster, is_crosscluster == "False")$mean_resid, cut(filter(avg_cluster, is_crosscluster == "False")$cum_trial, nbin), mean), 
                                             tapply(filter(avg_cluster, is_crosscluster == "True")$mean_resid, cut(filter(avg_cluster, is_crosscluster == "True")$cum_trial, nbin), mean)))


plot = ggplot(data=bin_data_cluster, aes(x=trial, y=mean_rt, color = is_crosscluster))
plot + geom_line(size=1) + ggtitle('RT over time, by transition') + 
  theme_minimal() + labs(x = 'Trial', y = 'RT (ms)') + scale_color_manual(values = c(rgb(215/255,190/255,123/255), rgb(33/255,67/255,104/255))) 
  ggsave(paste( 'experiment/data/preprocessed/images/resid_mTurk_bin_cc', ext,'.pdf', sep = ''))
  

# plot difference over time    
surprisal = data.frame( resid = unlist(bin_data_cluster[bin_data_cluster$"is_crosscluster" == "True",4] - bin_data_cluster[bin_data_cluster$"is_crosscluster" == "False",4]),
                       trials = unlist(bin_data_cluster[bin_data_cluster$"is_crosscluster" == "True",1]))
p = ggplot(data=surprisal, aes(x=trials, y=resid))
p + geom_line() + theme_minimal()
ggsave(paste( 'experiment/data/preprocessed/images/diff_resid', ext,'.png', sep = ''))


# probability
bin_data_prob = data_frame(trial = c(tapply(avg_prob$cum_trial, cut(avg_prob$cum_trial, nbin), mean), 
                                        tapply(avg_prob$cum_trial, cut(avg_prob$cum_trial, nbin), mean),
                                     tapply(avg_prob$cum_trial, cut(avg_prob$cum_trial, nbin), mean)),
                              mean_rt = c(tapply(filter(avg_prob, probability == 0.1)$mean_rt, cut(filter(avg_prob, probability == 0.1)$cum_trial, nbin), mean), 
                                          tapply(filter(avg_prob, probability == 0.25)$mean_rt, cut(filter(avg_prob, probability == 0.25)$cum_trial, nbin), mean),
                                          tapply(filter(avg_prob, probability == 0.7)$mean_rt, cut(filter(avg_prob, probability == 0.7)$cum_trial, nbin), mean)),
                           probability = c(rep("low", times = length(tapply(avg_prob$cum_trial, cut(avg_prob$cum_trial, nbin), mean))), 
                                          rep("med", times = length(tapply(avg_prob$cum_trial, cut(avg_prob$cum_trial, nbin), mean))),
                                           rep("high", times = length(tapply(avg_prob$cum_trial, cut(avg_prob$cum_trial, nbin), mean)))))


plot = ggplot(data=bin_data_prob, aes(x=trial, y=mean_rt, color = probability))
plot + geom_line(size=1) + ggtitle('RT over time, by Graph') +
  theme_minimal() + scale_color_manual(values = (wes_palette("Royal2"))) + labs(x = 'Trial', y = 'RT (ms)') + 
  ggsave(paste( 'experiment/data/preprocessed/images/rt_mTurk_bin_prob', ext,'.png', sep = ''))


df_correct = mutate(df_correct, ho_prob = paste(is_crosscluster, '_', probability, sep = ""))
avg_stats = df_correct %>%
  group_by(workerid,ho_prob) %>%
  dplyr::summarise(mean_rt = mean(rt))

plot = ggplot(data=avg_stats, aes(x=ho_prob, y=mean_rt, fill = ho_prob))
plot + #geom_violin(trim = TRUE, position = position_dodge(0.75)) + 
  geom_boxplot(width = 0.15,
               position = position_dodge(0.75)
  ) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.5, position = position_dodge(0.75)) + ggtitle('RT for different categories') +
  theme_minimal() + scale_fill_manual(values = (wes_palette("Chevalier1"))) + labs(x = 'Stat Category', y = 'RT (ms)') + 
  ggsave(paste( 'experiment/data/preprocessed/images/rt_mTurk_stat_cat', ext,'.png', sep = ''))


