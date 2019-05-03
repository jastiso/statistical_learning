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
library(afex)
setwd("/Users/stiso/Documents/Python/graphLearning/mTurk/")

df = read_csv('data/preprocessed/taskdata.csv.gz')
summary(df)


###########
# Check accuracy across subjects
acc = filter(df, nTries == 1, stage != 'demo') %>%
  group_by(workerid) %>%
  dplyr::summarise(total_acc = sum(correct)/10)
p<-ggplot(acc, aes(x=total_acc)) + 
  geom_histogram(fill='lightblue', color='white')
p
summary(acc$total_acc)

######
# Clean up - remove columns you dont need, only get correct trials, etc

df_clean = subset(df, select = -c(assignmentid,uniqueid,walk_id,finger_mapping,target,response,query,phase,node,event,X1))
df_clean = na.omit(df_clean)
summary(df_clean)

# remove demo
df_clean = dplyr::filter(df_clean, stage != 'demo')

#make factors
cat_vars = c('keyCode', 'correct','workerid','stage', 'is_lattice', 'hand','hand_transition','is_crosscluster', 'correct')
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

summary(df_clean)

# remove incorrect trials
df_correct = dplyr::filter(df_clean, correct == TRUE & nTries == 1)
df_correct = subset(df_correct, select = -c(correct, nTries))

#df_correct = filter(df_correct, stage != "walk_one")

df_acc = filter(df_clean, nTries == 1)
df_acc$correct = as.factor(as.integer(df_acc$correct))

df_modular = filter(df_correct, is_lattice == 0)
df_modular_acc = filter(df_acc, is_lattice == 0)

df_left = filter(df_modular, hand == 'left')
df_right = filter(df_modular, hand == 'right')
summary(df_correct)

#################
# LMER


stat_learn = lmer(data=df_correct, rt~scale(cum_trial) * is_lattice + scale(log10(trial)) + finger + hand + hand_transition + stage_num + scale(lag10) + 
                    scale(recency) + (1 + scale(cum_trial) + scale(lag10) + scale(recency)|workerid))
anova(stat_learn)

contrasts(df_correct$is_lattice) <- contr.helmert(2)/2
stat_graph = lmer(data=df_correct, rt~scale(cum_trial)*is_lattice + scale(log10(trial)) + finger + hand + hand_transition + stage_num + scale(lag10) + 
                    scale(recency) + (1 + scale(cum_trial)*is_lattice|workerid))
anova(stat_graph)

stat_surprisal = lmer(data=df_modular, rt~scale(cum_trial)*is_crosscluster + scale(log10(trial)) + stage_num + finger + hand + hand_transition + scale(lag10) + 
                        scale(recency) + (1 + scale(cum_trial)*is_crosscluster + scale(lag10) + scale(recency)|workerid))
anova(stat_surprisal)

contrasts(df_modular_acc$hand_transition) <- contr.helmert(2)/2
contrasts(df_modular_acc$hand) <- contr.helmert(2)/2
contrasts(df_modular_acc$correct) <- contr.helmert(2)/2
contrasts(df_modular_acc$finger) <- contr.helmert(5)
stat_acc = glmer(data=df_modular_acc, correct~scale(cum_trial)*is_crosscluster + scale(log10(trial)) + stage_num + finger + hand + 
                   hand_transition + (1|workerid), 
                 family = binomial(link = "logit"))
anova(stat_acc)

contrasts(df_correct$hand_transition) <- contr.helmert(2)/2
contrasts(df_correct$hand) <- contr.helmert(2)/2
df_correct$communicability = as.factor(df_correct$communicability)
contrasts(df_correct$communicability) <- contr.helmert(5)
contrasts(df_correct$finger) <- contr.helmert(5)
stat_comm = lmer(data=df_correct, rt~scale(cum_trial)*(communicability) + scale(log10(trial)) + stage_num + finger + hand + hand_transition + scale(lag10) + 
                   scale(recency) + (1 + scale(cum_trial)*(communicability)|workerid))
anova(stat_comm)

###################
# Some plots that are easier in R

avg_data = df_correct %>%
  group_by(cum_trial, keyCode) %>%
  dplyr::summarise(mean_rt = mean(rt), sd_rt = sd(rt))

plot = ggplot(data=avg_data, aes(x=cut(cum_trial), y=mean_rt, color=keyCode))
plot + geom_line(size=1) + ggtitle('RT over time, by Finger') +
  theme_minimal() + scale_color_manual(values = (brewer.pal(10, "Set3"))) + labs(x = 'Trial', y = 'RT (ms)')
ggsave(paste( 'data/preprocessed/images/rt_mTurk_finger.png', sep = ''))


avg_cluster = df_modular %>%
  group_by(cum_trial, is_crosscluster) %>%
  dplyr::summarise(mean_rt = mean(rt), sd_rt = sd(rt))
avg_cluster$is_crosscluster = factor(avg_cluster$is_crosscluster,levels(avg_cluster$is_crosscluster)[c(2,1)])

plot = ggplot(data=avg_cluster, aes(x=cum_trial, y=mean_rt, color=is_crosscluster))
plot + geom_line(size=1) + ggtitle('RT over time, by Transition') +
  theme_minimal() + labs(x = 'Trial', y = 'RT (ms)')
ggsave(paste( 'data/preprocessed/images/rt_mTurk_cc.png', sep = ''))


avg_graph = df_correct %>%
  group_by(cum_trial, is_lattice) %>%
  dplyr::summarise(mean_rt = mean(rt), sd_rt = sd(rt))

plot = ggplot(data=avg_graph, aes(x=cum_trial, y=mean_rt, color=is_lattice))
plot + geom_line(size=1) + ggtitle('RT over time, by Graph') +
  theme_minimal() + scale_color_manual(values = (wes_palette("Royal2"))) + labs(x = 'Trial', y = 'RT (ms)')
ggsave(paste( 'data/preprocessed/images/rt_mTurk_graph.png', sep = ''))

avg_rt = df_correct %>%
  group_by(cum_trial) %>%
  dplyr::summarise(mean_rt = mean(rt), sd_rt = sd(rt))

#bin
breaks = seq(from=1, to=1000, by=20)
bin_data = data_frame(trial = stats.bin(breaks, avg_rt$cum_trial, N=50), mean_rt = stats.bin(avg_rt$cum_trial, avg_rt$mean_rt, N=50, breaks))

plot = ggplot(data=avg_rt, aes(x=cum_trial, y=mean_rt))
plot + geom_line(size=1) + ggtitle('RT over time, by Graph') +
  theme_minimal() + labs(x = 'Trial', y = 'RT (ms)') + 
  geom_ribbon(aes(x=cum_trial, y=mean_rt, ymax=mean_rt+sd_rt, ymin = mean_rt-sd_rt), color = 'grey', alpha = 0.2)
ggsave(paste( 'data/preprocessed/images/rt_mTurk.png', sep = ''))

plot = ggplot(data=bin_data, aes(x=trial, y=mean_rt))
plot + geom_line(size=1) + ggtitle('RT over time, by Graph') +
  theme_minimal() + labs(x = 'Trial', y = 'RT (ms)') + 
ggsave(paste( 'data/preprocessed/images/rt_mTurk_bin.png', sep = ''))




plot = ggplot(data=avg_rt, aes(x=cum_trial, y=mean_rt))
plot + geom_point() + geom_smooth(method='lm') + ggtitle('RT over time, by Graph') +
  theme_minimal() + labs(x = 'Trial', y = 'RT (ms)')
ggsave(paste( 'data/preprocessed/images/rt_mTurk_scatter.png', sep = ''))

