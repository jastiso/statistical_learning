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

# load raw data (output of python script)
df = read.csv('behavior_preprocessed/group_behavior.csv')
demo = read.csv('behavioral_data_raw/demo.csv')
df = merge(df, demo, by='subj')
summary(df)


######
# Clean up - remove columns you dont need, etc

df_clean = subset(df, select = -c(ISI_raw,onset_raw,path,resp_raw,pID))
summary(df_clean)

# chage typing response into linear 1-2-3-4 instead of a,b,c,d
df_clean$typing_raw = as.numeric(mapvalues(df_clean$typing_raw, from = c('a', 'b', 'c', 'd'), to = c('1','2','3','4')))


# add "reward" column for subjects that had points feedback
df_clean$points = TRUE
df_clean$points[df_clean$subj %in% c(1,3,8)] = FALSE

# check accuracy
acc = dplyr::filter(df_clean) %>%
  group_by(subj) %>%
  dplyr::summarise(total_acc = mean(correct_raw))
mean(acc$total_acc)
sd(acc$total_acc)
tmp = (dplyr::filter(df_clean) %>%
         group_by(subj) %>%
         dplyr::summarise(med_rt = median(rt_raw)))
summary(tmp$med_rt)
sd(tmp$med_rt)
#make factors
cat_vars = c('graph', 'correct_raw', 'hand','subj','hand_transition','transition', 'walk', 'points')
for (var in cat_vars){
  df_clean[var] = 
    as.factor(unlist(df_clean[var]))
  
}

# log transfor RT (looks much more normal after)
p<-ggplot(df_clean, aes(x=rt_raw)) + 
  geom_histogram(fill='pink', color='white')
p
ggsave(paste( 'behavior_preprocessed/images/rt_ecog.png', sep = ''))
# make rts log
df_clean$rt = log10(df_clean$rt_raw)
p<-ggplot(df_clean, aes(x=rt)) + 
  geom_histogram(fill='lightblue', color='white')
p
ggsave(paste( 'behavior_preprocessed/images/rt_log_ecog.png', sep = ''))

summary(df_clean)

# add finger
finger = character(length = length(df_clean$resp))
for (t in 1:length(finger)){
  curr_key = df_clean$resp[t]
  if (curr_key == 'q' | curr_key == 'p'){
    finger[t] = 'pinky'
  } else if (curr_key == 'w' | curr_key == 'o'){
    finger[t] = 'ring'
  } else if (curr_key == 'e' | curr_key == 'i'){
    finger[t] = 'middle'
  } else if (curr_key == 'r' | curr_key == 'u'){
    finger[t] = 'index'
  } else if (curr_key == 'v' | curr_key == 'b'){
    finger[t] = 'thumb'
  } else {
    finger[t] = NA
  }
}
df_clean$finger = as.factor(finger)

# remove incorrect trials (but save for accuracy models)
df_correct = dplyr::filter(df_clean, correct_raw == 1)
df_correct = dplyr::filter(df_correct, rt_raw > 0.05)
df_acc = dplyr::filter(df_clean, rt_raw > 0.05)
df_acc$correct_raw = as.factor(df_acc$correct_raw)

# get recency as a factor
stop_num = 10
f = function(x) {
  if (x > stop_num) {
    y = stop_num
  } else {
    y = x
  }
  return(y)
}
recency_fact = lapply(df_correct$recency, f)
df_correct$recency_fact = unlist(recency_fact)
df_acc$recency_fact = unlist(lapply(df_acc$recency, f))
# data frmae of only modular graphs
df_modular = dplyr::filter(df_correct, graph == "modular")

save(df_correct, file = 'behavior_preprocessed/clean.RData')



#################
# LMER

# learn and graph
stat_learn = lmer(data=df_correct, rt~scale((order))*graph + sex + yob + finger + hand + typing_raw + hand_transition + scale(block)*graph + points + scale(log(recency_fact)) + scale(sess) + (1 + scale((order))*graph + scale(log(recency_fact)) |subj))
anova(stat_learn)

  # save residuals
df_correct$resid = resid(stat_learn)
write.csv(df_correct, file = 'behavior_preprocessed/residuals.csv')


### surprisal
contrasts(df_modular$hand_transition) <- contr.helmert(2)/2
contrasts(df_modular$hand) <- contr.helmert(2)/2
contrasts(df_modular$correct) <- contr.helmert(2)/2
contrasts(df_modular$finger) <- contr.helmert(5)
df_modular$order_sq = (scale((df_modular$order))^2)
stat_surprisal1 = lmer(data=df_modular, rt~scale((order))*transition + sex + scale(yob) + typing_raw + finger + points + hand + hand_transition +  scale(block) + scale(log(recency_fact)) + sess + 
                         (1 + scale((order))*transition |subj))


stat_surprisal2 = lmer(data=df_modular, rt~scale((order))*transition + sex + yob + finger + hand + typing_raw + hand_transition + scale(block) + points + scale(log(recency_fact)) + scale(sess) + (1 + scale((order))*transition + scale(log(recency_fact)) |subj))
# chi sq
anova(stat_surprisal2, stat_surprisal1, test="Chisq")

anova(stat_surprisal2)
summary(stat_surprisal2)


# accuracy 
contrasts(df_acc$hand_transition) <- contr.helmert(2)/2
contrasts(df_acc$hand) <- contr.helmert(2)/2
contrasts(df_acc$correct) <- contr.helmert(2)/2
contrasts(df_acc$finger) <- contr.helmert(5)
df_acc$order_sq = (scale((df_acc$order))^2)
stat_acc1 = glmer(data=df_acc, correct_raw~order_sq*graph + finger + hand_transition + scale(block)*graph + scale(log(recency_fact)) + scale(sess) +
                  (1 + order_sq |subj),
                family = binomial(link = "logit"))
stat_acc2 = glmer(data=df_acc, correct_raw~scale((order))*graph + finger + hand_transition + scale(block)*graph + scale(log(recency_fact)) + scale(sess) +
                    (1 + scale((order)) |subj),
                  family = binomial(link = "logit"))
anova(stat_acc1, stat_acc2, test="Chisq")

summary(stat_acc1)
summary(stat_acc2)



##################
# Plot

plot = ggplot(data = df_correct, aes(x=rt, fill=as.factor(typing_raw)))
plot + geom_histogram(aes(y=..density..),alpha=.8)

# average over participants for each trial, grouped by graph
avg_data = df_correct %>%
  group_by(order, graph) %>%
  dplyr::summarise(mean_rt = mean(rt_raw), sd_rt = sd(rt_raw)/sqrt(length(rt_raw)))
# same for acc
avg_acc = df_acc %>%
  group_by(order, graph) %>%
  dplyr::summarise(mean_acc = mean(as.numeric(correct_raw)-1), sd_rt = sd(as.numeric(correct_raw)-1)/sqrt(length(correct_raw)))

plot = ggplot(data=avg_data, aes(x=order, y=mean_rt))
plot + geom_line(size=1) + ggtitle('RT over time, by Graph') +
  theme_minimal() + labs(x = 'Trial', y = 'RT (ms)') + 
  geom_ribbon(aes(x=order, y=mean_rt, ymax=mean_rt+sd_rt, ymin = mean_rt-sd_rt), color = 'grey', alpha = 0.2)
ggsave(paste( 'behavior_preprocessed/images/rt_ECoG.png', sep = ''))

# averages grouped by CC
avg_cluster = filter(df_modular, sess == '1') %>%
  group_by(order, transition) %>%
  dplyr::summarise(mean_rt = mean(rt_raw), sd_rt = sd(rt_raw))

plot = ggplot(data=avg_cluster, aes(x=order, y=mean_rt, color=transition))
plot + geom_line(size=1) + ggtitle('RT over time, by Transition') +
  theme_minimal() + labs(x = 'Trial', y = 'RT (ms)')
ggsave(paste( 'behavior_preprocessed/images/rt_ECoG_cc.png', sep = ''))

# bin data to make the plots look a little smoother (grouped by cluster)
nbin = 40
bin_data_cluster = data_frame(trial = c(tapply(avg_cluster$order, cut(avg_cluster$order, nbin), mean), 
                                        tapply(avg_cluster$order, cut(avg_cluster$order, nbin), mean)),
                              mean_rt = c(tapply(filter(avg_cluster, transition == "False")$mean_rt, cut(filter(avg_cluster, transition == "False")$order, nbin), mean), 
                                          tapply(filter(avg_cluster, transition == "True")$mean_rt, cut(filter(avg_cluster, transition == "True")$order, nbin), mean)),
                              transition = c(rep("False", times = length(tapply(avg_cluster$order, cut(avg_cluster$order, nbin), mean))), 
                                                  rep("True", times = length(tapply(avg_cluster$order, cut(avg_cluster$order, nbin), mean)))))


plot = ggplot(data=bin_data_cluster, aes(x=trial, y=mean_rt, color = transition))
plot + geom_line(size=1) + ggtitle('RT over time, by Edge Type') +
  theme_minimal() + labs(x = 'Trial', y = 'RT (ms)') + scale_color_manual(values = c(rgb(215/255,190/255,123/255), rgb(33/255,67/255,104/255))) +
  ggsave(paste( 'behavior_preprocessed/images/rt_ECoG_bin_cc.pdf', sep = ''))

## bin data that was grouped by Graph
nbin = 40
bin_data_graph= data_frame(trial = c(tapply(avg_data$order, cut(avg_data$order, nbin), mean), 
                                        tapply(avg_data$order, cut(avg_data$order, nbin), mean)),
                              mean_rt = c(tapply(filter(avg_data, graph == "modular")$mean_rt, cut(filter(avg_data, graph == "modular")$order, nbin), mean), 
                                          tapply(filter(avg_data, graph == "lattice")$mean_rt, cut(filter(avg_data, graph == "lattice")$order, nbin), mean)),
                              std_rt = c(tapply(filter(avg_data, graph == "modular")$sd_rt, cut(filter(avg_data, graph == "modular")$order, nbin), mean), 
                                       tapply(filter(avg_data, graph == "lattice")$sd_rt, cut(filter(avg_data, graph == "lattice")$order, nbin), mean)),
                              transition = c(rep("modular", times = length(tapply(avg_data$order, cut(avg_data$order, nbin), mean))), 
                                             rep("lattice", times = length(tapply(avg_data$order, cut(avg_data$order, nbin), mean)))))


plot = ggplot(data=bin_data_graph, aes(x=trial, y=mean_rt, color = transition))
plot + geom_line(size=1) + ggtitle('RT over time, by Graph') +
  geom_ribbon(aes(x=trial, y=mean_rt, ymax=mean_rt+std_rt, ymin = mean_rt-std_rt, fill=transition), alpha = 0.2) +
  theme_minimal() + labs(x = 'Trial', y = 'RT (ms)') + scale_color_manual(values = c(rgb(101/255,111/255,147/255), rgb(125/255,138/255,95/255))) +
  scale_fill_manual(values = c(rgb(101/255,111/255,147/255), rgb(125/255,138/255,95/255)))+
  ggsave(paste( 'behavior_preprocessed/images/rt_ECoG_bin_graph.pdf', sep = ''))

## bin data grouped by graph for acc
nbin = 40
bin_data_graph_acc= data_frame(trial = c(tapply(avg_acc$order, cut(avg_acc$order, nbin), mean), 
                                     tapply(avg_acc$order, cut(avg_acc$order, nbin), mean)),
                           mean_rt = c(tapply(filter(avg_acc, graph == "modular")$mean_acc, cut(filter(avg_acc, graph == "modular")$order, nbin), mean), 
                                       tapply(filter(avg_acc, graph == "lattice")$mean_acc, cut(filter(avg_acc, graph == "lattice")$order, nbin), mean)),
                           std_rt = c(tapply(filter(avg_acc, graph == "modular")$sd_rt, cut(filter(avg_acc, graph == "modular")$order, nbin), mean), 
                                      tapply(filter(avg_acc, graph == "lattice")$sd_rt, cut(filter(avg_acc, graph == "lattice")$order, nbin), mean)),
                           transition = c(rep("modular", times = length(tapply(avg_acc$order, cut(avg_acc$order, nbin), mean))), 
                                          rep("lattice", times = length(tapply(avg_acc$order, cut(avg_acc$order, nbin), mean)))))


plot = ggplot(data=bin_data_graph_acc, aes(x=trial, y=mean_rt, color = transition))
plot + geom_line(size=1) + ggtitle('Accuracy over time, by Graph') +
  geom_ribbon(aes(x=trial, y=mean_rt, ymax=mean_rt+std_rt, ymin = mean_rt-std_rt, fill=transition), alpha = 0.2) +
  theme_minimal() + labs(x = 'Trial', y = 'Accuracy (%)') + scale_color_manual(values = c(rgb(101/255,111/255,147/255), rgb(125/255,138/255,95/255))) +
  scale_fill_manual(values = c(rgb(101/255,111/255,147/255), rgb(125/255,138/255,95/255)))+
  ggsave(paste( 'behavior_preprocessed/images/rt_ECoG_bin_graph_acc.pdf', sep = ''))


