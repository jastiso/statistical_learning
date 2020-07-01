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

df = read.csv('behavior_preprocessed/group_behavior.csv')
summary(df)


######
# Clean up - remove columns you dont need, only get correct trials, etc

df_clean = subset(df, select = -c(ISI_raw,onset_raw,path,resp_raw,pID))
summary(df_clean)

# chage typing response into linear 1-2-3-4 instead of a,b,c,d
df_clean$typing_raw = as.numeric(mapvalues(df_clean$typing_raw, from = c('a', 'b', 'c', 'd'), to = c('1','2','3','4')))


# add "reward" column for subjects that had points feedback
df_clean$points = TRUE
df_clean$points[df_clean$subj %in% c(1,3,8)] = FALSE

#make factors
cat_vars = c('graph', 'correct_raw', 'hand','subj','hand_transition','transition', 'walk', 'points')
for (var in cat_vars){
  df_clean[var] = 
    as.factor(unlist(df_clean[var]))
  
}

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

# remove incorrect trials
df_correct = dplyr::filter(df_clean, correct_raw == 1)
df_correct = dplyr::filter(df_correct, rt_raw > 0.05)

# add finger
# log order, and get continuous trial
finger = character(length = length(df_correct$resp))
for (t in 1:length(finger)){
  curr_key = df_correct$resp[t]
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
df_correct$finger = as.factor(finger)

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

df_modular = dplyr::filter(df_correct, graph == "modular")

save(df_correct, file = 'behavior_preprocessed/clean.RData')



#################
# LMER

# test if points makes a difference

## learn and graph
stat_learn = lmer(data=df_correct, rt~scale(log10(order))*graph + finger + typing_raw + hand_transition + block + scale(log(recency_fact)) + sess + (1 + scale(log10(order))*graph + scale(log(recency_fact)) |subj))
anova(stat_learn)

# save residuals
df_correct$resid = resid(stat_learn)
write.csv(df_correct, file = 'behavior_preprocessed/residuals.csv')


## graph
stat_graph = lmer(data=df_correct, rt~scale(log10(order))*graph + finger + typing_raw + hand_transition +  block + scale(log(recency_fact)) + sess + (1 + scale(log10(order))*graph|subj))
anova(stat_graph)


### surprisal
stat_surprisal1 = lmer(data=df_modular, rt~scale(log10(order))*transition + typing_raw + finger + hand + hand_transition +  block + scale(log(recency_fact)) + sess + 
                         (1 + scale(log10(order))*transition |subj))


stat_surprisal2 = lmer(data=df_modular, rt~scale(log10(order))*transition + finger + typing_raw + hand + hand_transition +  block + scale(log(recency_fact)) + sess + 
                         (1 + scale(log10(order))*transition + scale(log(recency_fact)) |subj))
# chi sq
anova(stat_surprisal2, stat_surprisal1, test="Chisq")

anova(stat_surprisal2)
summary(stat_surprisal2)

# permutation test: permute transition index within subject
nSim = 500
effect_sizes
for (i in seq(1,nSim)){
  
}


# no pooling
stat_no_pool = lm.beta(lm(data=df_modular, rt~scale(order)*transition + finger + hand_transition +  block + scale(log(recency_fact)) +sess +  subj))
anova(stat_no_pool)
summary(stat_no_pool)

# full pooling
stat_pool = lm.beta(lm(data=filter(df_modular, subj == '8'), rt~scale(order)*transition + finger + hand_transition +  block + scale(log(recency_fact)) ))
anova(stat_pool)
summary(stat_pool)



##################
# Plot

plot = ggplot(data = df_correct, aes(x=rt, fill=as.factor(typing_raw)))
plot + geom_histogram(aes(y=..density..),alpha=.8)

avg_data = df_correct %>%
  group_by(order, graph) %>%
  dplyr::summarise(mean_rt = mean(rt_raw), sd_rt = sd(rt_raw))

plot = ggplot(data=avg_data, aes(x=order, y=mean_rt))
plot + geom_line(size=1) + ggtitle('RT over time, by Graph') +
  theme_minimal() + labs(x = 'Trial', y = 'RT (ms)') + 
  geom_ribbon(aes(x=order, y=mean_rt, ymax=mean_rt+sd_rt, ymin = mean_rt-sd_rt), color = 'grey', alpha = 0.2)
ggsave(paste( 'behavior_preprocessed/images/rt_ECoG.png', sep = ''))


avg_cluster = filter(df_modular, sess == '1') %>%
  group_by(order, transition) %>%
  dplyr::summarise(mean_rt = mean(rt_raw), sd_rt = sd(rt_raw))

plot = ggplot(data=avg_cluster, aes(x=order, y=mean_rt, color=transition))
plot + geom_line(size=1) + ggtitle('RT over time, by Transition') +
  theme_minimal() + labs(x = 'Trial', y = 'RT (ms)')
ggsave(paste( 'behavior_preprocessed/images/rt_ECoG_cc.png', sep = ''))


nbin = 40
bin_data_cluster = data_frame(trial = c(tapply(avg_cluster$order, cut(avg_cluster$order, nbin), mean), 
                                        tapply(avg_cluster$order, cut(avg_cluster$order, nbin), mean)),
                              mean_rt = c(tapply(filter(avg_cluster, transition == "False")$mean_rt, cut(filter(avg_cluster, transition == "False")$order, nbin), mean), 
                                          tapply(filter(avg_cluster, transition == "True")$mean_rt, cut(filter(avg_cluster, transition == "True")$order, nbin), mean)),
                              transition = c(rep("False", times = length(tapply(avg_cluster$order, cut(avg_cluster$order, nbin), mean))), 
                                                  rep("True", times = length(tapply(avg_cluster$order, cut(avg_cluster$order, nbin), mean)))))


plot = ggplot(data=bin_data_cluster, aes(x=trial, y=mean_rt, color = transition))
plot + geom_line(size=1) + ggtitle('RT over time, by Graph') +
  theme_minimal() + labs(x = 'Trial', y = 'RT (ms)') + scale_color_manual(values = c(rgb(215/255,190/255,123/255), rgb(33/255,67/255,104/255))) +
  ggsave(paste( 'behavior_preprocessed/images/rt_ECoG_bin_cc.pdf', sep = ''))

## Graph
nbin = 40
bin_data_graph= data_frame(trial = c(tapply(avg_data$order, cut(avg_data$order, nbin), mean), 
                                        tapply(avg_data$order, cut(avg_data$order, nbin), mean)),
                              mean_rt = c(tapply(filter(avg_data, graph == "modular")$mean_rt, cut(filter(avg_data, graph == "modular")$order, nbin), mean), 
                                          tapply(filter(avg_data, graph == "lattice")$mean_rt, cut(filter(avg_data, graph == "lattice")$order, nbin), mean)),
                              transition = c(rep("modular", times = length(tapply(avg_data$order, cut(avg_data$order, nbin), mean))), 
                                             rep("lattice", times = length(tapply(avg_data$order, cut(avg_data$order, nbin), mean)))))


plot = ggplot(data=bin_data_graph, aes(x=trial, y=mean_rt, color = transition))
plot + geom_line(size=1) + ggtitle('RT over time, by Graph') +
  theme_minimal() + labs(x = 'Trial', y = 'RT (ms)') + scale_color_manual(values = c(rgb(215/255,190/255,123/255), rgb(33/255,67/255,104/255))) +
  ggsave(paste( 'behavior_preprocessed/images/rt_ECoG_bin_graph.pdf', sep = ''))

