# Simulate data 
# lmer(data=df_modular, rt~scale(log10(order))*transition + finger + points + typing_raw + hand + hand_transition +  block + scale(log(recency_fact)) + sess + 
# (1 + scale(log10(order))*transition + scale(log(recency_fact)) |subj))
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

# add all betas for fixed effects
order_beta = -20
transition_beta = 11
points_beta = -10
finger_ind = -50
finger_mid = -40
finger_rin = -30
finger_pin = -20
typing_beta = 10
hand_r_beta = 10
hand_tran_beta = 50
block_beta = -10
recency_beta = 40
sess_beta = -30
trans_order_beta = 1

# for random effects
subj_int_mu = 0
sub_int_sigma = 100
subj_trans_mu = 0
subj_trans_sigma = 5
subj_rec_mu = 0
subj_rec_sigma = 7

# load data
setwd("/Users/stiso/Documents/Code/graph_learning/ECoG_data/")
load('behavior_preprocessed/clean.RData')
df_correct = dplyr::filter(df_correct, graph == "modular")
  
get_rt <- function(df,mean_rt, subj_int, subj_trans, subj_rec){
  if (df$finger == 'index') {
    curr_finger = finger_ind
  } else if (df$finger == 'middle') {
    curr_finger = finger_mid
  } else if (df$finger == 'ring') {
    curr_finger = finger_rin
  } else if (df$finger == 'pinky') {
    curr_finger = finger_pin
  } else {
    curr_finger = 0
  }
  
  if (df$transition == 'True'){
    curr_order_beta = order_beta+trans_order_beta
    curr_transition = transition_beta
  } else {
    curr_order_beta = order_beta
    curr_transition = 0
  }
  
  if (df$points == TRUE) {
    curr_points = points_beta
  } else {
    curr_points = 0
  }
  
  if (df$hand == 'right'){
    curr_hand = hand_r_beta
  } else {
    curr_hand = 0
  }
  
  if (df$hand_transition == TRUE){
    curr_hand_tran = hand_tran_beta
  } else {
    curr_hand_tran = 0
  }
  s = unique(df_correct$subj)
  curr_offset = subj_int[which(df$subj == s)]
  curr_trans_offset = subj_trans[which(df$subj == s)]
  curr_rec_offset = subj_rec[which(df$subj == s)]
  rt = mean_rt + log(df$order)*curr_order_beta + (curr_transition+curr_trans_offset) + 
    curr_finger + curr_points + df$typing_raw*typing_beta + curr_hand + curr_hand_tran + df$block*block_beta + 
    log10(df$recency_fact)*(recency_beta+curr_rec_offset) + curr_offset
  return(rt)
}

subj_int = rnorm(length(unique(df_correct$subj)),subj_int_mu, sub_int_sigma)
subj_trans = rnorm(length(unique(df_correct$subj)),subj_trans_mu, subj_trans_sigma)
subj_rec = rnorm(length(unique(df_correct$subj)),subj_rec_mu, subj_rec_sigma)

rt = list()
for (row in 1:nrow(df_correct)){
  curr = get_rt(df_correct[row,],800, subj_int, subj_trans, subj_rec)
  rt = c(rt,curr)
}
df_correct$sim_rt = unlist(rt)

avg_cluster = filter(df_correct, sess == '1') %>%
  group_by(order, transition) %>%
  dplyr::summarise(mean_rt = mean(sim_rt), sd_rt = sd(sim_rt))

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
  ggsave(paste( 'behavior_preprocessed/images/rt_ECoG_bin_cc_sim.png', sep = ''))




stat_surprisal = lmer(data=df_correct, sim_rt~scale(log10(order))*transition + finger + points + typing_raw + hand + hand_transition + block + scale(log(recency_fact)) + sess + 
                         (1 + scale(log10(order))*transition + scale(log(recency_fact)) |subj))
anova(stat_surprisal)



