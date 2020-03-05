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
library(R.matlab)
setwd("/Users/stiso/Documents/Python/graphLearning/old_tasks/mTurk-10-node-breaks")
ext = '1'

df = read.csv('/Users/stiso/Documents/Python/graphLearning/old_tasks/mod_rand.csv')
summary(df)

###########
# Check accuracy across subjects
acc = dplyr::filter(df, nTries == 1, stage == 'walk_one') %>%
  group_by(workerid) %>%
  dplyr::summarise(total_acc = sum(correct == "True")/1500)
p<-ggplot(acc, aes(x=total_acc)) + 
  geom_histogram(fill='lightblue', color='white')
p
summary(acc$total_acc)

subj_rm = df$workerid == acc$workerid[acc$total_acc < .80]
df = df[!subj_rm,]

# remove incorrect trials
df_clean = dplyr::filter(df, stage == 'walk_one')
df_clean = dplyr::filter(df_clean, correct == "True")

####################
# get lag 10
stop_num=10
recency_fun = function(x) {
  if (x > stop_num) {
    y = stop_num+1
  } else {
    y = x
  }
  return(y)
}

f <- function(curr_row, lag10) {
  if (curr_row['trial'] == 0) { # if trial is 0
    # initialize
    lag10 = rep(-1, times = 10)
    # get current node
    curr_node = unlist(curr_row['node']) # node
    
    # lag10
    curr_lag = sum(lag10 == curr_node)
    # update
    lag10 = c(lag10[2:10], curr_node)
    
  } else {
    # get current node
    curr_node = unlist(curr_row['node']) # node
    
    # lag10
    curr_lag = sum(lag10 == curr_node)
    # update
    lag10 = c(lag10[2:10], curr_node)
    
  }
  
  return(list(curr_lag, lag10))

  
}

lag10 = rep(-1, times = 10)
all_lag = rep(0, times = nrow(df_clean))
for(r in 1:nrow(df_clean)) {
  output = f(df_clean[r,], lag10)
  lag10 = unlist(output[[2]])
  all_lag[r] = unlist(output[[1]])
  }
df_clean$lag10 = all_lag

# might need to change this depending on the task
df_clean = subset(df_clean, select = -c(target,phase,event))
df_clean = na.omit(df_clean)
summary(df_clean)

#make factors
cat_vars = c('keyCode', 'workerid', 'correct', 'graph')
for (var in cat_vars){
  df_clean[var] = 
    as.factor(unlist(df_clean[var]))
  
}

# remove RT > 2s
df_clean = filter(df_clean, rt < 2000)
df_clean = filter(df_clean, rt > 50)


finger = character(length = length(df_clean$keyCode))
for (t in 1:length(finger)){
  curr_key = df_clean$keyCode[t]
  if (curr_key == 186 || curr_key == 59){
    finger[t] = 'pinky'
  } else if (curr_key == 76){
    finger[t] = 'ring'
  } else if (curr_key == 75){
    finger[t] = 'middle'
  } else if (curr_key == 74){
    finger[t] = 'index'
  } else if (curr_key == 32){
    finger[t] = 'thumb'
  } else {
    finger[t] = NA
  }
}
df_clean$finger = as.factor(finger)

p<-ggplot(df_clean, aes(x=rt)) + 
  geom_histogram(fill='pink', color='white')
p
# make rts log
df_clean$rt = log10(df_clean$rt)
p<-ggplot(df_clean, aes(x=rt)) + 
  geom_histogram(fill='lightblue', color='white')
p
df_clean = na.omit(df_clean)
df_clean$trial = df_clean$trial + 1
summary(df_clean)


################################################
# effect of lag 10
################################################
nuissance_reg = lmer(data = df_clean, rt~scale(log10(trial))*graph + finger + (1 + scale(log10(trial)) | workerid))
recency_fact = lapply(df_clean$recency, recency_fun)
recency_exp = lapply(df_clean$recency, function(x) {-exp(-x)})
df_clean$recency_fact = unlist(recency_fact)
df_clean$recency_exp = unlist(recency_exp)
recency_data = data.frame( rt = resid(nuissance_reg), rt_raw = df_clean$rt, lag10 = as.factor(df_clean$lag10), recency = df_clean$recency, graph = as.factor(df_clean$graph), 
                           recency_exp = df_clean$recency_exp, subj = df_clean$workerid, recency_fact = unlist(recency_fact))
p = ggplot(data=recency_data, aes(y=rt, x=lag10, color = graph)) + geom_jitter(alpha=0.3) + theme_minimal()
p
ggsave(paste( 'experiment/data/preprocessed/images/mod_rand.png', sep = ''))


recency_avg <- recency_data %>%
  group_by(subj, recency_exp,graph) %>%
  dplyr::summarise(mean_rt = mean(rt)) 
p = ggplot(data=recency_avg, aes(y=mean_rt, x=recency_exp, color = graph)) + geom_jitter(width = 0.1, alpha=0.5, size=3)  + theme_minimal()
p + geom_smooth(method = 'lm') #+ geom_line(aes(group=subj))#+ theme(legend.position = "none")
ggsave(paste( 'experiment/data/preprocessed/images/mod_rand_rec_exp.png', sep = ''))



# fit
fit = lm(data=recency_data, rt~log(recency))
summary(fit)

p = ggplot(data=recency_avg, aes(y=mean_rt, x=(recency_fact), color = graph)) + geom_jitter(alpha=0.5, size=3)  + theme_minimal() + geom_smooth(method = 'loess')
p




# Compare fits
recency_fit = lmer(data = df_clean, rt~scale(log10(trial))*graph + scale(log(recency_fact)) + finger + (1 + scale(log10(trial)) | workerid))
summary(recency_fit)

full_recency_fit = lmer(data = df_clean, rt~scale(log10(trial))*graph + scale(log(recency)) + finger + (1 + scale(log10(trial)) | workerid))
summary(full_recency_fit)

recency_exp_fit = lmer(data = df_clean, rt~scale(log10(trial))*graph + scale(recency_exp) + finger + (1 + scale(log10(trial)) | workerid))
summary(recency_exp_fit)

orig_fit = lmer(data = df_clean, rt~scale(log10(trial))*graph + scale(recency) + scale(lag10) + finger + (1 + scale(log10(trial)) | workerid))
summary(orig_fit)

anova(orig_fit, recency_fit, test="Chisq")
anova(full_recency_fit, recency_fit, test="Chisq")
anova(recency_exp_fit, recency_fit, test="Chisq")


