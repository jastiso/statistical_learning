library(ggplot2)
library(R.matlab)
library(dplyr)
library(coin)
library(lmPerm)
library(car)
library(aplpack)
library(lmerTest)
library(RColorBrewer)
library(wesanderson)
library(ez)
setwd("/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_analysis/")

subj = '2'
df = readMat(paste('subj', subj, '/psds.mat', sep = ''))
psd = data.frame(power = unlist(df$psds.vect), freq = df$freq.id, elec = df$elec.id, trans = df$trans.exp, alt = df$alt.exp, trial = df$trial)

theta = c(4,10)
beta = c(24,30)

psd = mutate(psd, type = paste(alt,trans))
psd = dplyr::filter(psd, type == "0 1" | type == "1 0")
# get only hipp
psd = dplyr::filter(psd, elec >= 24 & elec <= 35)

psd_theta = dplyr::filter(psd, freq > theta[1] & freq <= theta[2])
psd_theta = psd_theta %>%
  group_by(elec, trial, type) %>%
  dplyr::summarise(mean_power = mean(power))
psd_beta = filter(psd, freq > beta[1] & freq <= beta[2])
psd_beta = psd_beta %>%
  group_by(elec, trial, type) %>%
  dplyr::summarise(mean_power = mean(power))



plot = ggplot(psd_theta, aes(x = as.factor(elec), y = mean_power, fill = type) )
plot + geom_violin(trim = FALSE, position = position_dodge(0.75)) + 
  geom_boxplot(width = 0.15,
               position = position_dodge(0.75)
  ) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5) + 
  labs(x = 'Type', y = 'log theta power')  + theme_minimal() #+ scale_fill_manual(values = wes_palette("Royal1"))
ggsave(paste('theta_type', '.png', sep = ''))

plot = ggplot(psd_beta, aes(x = as.factor(elec), y = mean_power, fill = as.factor(type)) )
plot + geom_violin(trim = FALSE, position = position_dodge(0.75)) + 
  geom_boxplot(width = 0.15,
               position = position_dodge(0.75)
  ) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.05, position = position_dodge(0.75), alpha = 0.5) + 
  labs(x = 'elec', y = 'log beta power')  + theme_minimal() #+ scale_fill_manual(values = wes_palette("Royal1"))
ggsave(paste('beta_type', '.png', sep = ''))


## single elec
plot = ggplot(dplyr::filter(psd_theta, elec == 33), aes(x = type, y = mean_power, fill = type) )
plot + geom_violin(trim = FALSE, position = position_dodge(0.75)) + 
  geom_boxplot(width = 0.15,
               position = position_dodge(0.75)
  ) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.5) + 
  labs(x = 'Type', y = 'log theta power')  + theme_minimal() + scale_fill_manual(values = c(rgb(33/255,67/255,104/255), rgb(215/255,190/255,123/255)))
ggsave(paste('theta_single', '.pdf', sep = ''))



########################################
## Stats from each elec

df = readMat(paste('subj', subj, '/stats.mat', sep = ''))
stat = data.frame(p = c(df$adj.p, df$adj.pb), freq = c(rep('theta', times = length(df$adj.p)), rep('beta', times = length(df$adj.pb))), elec = c(df$e, df$e))
#stat = dplyr::filter(stat, elec >=24 & elec <= 35)

plot = ggplot(stat, aes(x = elec, y = (p), fill = freq, group = freq))
plot + geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c(rgb(53/255,57/255,61/255), rgb(149/255,192/255,76/255))) +
  geom_hline(yintercept = (0.05), linetype = "dashed", color = "red") +
  #geom_errorbar(aes(ymin = mean_p - sd_p, ymax = (mean_p) + (sd_p)), width=0.2, position=position_dodge(.9)) +
  labs(x = 'Elec', y = 'FDR adjusted P-value')  + theme_minimal() 
ggsave(paste('theta_beta_p.pdf', sep = ''))


