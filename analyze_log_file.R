data = readMat('/Users/stiso/Documents/Python/graphLearning/behavioral data/subj2_beh.mat')
df = data.frame(rt = data$rt[,1], lag10 = data$lag10[,1], recency = data$recency[,1], target = unlist(data$stim.id), order = data$order[,1], edgeType = data$trans.idx[,1])

fit = lm(rt ~ target + edgeType, df)
anova(fit)
summary(fit)
