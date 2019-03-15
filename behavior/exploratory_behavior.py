#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 11:53:01 2018

@author: stiso
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# read data
top_dir = '/Users/stiso/Documents/Python/graphLearning/'
subj = '999'
data_raw = pd.read_csv("".join((top_dir, 'behavioral data/subj', subj, '_log_motor_run1_10node.csv')))

# cut down to only relevant fields
# get only correct responses
data = data_raw.loc[data_raw.correct_raw == 1]
data = data.dropna(how='any')
# convert rt from string
#data.rt_raw = data.rt_raw.apply(lambda x: x.strip('\''))
data.rt_raw = data.rt_raw.apply(float)
# remove rt's greater than 2 seconds
data = data.loc[data.rt_raw < 2]

# plot
num_bins = 50
fig = plt.figure()
# the histogram of the data
plot_data = data.rt_raw.tolist()
n, bins, patches = plt.hist(plot_data, num_bins, facecolor='blue', alpha=0.5)
plt.xlabel('RT')
plt.ylabel('Frequency')
fig.savefig("".join([top_dir, 'behavioral data/subj', subj, '_rt_hist.png']))

fig = plt.figure()
plt.plot(plot_data)
plt.xlabel('Trial')
plt.ylabel('RT')
fig.savefig("".join([top_dir, 'behavioral data/subj', subj, '_rt.png']))

#%% Get transition data

#transitions = [0,4,5,9,10,14]
transitions = [0,4,5,9]
trans_idx = []
for i in data.walk.tolist():
    trans_idx.append(i in transitions)
    
transition_data = data.loc[trans_idx]

# append transition to data
data['transition'] = trans_idx
data['transition'] = data['transition'].astype('category')

# plot
num_bins = 50
fig = plt.figure()
# the histogram of the data
plot_data = transition_data.rt_raw.tolist()
n, bins, patches = plt.hist(plot_data, num_bins, facecolor='orange', alpha=0.5)
plt.xlabel('RT')
plt.ylabel('Frequency')
fig.savefig("".join([top_dir, 'behavioral data/subj', subj, '_rt_hist_trans.png']))

fig = plt.figure()
plt.plot(plot_data)
plt.xlabel('Trial')
plt.ylabel('RT')
fig.savefig("".join([top_dir, 'behavioral data/subj', subj, '_rt_trans.png']))


#%% Not transition data

trans_idx = []
for i in data.walk.tolist():
    trans_idx.append(i not in transitions)
    
within_data = data.loc[trans_idx]

# plot
num_bins = 50
fig = plt.figure()
# the histogram of the data
plot_data = within_data.rt_raw.tolist()
n, bins, patches = plt.hist(plot_data, num_bins, facecolor='red', alpha=0.5)
plt.xlabel('RT')
plt.ylabel('Frequency')
fig.savefig("".join([top_dir, 'behavioral data/subj', subj, '_rt_hist_within.png']))

fig = plt.figure()
plt.plot(plot_data)
plt.xlabel('Trial')
plt.ylabel('RT')
fig.savefig("".join([top_dir, 'behavioral data/subj', subj, '_rt_within.png']))


#%% Stats

print(transition_data.rt_raw.mean())
print(within_data.rt_raw.mean())
print(data.rt_raw.mean())

# make predictors for recency and priming, according to karuza paper
# lag10
lag = 10
walk = data.walk.astype(float)
walk = walk.tolist()
# append zeros to the beginning
for i in range(lag):
   walk.insert(0,np.nan)
lag10 = [];
for i,n in enumerate(walk):
    if i > 9:
        curr = walk[(i-10):(i-1)]
        lag10.append(curr.count(n))
    
# recency
walk = data.walk.astype(float)
walk = walk.tolist()
recency = [np.nan]
for i,n in enumerate(walk):
    if i > 0:
        curr = walk[0:(i-1)]
        if n in curr:
            idx = i - max(loc for loc, val in enumerate(curr) if val == n)
        else:
            idx = np.nan
        recency.append(idx)

import statsmodels.api as sm

data['lag10'] = lag10
data['recency'] = recency
data['resp_raw'] = data['resp_raw'].astype('category')
model_data = data.dropna(how='any')
Y = model_data["rt_raw"]
X = model_data[["transition", "order", "walk"]]
X = sm.add_constant(X)
#X["resp_raw"].cat.categories = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
#X["resp_raw"].cat.categories = [1,2,3,4,5,6,7,8,9,10]

# only first 500 trials
#X = X[:500:]
#Y = Y[:500]
# Note the difference in argument order
model = sm.OLS(Y, X.astype(float)).fit()
predictions = model.predict(X) # make the predictions by the model

# Print out the statistics
model.summary()


#%% Check for colinearity

#from libraries.settings import *
from scipy.stats.stats import pearsonr
import itertools

correlations = {}
columns = X.columns.tolist()

for col_a, col_b in itertools.combinations(columns, 2):
    correlations[col_a + '__' + col_b] = pearsonr(X.loc[:, col_a], X.loc[:, col_b])

result = pd.DataFrame.from_dict(correlations, orient='index')
result.columns = ['PCC', 'p-value'] 

print(result)
