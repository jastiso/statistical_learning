#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
linear mixed effects model accross subjects

Created on Thu Mar 14 17:12:54 2019

@author: stiso
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import statsmodels.api as sm
import statsmodels.formula.api as smf

# read data
top_dir = '/Users/stiso/Documents/Python/graphLearning/ECoG data/behavioral_data_raw/'
subj = [1,2]

#%% Load and concatenate data

data = pd.DataFrame()
for i in subj:
    curr = pd.read_csv("".join([top_dir, 'subj', str(i), '/subj', str(i), '_clean_data']))
    curr['subj'] = i
    data = data.append(curr)
    
    
#%% Linear Mixed effects model
    
 md = smf.mixedlm("rt_raw ~ transition*order*graph + lag10 + recency + walk + hand + hand_transition", data, groups=data["subj"])

mdf = md.fit()   

mdf.summary()