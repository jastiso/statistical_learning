#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 14:13:15 2018

For every subject, make csvs with all the info relevant for generating walks for their experiment

This includes: 
    image_mapping.csv - mapping from nodes to pngs
    prac_walk - practice order with images
    exposure_walk - full task with images

@author: stiso
"""
import random
import csv
import os
import math

nNode = 10
nSubj = 49

# Key presses
keysR = ['space', 'j', 'k', 'l', 'semicolon']
keysL = ['a', 's', 'd', 'f', 'space']

# create stimuli png list
stim = []
if nNode == 15:
    resp = [[0],[1],[2],[3],[4],[0,1],[0,2],[0,3],[0,4],[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]
    resp_idx = range(nNode)
else:
    resp = [[0,1],[0,2],[0,3],[0,4],[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]
    resp_idx = [5,6,7,8,9,10,11,12,13,14]
respR = []
respL = []
for i in range(nNode):
    stim.append("".join(('img/target_', str(resp_idx[i]+1), '.png')))
    tmpR = []
    tmpL = []
    for j in resp[i]:
        tmpR.append(keysR[j])
        tmpL.append(keysL[j])
    respR.append(tmpR)
    respL.append(tmpL)

#%% Loop through subjects
    
for i in range(nSubj-1):
    # save dir: consistent with Steve's naming scheme
    save_dir = "".join(['/Users/stiso/Documents/Python/graphLearning/subj', str(i+1), '/'])
    if not os.path.exists(save_dir):
            os.makedirs(save_dir)
            
    # add first row of all CSVs
    csvfile = "".join([save_dir, 'exposure_walk.csv'])
    with open(csvfile, "a") as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(('pID', 'trialNum', 'graph', 'walk', 'respR','respL','path'))
    csvfile = "".join([save_dir, 'image_mapping.csv'])
    with open(csvfile, "a") as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(('sID','imageNums', 'imageFilenames'))
    csvfile = "".join([save_dir, 'prac_walk.csv'])
    with open(csvfile, "a") as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(('pID', 'trialNum', 'walk', 'respR','respL','path'))
    
    # create random mapping: randomize list of 15 stimuli
    all_stim = list(zip(stim,respR,respL))
    random.shuffle(all_stim)
    stim,respR,respL = zip(*all_stim)
    
    # write to csv
    csvfile = "".join([save_dir, 'image_mapping.csv'])
    sID = [i+1]*nNode
    all_mapping = [sID, range(1, nNode+1), stim]
    l = zip(*all_mapping)
    with open(csvfile, "a") as fp:
        wr = csv.writer(fp, dialect='excel')
        for entry in l:
            print(entry)
            wr.writerow(entry)
        
    # practice trials just go 1:15 for now
    trialNum = range(1, nNode+1)
    walk = range(1, nNode+1)
    
    # save to csv
    all_prac = [sID, trialNum, walk, respR, respL, stim]
    l = zip(*all_prac)
    csvfile = "".join([save_dir, 'prac_walk.csv'])
    with open(csvfile, "a") as fp:
        wr = csv.writer(fp, dialect='excel')
        for entry in l:
            wr.writerow(entry)
        
    if i%2 == 0:
        
        # if even, use ring lattice
        if nNode == 10:
            filename = "/Users/stiso/Documents/Python/graphLearning/ring_lattice10.csv"
        else:
            filename = "/Users/stiso/Documents/Python/graphLearning/ring_lattice.csv"
        with open(filename) as f:
            reader = csv.reader(f)
            cnt = 0
            for row in reader:
                if cnt == math.floor(i/2):
                    walks = row
                    break
                cnt += 1
        # if walks is a list of 1, break up at semi colon
        if len(walks) == 1:
            walks = walks[0].split(';')
        
        # get mapping to image files and responses
        full_stim = [];
        full_respR = []
        full_respL = []
        for j in walks:
            full_stim.append(stim[int(j)-1])
            full_respR.append(respR[int(j)-1])
            full_respL.append(respL[int(j)-1])
        
        # write to csv
        graph = ['lattice']*len(full_stim)
        sID = [i+1]*len(full_stim)
        all_exp = [sID, range(1,len(full_stim)+1), graph, walks, full_respR, full_respL, full_stim]
        l = zip(*all_exp)
        csvfile = "".join([save_dir, 'exposure_walk.csv'])
        with open(csvfile, "a") as fp:
            wr = csv.writer(fp, dialect='excel')
            for entry in l:
                wr.writerow(entry)
    else:
       # if odd, use modular
        if nNode == 10:
            filename = "/Users/stiso/Documents/Python/graphLearning/modular10.csv"
        else:
            filename = "/Users/stiso/Documents/Python/graphLearning/modular.csv"
        with open(filename) as f:
            reader = csv.reader(f)
            cnt = 0
            for row in reader:
                if cnt == math.floor(i/2):
                    walks = row
                    break
                cnt += 1
        # if walks is a list of 1, break up at semi colon
        if len(walks) == 1:
            walks = walks[0].split(';')
        
        # get mapping to image files and responses
        full_stim = [];
        full_respR = []
        full_respL = []
        for j in walks:
            full_stim.append(stim[int(j)-1])
            full_respR.append(respR[int(j)-1])
            full_respL.append(respL[int(j)-1])
        
        # write to csv
        graph = ['modular']*len(full_stim)
        sID = [i+1]*len(full_stim)
        all_exp = [sID, range(1,len(full_stim)+1), graph, walks, full_respR, full_respL, full_stim]
        l = zip(*all_exp)
        csvfile = "".join([save_dir, 'exposure_walk.csv'])
        with open(csvfile, "a") as fp:
            wr = csv.writer(fp, dialect='excel')
            for entry in l:
                wr.writerow(entry)