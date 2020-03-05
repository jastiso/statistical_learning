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
keys = ['q', 'w', 'e', 'r', 'v', 'b', 'u', 'i', 'o', 'p']

# create stimuli png list
stim = []
stim_prac = []
resp = [0,1,2,3,4,5,6,7,8,9]
for i in range(nNode):
    stim.append("".join(('img/target_', str(resp[i]+1), '.png')))
    stim_prac.append("".join(('img/target_', str(resp[i]+1), '_prac.png')))
    

#%% Loop through subjects
    
for i in range(nSubj-1):
    # save dir: consistent with Steve's naming scheme
    save_dir = "".join(['/Users/stiso/Documents/Python/graphLearning/subj', str(i+1), '/'])
    if not os.path.exists(save_dir):
            os.makedirs(save_dir)
    else:
        print("You have already made these files! delete them first and rerun")
        break
            
    # add first row of all CSVs
    csvfile = "".join([save_dir, 'exposure_walk.csv'])
    with open(csvfile, "a") as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(('pID', 'trialNum', 'graph', 'walk', 'resp','path'))
    csvfile = "".join([save_dir, 'image_mapping.csv'])
    with open(csvfile, "a") as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(('sID','imageNums', 'imageFilenames'))
    csvfile = "".join([save_dir, 'prac_walk.csv'])
    with open(csvfile, "a") as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(('pID', 'trialNum', 'walk', 'resp','path'))
    
    # create random mapping: randomize list of 15 stimuli
    all_stim = list(zip(stim,keys,stim_prac))
    random.shuffle(all_stim)
    stim,keys,stim_prac = zip(*all_stim)
    
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
    all_prac = [sID, trialNum, walk, keys, stim_prac]
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
        full_resp = []
        for j in walks:
            full_stim.append(stim[int(j)-1])
            full_resp.append(keys[int(j)-1])
        
        # write to csv
        graph = ['lattice']*len(full_stim)
        sID = [i+1]*len(full_stim)
        all_exp = [sID, range(1,len(full_stim)+1), graph, walks, full_resp, full_stim]
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
        full_resp = []
        for j in walks:
            full_stim.append(stim[int(j)-1])
            full_resp.append(keys[int(j)-1])
        
        # write to csv
        graph = ['modular']*len(full_stim)
        sID = [i+1]*len(full_stim)
        all_exp = [sID, range(1,len(full_stim)+1), graph, walks, full_resp, full_stim]
        l = zip(*all_exp)
        csvfile = "".join([save_dir, 'exposure_walk.csv'])
        with open(csvfile, "a") as fp:
            wr = csv.writer(fp, dialect='excel')
            for entry in l:
                wr.writerow(entry)