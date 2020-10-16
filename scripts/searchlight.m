%% Compare theta power at within vs between cluster trantisions

clear
clc

addpath(genpath('/Users/stiso/Documents/MATLAB/IRASA/'))
addpath(genpath('/Users/stiso/Documents/MATLAB/Colormaps/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
addpath(genpath('/Users/stiso/Documents/MATLAB/BrainNetViewer_20171031/'))
addpath(genpath('/Users/stiso/Documents/Code/graph_learning/functions/'))

% define variables
subjs = [{'2'},{'1'}, {'4'}, {'6'}, {'8'}, {'10'},{'12'}, {'3'}];
feat_type = 'lfp'; % pow or lfp
freqs = logspace(log10(3), log10(150), 50);

trans1 = [253,224,239]./255;
trans2 = [230,245,208]./255;
within1 = [233,163,201]./255;
within2 = [161,215,106]./255;
center1 = [197,27,125]./255;
center2 = [77,146,33]./255;
cluster_colors = [trans1; within1; center1; within1; trans1; trans2; within2; center2; within2; trans2];

parfor subj_idx = 1:numel(subjs)
    subj = subjs{subj_idx};
    
    searchlight_rsa(subj, subj_idx, cluster_colors, 100, feat_type)
end

