%% Compare theta power at within vs between cluster trantisions

clear
clc

addpath(genpath('/Users/stiso/Documents/MATLAB/IRASA/'))
addpath(genpath('/Users/stiso/Documents/MATLAB/Colormaps/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
addpath(genpath('/Users/stiso/Documents/MATLAB/BrainNetViewer_20171031/'))
addpath(genpath('/Users/stiso/Documents/Code/graph_learning/functions/'))
save_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_raw/';
r_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_analysis/';

% define variables
subjs = [{'1'}, {'2'}, {'3'}, {'4'}, {'5'}, {'6'}, {'7'}, {'8'},{'10'}, {'12'}];
A_hat_order = load([r_dir, 'ahat_order.mat']);
A_hat_order = A_hat_order.subjs;
feat_type = 'lfp_end'; % pow or lfp (_end. _mid)
freqs = logspace(log10(3), log10(150), 50);

trans1 = [253,224,239]./255;
trans2 = [230,245,208]./255;
within1 = [233,163,201]./255;
within2 = [161,215,106]./255;
center1 = [197,27,125]./255;
center2 = [77,146,33]./255;
cluster_colors = [trans1; within1; center1; within1; trans1; trans2; within2; center2; within2; trans2];

get_stim_mapping(subjs)

parfor s = 1:numel(subjs)

    subj = subjs{s};
    subj_idx = find(cellfun(@(x) strcmp(subj,x),A_hat_order));
    
    searchlight_rsa(subj, subj_idx, cluster_colors, 100, feat_type)
end

