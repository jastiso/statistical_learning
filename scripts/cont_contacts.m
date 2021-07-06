%% Count removed elecs

addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current/'))
addpath(genpath('/Users/stiso/Documents/Code/graph_learning/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')

% define variables
subj = '18';
sess = '_sess1';
save_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_raw/';

load([save_dir, subj, '/header', sess,'.mat'], 'elec_labels', 'srate', 'HUP_ID', 'subj')
origc = elec_labels;
load([save_dir, subj, '/header_clean.mat'], 'elec_labels', 'srate', 'HUP_ID', 'subj', 'regions', 'sessions')

fprintf('Total contacs: %d, percent removed %d', numel(elec_labels), round((numel(origc) - numel(elec_labels))/numel(elec_labels)*100,2))