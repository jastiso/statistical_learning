%% Fit maximum entropy model

clear

addpath(genpath('/Users/stiso/Documents/MATLAB/IRASA/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
% define variables
save_dir = '/Users/stiso/Documents/Python/graphLearning/ECoG data/behavior_preprocessed/';
ephys_dir = '/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_raw/';
img_dir = '/Users/stiso/Documents/Python/graphLearning/ECoG data/behavior_preprocessed/images/';


% load stuff
data = readtable([save_dir, 'residuals.csv']);


%% Get model

subjs = [{'2'}, {'4'}, {'6'}];
nSubj = numel(subjs);
beta = zeros(nSubj,1);
r0 = zeros(nSubj,1);
r1 = zeros(nSubj,1);
E = zeros(nSubj,1);
diff = zeros(nSubj,1);

for s = 1:nSubj
    % select subject
    curr = data(strcmpi(data.subj,subjs{s}),:);
    
    rt = data.resid;
    trials = data.order;
    
    load([ephys_dir, subjs{s}, '/task_data.mat'])
    
    [beta(s), r0(s), r1(s), E(s), diff(s)] = learn_linear_real_full(walk, rt, trials);
end
save([save_dir, 'max_ent.mat'], 'beta', 'r0', 'r1', 'E', 'diff')


%% plot
figure(1); clf
histogram(log10(beta), 'facecolor', rgb('steelblue'), 'facealpha', 0.6, 'edgecolor', 'white')
xlabel('log( beta )')
ylabel('Count')
saveas(gca, [img_dir, 'mTurk_betas.png'], 'png')