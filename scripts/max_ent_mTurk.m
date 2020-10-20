%% Fit maximum entropy model

clear

addpath(genpath('/Users/stiso/Documents/MATLAB/IRASA/'))
addpath(genpath('/Users/stiso/Documents/Code/graph_learning/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
% define variables
save_dir = '/Users/stiso/Documents/Python/graphLearning/old_tasks/mTurk-10-node-breaks/data/preprocessed/';
ephys_dir = '/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_raw/';
img_dir = '/Users/stiso/Documents/Python/graphLearning/images/';

% load stuff
data = readtable([save_dir, 'residuals.csv']);
%data = data(data.is_lattice == 0,:);
walks = readtable('/Users/stiso/Documents/Python/graphLearning/old_tasks/mTurk-10-node-breaks/experiment/stims/nodes_random.csv');

%% Get model

subjs = unique(data.workerid);
nSubj = numel(subjs);
beta = zeros(nSubj,1);
r0 = zeros(nSubj,1);
r1 = zeros(nSubj,1);
E = zeros(nSubj,1);
diff = zeros(nSubj,1);

for s = 1:nSubj
    % select subject
    curr = data(strcmpi(data.workerid,subjs{s}),:);
    % get rt (residuals) and trial
    rt = curr.resid;
    trials = curr.trial;
    % get walk from csv
    idx = curr.walk_id(1);
    walk = walks(idx,:).Variables;
    
    [beta(s), r0(s), r1(s), E(s), diff(s)] = learn_linear_real_full(walk, rt, trials);
end
save([save_dir, 'max_ent.mat'], 'beta', 'r0', 'r1', 'E', 'diff')

%% Plot

plot_beta = beta;
plot_beta(beta == 1000) = 4;
plot_beta(beta == 0) = -1;

figure(1); clf
histogram((plot_beta), 20, 'facecolor', rgb('steelblue'), 'facealpha', 0.6, 'edgecolor', 'white')
xlabel('log( beta )')
ylabel('Count')
saveas(gca, [img_dir, 'mTurk_betas.png'], 'png')
