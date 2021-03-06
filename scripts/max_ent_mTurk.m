%% Fit maximum entropy model

clear

addpath(genpath('/Users/stiso/Documents/MATLAB/IRASA/'))
addpath(genpath('/Users/stiso/Documents/Code/graph_learning/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
% define variables
save_dir = '/Users/stiso/Documents/Python/graphLearning/old_tasks/mTurk-10-node-breaks/data/preprocessed/';
ephys_dir = '/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_raw/';
img_dir = '/Users/stiso/Documents/Python/graphLearning/images/';
r_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_analysis/';

% load stuff
data = readtable([save_dir, 'residuals.csv']);
%data = data(data.is_lattice == 0,:);
walks = readtable('/Users/stiso/Documents/Python/graphLearning/old_tasks/mTurk-10-node-breaks/experiment/stims/nodes_random.csv');

nNode = 10;

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
    trials = curr.cum_trial;
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

%% Get A_hat in blocks

load([save_dir, 'max_ent.mat'], 'beta', 'A_hat')
shift = 100; % how much to slide windows
win = 500; % number of trials per window
nBlock = (1000 - win)/shift + 1;
blocks = zeros(nBlock,2);
st = 1;
for b = 1:nBlock
    blocks(b,:) = [st,st+(win)-1];
    st = st + shift;
end
A_hat_block = zeros(nBlock,nSubj,nNode,nNode);
mat_dist = zeros(nBlock,nSubj);
beta_block = zeros(nBlock,nSubj);
tri_mask = logical(triu(ones(nNode),1));

for s = 1:nSubj
    subj_data = data(strcmp(data.workerid,subjs{s}),:);
    for b = 1:nBlock
        curr = subj_data((subj_data.cum_trial >= blocks(b,1)) & (subj_data.cum_trial <= blocks(b,2)),:);
        rt = curr.resid;
        trials = curr.trial;
        idx = curr.walk_id(1);
        walk = walks(idx,:).Variables;

        [beta_block(b,s), r0(s), r1(s), E(s), diff(s)] = learn_linear_real_full(walk, rt, trials);


    end
end

boxplot(log(mat_dist)')
beta_block(beta_block > 1000) = 1000;
beta_block(beta_block < 0) = 0;
block_idx = [];
beta_list = [];
subj_list = [];
    for m = 1:nBlock
        block_idx = [block_idx; ones(nSubj,1)+(m-1)];
        beta_list = [beta_list; beta];
        subj_list = [subj_list; subjs];
    end
ahat_block_data = table(subj_list, block_idx,...
    reshape(beta_block', [] ,1), beta_list, ...
    'VariableNames', [{'subj'}, {'block'},  {'beta_block'}, {'beta'}]);
writetable(ahat_block_data, [r_dir, 'ahat_block_mturk.csv']);

