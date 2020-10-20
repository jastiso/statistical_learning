%% Plot mTurk and iEEG betas together

mTurk_dir = '/Users/stiso/Documents/Python/graphLearning/old_tasks/mTurk-10-node-breaks/data/preprocessed/';
ieeg_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/behavior_preprocessed/';

%% Load betas

load([mTurk_dir, 'max_ent.mat'], 'beta')
beta_mturk = beta;
load([ieeg_dir, 'max_ent.mat'], 'beta')

figure(1); clf
histogram(beta_mturk(beta_mturk < 1000),10)
