%% Plot mTurk and iEEG betas together

clear
mTurk_dir = '/Users/stiso/Documents/Python/graphLearning/old_tasks/mTurk-10-node-breaks/data/preprocessed/';
ieeg_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/behavior_preprocessed/';

%% Load betas

load([mTurk_dir, 'max_ent.mat'], 'beta')
beta_mturk = beta;
load([ieeg_dir, 'max_ent.mat'], 'beta')

figure(1); clf
ksdensity(beta_mturk(beta_mturk < 1000 & beta_mturk > 0)); hold on
for b = 1:numel(beta)
  if (beta(b) > 0) && (beta(b) < 1000) 
      plot([beta(b), beta(b)], [0,0.5], 'r', 'linewidth', 2)
  end
end
saveas(gca, [ieeg_dir, 'images/all_beta.png'], 'png')