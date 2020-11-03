%% Plot mTurk and iEEG betas together

clear
mTurk_dir = '/Users/stiso/Documents/Python/graphLearning/old_tasks/mTurk-10-node-breaks/data/preprocessed/';
ieeg_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/behavior_preprocessed/';

%% Load betas

load([mTurk_dir, 'max_ent.mat'], 'beta')
beta_mturk = beta;
load([ieeg_dir, 'max_ent.mat'], 'beta')

figure(1); clf
[f,xi,bw]=ksdensity(log(beta_mturk(beta_mturk < 1000 & beta_mturk > 0))); hold on
plot(xi,f,'linewidth',4)
for b = 1:numel(beta)
  if (beta(b) > 0) && (beta(b) < 1000) 
      plot([log(beta(b)), log(beta(b))], [0,0.1], 'r', 'linewidth', 2)
  end
end
saveas(gca, [ieeg_dir, 'images/all_beta.png'], 'png')

fprintf('\nMean for mTurk: %d std for mTurk %d', (mean((beta_mturk(beta_mturk < 1000 & beta_mturk > 0)))), ...
    (std((beta_mturk(beta_mturk < 1000 & beta_mturk > 0)))))
fprintf('\nMean for iEEG: %d std for iEEG %d', (mean((beta(beta < 1000 & beta > 0)))), ...
    (std(log10(beta(beta < 1000 & beta > 0)))))
fprintf('\nNo TD for mTurk: %d ', sum(~(beta_mturk < 1000 & beta_mturk > 0))/numel(beta_mturk))
fprintf('\nNo TD for iEEG: %d ', sum(~(beta < 1000 & beta > 0))/numel(beta))

permtest(beta_mturk(beta_mturk < 1000 & beta_mturk > 0), beta(beta > 1000 & beta < 0))
