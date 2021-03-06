%% Plot mTurk and iEEG betas together

clear
mTurk_dir = '/Users/stiso/Documents/Python/graphLearning/old_tasks/mTurk-10-node-breaks/data/preprocessed/';
ieeg_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/behavior_preprocessed/';

%% Load betas

load([mTurk_dir, 'max_ent.mat'], 'beta')
beta_mturk = beta;
load([ieeg_dir, 'max_ent.mat'], 'beta')
data = readtable('/Users/stiso/Documents/Python/graphLearning/old_tasks/mTurk-10-node-breaks/data/preprocessed/residuals.csv');
subjs = [{'1'}, {'2'}, {'3'}, {'4'}, {'5'}, {'6'}, {'7'}, {'8'}, {'10'}, {'12'},{'18'}];
% get mturk graph idx
mturk_subjs = unique(data.workerid);
mturk_idx = cellfun(@(x) any(data.is_lattice(strcmp(data.workerid,x))), mturk_subjs);

figure(1); clf
[f,xi,bw]=ksdensity(log(beta_mturk(beta_mturk < 1000 & beta_mturk > 0))); hold on
plot(xi,f,'k','linewidth',4)
for b = 1:numel(beta)
    if (beta(b) > 0) && (beta(b) < 1000)
        if mod(str2double(subjs{b}),2) ~= 0
            c = [101/255,111/255,147/255];
        else
            c = [125/255,138/255,95/255];
        end
        plot([log(beta(b)), log(beta(b))], [0,0.1], 'color' ,c, 'linewidth', 2)
    end
end
saveas(gca, [ieeg_dir, 'images/all_beta.pdf'], 'pdf')

fprintf('\nMean for mTurk: %d std for mTurk %d', (mean((beta_mturk(beta_mturk < 1000 & beta_mturk > 0)))), ...
    (std((beta_mturk(beta_mturk < 1000 & beta_mturk > 0)))))
fprintf('\nMean for iEEG: %d std for iEEG %d', (mean((beta(beta < 1000 & beta > 0)))), ...
    (std(log10(beta(beta < 1000 & beta > 0)))))
fprintf('\nNo TD for mTurk: %d ', sum((beta_mturk < 1000 & beta_mturk > 0))/numel(beta_mturk))
fprintf('\nNo TD for iEEG: %d ', sum(~(beta < 1000 & beta > 0))/numel(beta))

[h,p] = ttest2(log10(beta_mturk(beta_mturk < 1000 & beta_mturk > 0)), log10(beta(beta < 1000 & beta > 0)))

%% Split mturk by graph

beta_mod = beta_mturk(~mturk_idx);
beta_lat = beta_mturk(mturk_idx);
figure(1); clf
[f,xi,bw]=ksdensity(log(beta_mod(beta_mod < 1000 & beta_mod > 0))); hold on
plot(xi,f,'color',[125/255,138/255,95/255],'linewidth',4)
[f,xi,bw]=ksdensity(log(beta_lat(beta_lat < 1000 & beta_lat > 0))); 
plot(xi,f,'color',[101/255,111/255,147/255],'linewidth',4)
for b = 1:numel(beta)
    if (beta(b) > 0) && (beta(b) < 1000)
        if mod(str2double(subjs{b}),2) ~= 0
            c = [101/255,111/255,147/255];
        else
            c = [125/255,138/255,95/255];
        end
        plot([log(beta(b)), log(beta(b))], [0,0.1], 'color' ,c, 'linewidth', 2)
    end
end
saveas(gca, [ieeg_dir, 'images/all_beta.pdf'], 'pdf')
[h,p] = ttest2(log10(beta_mod(beta_mod < 1000 & beta_mod > 0)), log10(beta_lat(beta_lat < 1000 & beta_lat > 0)))
