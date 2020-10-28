%% Fit maximum entropy model

clear

addpath(genpath('/Users/stiso/Documents/Code/graph_learning/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
% define variables
save_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/behavior_preprocessed/';
ephys_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_raw/';
img_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/behavior_preprocessed/images/';
r_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_analysis/';

% load stuff
data = readtable([save_dir, 'residuals.csv']);


%% Get model

subjs = [{'1'}, {'2'}, {'3'}, {'4'}, {'6'}, {'8'}, {'10'}, {'12'},{'18'}];
subjs = [{'18'}];
nSubj = numel(subjs);
nNode = 10;
beta = zeros(nSubj,1);
r0 = zeros(nSubj,1);
r1 = zeros(nSubj,1);
E = zeros(nSubj,1);
diff = zeros(nSubj,1);
M = [0 1 1 1 0 0 0 0 0 1;
    1 0 1 1 1 0 0 0 0 0;
    1 1 0 1 1 0 0 0 0 0;
    1 1 1 0 1 0 0 0 0 0;
    0 1 1 1 0 1 0 0 0 0;
    0 0 0 0 1 0 1 1 1 0;
    0 0 0 0 0 1 0 1 1 1;
    0 0 0 0 0 1 1 0 1 1;
    0 0 0 0 0 1 1 1 0 1
    1 0 0 0 0 0 1 1 1 0].*0.25;
L = [0 1 1 0 0 0 0 0 1 1;
    1 0 1 1 0 0 0 0 0 1;
    1 1 0 1 1 0 0 0 0 0;
    0 1 1 0 1 1 0 0 0 0;
    0 0 1 1 0 1 1 0 0 0;
    0 0 0 1 1 0 1 1 0 0;
    0 0 0 0 1 1 0 1 1 0;
    0 0 0 0 0 1 1 0 1 1;
    1 0 0 0 0 0 1 1 0 1
    1 1 0 0 0 0 0 1 1 0].*0.25;
A_hat = zeros(nSubj, nNode, nNode);

for s = 1:nSubj
    if ~exist([r_dir, 'subj', subjs{s}],'dir')
       mkdir([r_dir, 'subj', subjs{s}]); 
    end
    
    % select subject
    fprintf('\n************** Subj %s ************\n', subjs{s});
    curr = data(data.subj==str2double(subjs{s}),:);
    
    if mod(str2double(subjs{s}),2) == 0
        A = M;
    else
        A = L;
    end
    
    rt = curr.resid;
    trials = curr.order;
    
    load([ephys_dir, subjs{s}, '/task_data.mat'])
    
    [beta(s), r0(s), r1(s), E(s), diff(s)] = learn_linear_real_full(walk, rt, trials);
    
    % get A_hat
    A_hat(s, :, :) = (1 - exp(-beta(s)))*A*(eye(nNode) - exp(-beta(s))*A)^(-1);
    
    figure(1); clf
    imagesc(squeeze(A_hat(s,:,:))); colorbar
    saveas(gca, [img_dir, 'A_hat_', num2str(s), '.png'], 'png')

end
save([save_dir, 'max_ent.mat'], 'beta', 'r0', 'r1', 'E', 'diff', 'A_hat')


%% plot

figure(1); clf
histogram((beta(beta < 900)), 30, 'facecolor', rgb('steelblue'), 'facealpha', 0.6, 'edgecolor', 'white')
xlabel('log( beta )')
ylabel('Count')
saveas(gca, [img_dir, 'ecog_betas.png'], 'png')


%% Save A_hat contrast

for s = 1:nSubj
    % select subject
    fprintf('\n************** Subj %s ************\n', subjs{s});
    curr = squeeze(A_hat(s,:,:));
    
    load([ephys_dir, subjs{s}, '/task_data.mat'])
    load([ephys_dir, subjs{s}, '/good_events.mat'])
    
    max_ent_cont = zeros(size(walk));
    
    transitions = [walk; circshift(walk,1)] + 1; % offset 0 index
    for i = 2:numel(walk)
        max_ent_cont(i) = curr(transitions(1,i),transitions(2,i));
    end
    max_ent_cont = max_ent_cont(good_trials);
    save([r_dir, 'subj', subjs{s}, '/max_ent_contrast.mat'], 'max_ent_cont')
end

