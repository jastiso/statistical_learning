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

subjs = [{'1'}, {'2'}, {'3'}, {'4'}, {'5'}, {'6'}, {'8'}, {'10'}, {'12'},{'18'}];
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

%% get model
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

%% Get A_hat in blocks

load([save_dir, 'max_ent.mat'], 'beta', 'A_hat')
nBlock = 2;
blocks = [1, 2; 3,4];
A_hat_block = zeros(nBlock,nSubj,nNode,nNode);
mat_dist = zeros(nBlock,nSubj);
beta_block = zeros(nBlock,nSubj);
tri_mask = logical(triu(ones(nNode),1));

for s = 1:nSubj
    subj_data = data(data.subj==str2double(subjs{s}),:);
    for b = 1:nBlock
        curr = subj_data((subj_data.block == blocks(b,1)) | (subj_data.block == blocks(b,2)),:);
        if mod(str2double(subjs{s}),2) == 0
            A = M;
        else
            A = L;
        end

        rt = curr.resid;
        trials = curr.order;

        load([ephys_dir, subjs{s}, '/task_data.mat'])

        [beta_block(b,s), r0(s), r1(s), E(s), diff(s)] = learn_linear_real_full(walk, rt, trials);

        % get A_hat
        A_hat_block(b, s, :, :) = (1 - exp(-beta_block(b,s)))*A*(eye(nNode) - exp(-beta_block(s))*A)^(-1);

        figure(1); clf
        imagesc(squeeze(A_hat_block(b,s,:,:))); colorbar
        saveas(gca, [img_dir, 'A_hat_', num2str(s), '_', num2str(b), '.png'], 'png')

        %a_block = squeeze(A_hat_block(b,s,:,:)) - diag(diag(squeeze(A_hat_block(b,s,:,:))));
        %a = squeeze(A_hat(s,:,:)) - diag(diag(squeeze(A_hat(s,:,:))));
        mat_dist(b,s) = norm(reshape(A_hat(s,:,:),[],1) - reshape(A_hat_block(b,s,:,:),[],1));

    end
end

boxplot(log(mat_dist)')
beta_block(beta_block > 1000) = 1000;
beta_block(beta_block < 0) = 0;
ahat_block_data = table([subjs,subjs]', [repmat(ones,nSubj,1);repmat(ones,nSubj,1)+1],...
    [mat_dist(1,:)'; mat_dist(2,:)'], [beta_block(1,:)'; beta_block(2,:)'], [beta;beta], ...
    'VariableNames', [{'subj'}, {'block'}, {'dist'}, {'beta_block'}, {'beta'}]);
writetable(ahat_block_data, [r_dir, 'ahat_block.csv']);

%% Get A_hat time varying

plt=0;
load([save_dir, 'max_ent.mat'], 'beta', 'A_hat')
Ahat_seq = zeros(nNode*nNode,nSubj);
Ahat_conv = cell(nSubj,1);
for s = 1:nSubj
    if ~exist([r_dir, 'subj', subjs{s}],'dir')
       mkdir([r_dir, 'subj', subjs{s}]); 
    end
    
    % select subject
    fprintf('\n************** Subj %s ************\n', subjs{s});
    load([ephys_dir, subjs{s}, '/task_data.mat'])
    curr_A = A_hat(s,:,:);
    
    %get ahat at every time point
    [curr_seq] = Ahat_sequence(walk+1, beta(s));
    
    %optional visualization
    if plt
    for i = 1:size(curr_seq,3)
        imagesc(curr_seq(:,:,i));
        title(num2str(i));
        pause(0.0001)
    end
    end
    
    %get similarity to final ahat
    curr = data(data.subj==str2double(subjs{s}),:);
    curr_v_seq = zeros(1,size(curr_seq,3));
    for i = 1:size(curr_seq,3)
        curr_v_seq(i) = corr(reshape(A_hat(s,:,:),[],1),reshape(curr_seq(:,:,i),[],1));
    end
    Ahat_conv{s} = [1:size(curr_seq,3);curr_v_seq];
end

figure(1); clf
all_trial = [];
all_corr = [];
all_subj = {}; all_beta = [];
cnt = 1;
for s = 1:nSubj
    if mod(str2double(subjs{s}),2) == 0
        color = 'k';
    else
        color = 'r';
    end
   plot(Ahat_conv{s}(1,:),Ahat_conv{s}(2,:), color); hold on 
   all_corr = [all_corr;Ahat_conv{s}(2,:)'];
   all_trial = [all_trial;Ahat_conv{s}(1,:)'];
   all_beta= [all_beta; repmat(beta(s), numel(Ahat_conv{s}(1,:)),1)];
   all_subj(cnt:(cnt + numel(Ahat_conv{s}(1,:))-1),1) = {subjs{s}};
   cnt = cnt + numel(Ahat_conv{s}(1,:));
end

ahat_seq_data = table(all_subj, all_trial, all_corr, all_beta, ...
    'VariableNames', [{'subj'}, {'trial'}, {'dist'}, {'beta'}]);
writetable(ahat_seq_data, [r_dir, 'ahat_seq.csv']);

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

