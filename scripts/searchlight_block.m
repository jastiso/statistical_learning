%% Get similarity matrices for smaller blocks of trials

clear
clc

addpath(genpath('/Users/stiso/Documents/MATLAB/IRASA/'))
addpath(genpath('/Users/stiso/Documents/MATLAB/Colormaps/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
addpath(genpath('/Users/stiso/Documents/MATLAB/BrainNetViewer_20171031/'))
addpath(genpath('/Users/stiso/Documents/Code/graph_learning/functions/'))

% define variables
subjs = [ {'1'},{'2'},{'3'},{'4'},{'5'},{'6'},{'7'},{'8'},{'10'},{'12'}];
A_hat_order = load([r_dir, 'ahat_order.mat']);
A_hat_order = A_hat_order.subjs;
feature = 'lfp'; % pow or lfp
freqs = logspace(log10(3), log10(150), 50);
shift = 100; % how much to slide windows
win = 500; % number of trials per window
nBlock = (1000 - win)/shift + 1;
block_indices = zeros(nBlock,2);
st = 1;
for b = 1:nBlock
    block_indices(b,:) = [st,st+(win)-1];
    st = st + shift;
end
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
% load null model
load('/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_analysis/RSA_dist_null.mat');
D_null = D;
clear D

trans1 = [239,169,186]./255;
trans2 = [177,191,146]./255;
within1 = [174,116,133]./255;
within2 = [126,138,96]./255;
cluster_colors = [trans1; within1; within1; within1; trans1; trans2; within2; within2; within2; trans2];

% initialize for stats
varnames = {'subj','space','block','elec','corr', 'full_corr'};
block_data = cell2table(cell(0,6), 'VariableNames',varnames);

for subj_idx = 1:numel(subjs)
    subj = subjs{subj_idx};
    ahat_idx = find(cellfun(@(x) strcmp(subj,x),A_hat_order));
    save_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_raw/';
    r_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_analysis/';
    img_dir = ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/subj', subj];
    
    
    % make diractories
    if ~exist(img_dir, 'dir')
        mkdir(img_dir);
    end
    
    % load stuff
    load([save_dir, subj, '/ft_data.mat'])
    load([save_dir, subj, '/header_clean.mat'])
    load([save_dir, subj, '/good_events.mat']) % in samples
    load([r_dir, 'subj' subj, '/searchlight_corrs.mat'],'sig_idx', 'A_hat_corr')
    full_ahat_corr = A_hat_corr;
    load([save_dir, subj, '/order.mat'],'order')
    D_null = D_null(order,order);
    
    try
        fname = dir([save_dir, subj, '/*/electrodenames_coordinates_mni.csv']);
        coords = readtable([fname.folder, '/', fname.name]);
        coords.Var1 = fix_names(coords.Var1);
        coord_idx = ismember(cell2mat(coords.Var1),cell2mat(elec_labels),'rows');
        coords = coords(coord_idx,:);
        flag = true;
    catch
        fprintf('No task data')
        flag = false;
    end
    load([save_dir, subj, '/task_data.mat'])
    
    % analysis varaibles
    nElec = size(ft_data.trial{1},1);
    nNode = numel(unique(walk));
    tri_mask = logical(triu(ones(nNode),1));
    
    % timing varaibles
    exp_st = good_events(1);
    exp_en = good_events(end);
    
    % benchmarks
    if mod(str2double(subjs{subj_idx}),2) == 0
        A = M;
        colors = cluster_colors;
    else
        A = L;
        colors = viridis(nNode);
    end
    G = expm(A);
    % load A_hat
    load('/Users/stiso/Documents/Code/graph_learning/ECoG_data/behavior_preprocessed/max_ent.mat', 'A_hat')
    A_hat = squeeze(A_hat(ahat_idx,:,:));
    
    % initialize data structure
    G_corr = nan(nElec,nBlock);
    A_corr = nan(nElec,nBlock);
    A_hat_corr = nan(nElec,nBlock);
    N_corr = nan(nElec,nBlock);
    curr_trials = trial(good_trials);
    % need at least 200 ms for power analyses, also anything short doesnt
    % really make sense
    min_dur = min(good_events(:,2) - good_events(:,1));
    good_walk = walk(good_trials) + 1; % sswitch to matlab indexing
    keep_idx = true(numel(good_walk),1);
    if min_dur < .2*srate
        keep_idx = keep_idx & (good_events(:,2) - good_events(:,1) > .2*srate);
    end
    
    
    % prewhiten
    cfg = [];
    cfg.derivative = 'yes';
    ft_data = ft_preprocessing(cfg, ft_data);
    
    for b = 1:nBlock
        clear feats curr_data
        
        % get trial indicies for current block
        ind = block_indices(b,:);
        curr_idx = ((curr_trials > ind(1)) & (curr_trials <= ind(2)))' & keep_idx;
        curr_walk = good_walk(curr_idx);
        
        if sum(curr_idx) == 0
            continue
        end
        %% get features
        nTrial = sum(curr_idx);
        
        if strcmp(feature, 'pow')
            % get frequency band feats
            cfg = [];
            cfg.method = 'mtmfft';
            cfg.channel = elec_labels;
            cfg.taper = 'dpss';
            cfg.ouput = 'pow';
            cfg.pad = 'nextpow2';
            cfg.trials = curr_idx;
            cfg.foi = freqs; % this includes notch filtered freqs within the band!!
            cfg.keeptrials = 'yes';
            cfg.tapsmofrq = 4; %smoothing index - check if same effect is present for others
            
            pow = ft_freqanalysis(cfg,ft_data);
            feats = log10(pow.powspctrm);
            curr_data = ft_data;
            
        elseif strcmp(feature, 'lfp')
            % cut to same number of timepoints
            cfg = [];
            cfg.trl = [good_events(:,1), good_events(:,1) + min_dur, zeros(size(good_events,1),1)];
            % if multiple sessions, add that
            if isfield(ft_data, 'trialinfo')
                cfg.trl = [cfg.trl, ft_data.trialinfo];
            end
            cfg.trl = cfg.trl(curr_idx,:);
            curr_data = ft_redefinetrial(cfg,ft_data);
            % reshape into Trial x timepoint x elec
            feats = zeros(nTrial, nElec, size(curr_data.trial{1},2));
            for i = 1:nTrial
                feats(i,:,:) = curr_data.trial{i};
            end
        end
        
        %% get CV normalized euclidean distance
        % same as Mahalanobis Distance
        
        
        for e = 1:nElec
            curr_feats = squeeze(feats(:,e,:));
            elec = curr_data.label{e};
            
            D = zeros(nNode);
            N = zeros(nNode);
            % leave one out cv
            k = nTrial;
            
            for i = 1:k
                
                % split
                train = true(nTrial,1);
                train(i) = false;
                test = ~train;
                
                % get dist
                [d,m] = get_rdm(curr_feats, train, test, curr_walk, nNode);
                D = D + d;
                N = N + m;
            end
            D = D./N;
            D(logical(eye(nNode))) = NaN;
            
            %% correlations
            
            % correlations
            G_corr(e,b) = -corr(reshape(G(tri_mask),[],1), reshape(D(tri_mask),[],1));
            A_corr(e,b) = -corr(reshape(A(tri_mask),[],1), reshape(D(tri_mask),[],1));
            A_hat_corr(e,b) = -corr(reshape(A_hat(tri_mask),[],1), reshape(D(tri_mask),[],1));
            N_corr(e,b) = corr(reshape(D_null(tri_mask), [], 1), reshape(D(tri_mask),[],1));
            
            
        end
        
        %% MDS
        if b == 1
            feats = feats(:,sig_idx,:);
            feats = reshape(feats, nTrial, sum(sig_idx)*size(feats,3), []);
            
            D = zeros(nNode);
            N = zeros(nNode);
            % leave one out cv
            k = nTrial;
            
            for i = 1:k
                
                % split
                train = true(nTrial,1);
                train(i) = false;
                test = ~train;
                
                % get dist
                [d,m] = get_rdm(feats, train, test, curr_walk, nNode);
                D = D + d;
                N = N + m;
            end
            D = D./N;
            D(logical(eye(nNode))) = 0;
            [Y, ~] = cmdscale(D,2);
            
            
            figure(2); clf
            scatter(Y(:,1), Y(:,2), 10000, colors, '.', 'MarkerFaceAlpha', 0.4)
            title([subj, ' block', num2str(b)])
            saveas(gca, [img_dir, '/MDS_', feature, '_block', num2str(b), '.png'], 'png')
        end
        
    end
    A_corr = A_corr(sig_idx,:);
    A_hat_corr = A_hat_corr(sig_idx,:);
    N_corr = N_corr(sig_idx,:);
    full_ahat_corr = full_ahat_corr(sig_idx);
    
    save([r_dir, 'subj' subj, '/searchlight_corrs_block.mat'], 'G_corr', 'A_corr', 'A_hat_corr')
    num = (sum(sig_idx)*nBlock*2);
    block_idx = [];
    for m = 1:nBlock
        block_idx = [block_idx; ones(sum(sig_idx),1)+(m-1)];
    end
    block_idx = [block_idx; block_idx];
    block_data = [block_data; table(repmat({subj}, num,1), ...
        [repmat({'latent'}, sum(sig_idx)*nBlock,1); repmat({'euclid'}, sum(sig_idx)*nBlock,1)],...
        num2cell(block_idx),...
        repmat(elec_labels(sig_idx,1), nBlock*2,1),num2cell([reshape(A_hat_corr,[],1); reshape(N_corr,[],1)]),...
        repmat(full_ahat_corr, nBlock*2,1),...
        'VariableNames', varnames)];
    
    %plot
    figure(1); clf
    plot(A_corr', 'linewidth', 2); hold on
    plot(1:nBlock, mean(N_corr), 'k', 'linewidth', 3)
    shade_plot(1:nBlock, mean(N_corr), std(N_corr, [], 1)/sqrt(nElec), 'k', 0.4)
    title('A')
    saveas(gca, [img_dir, '/block_A.png'], 'png')
    
    figure(2); clf
    plot(mean(A_hat_corr)', 'y', 'linewidth', 2); hold on
    shade_plot(1:nBlock, mean(A_hat_corr), std(A_hat_corr, [], 1)/sqrt(sum(sig_idx)),rgb('gold'), 0.4)
    plot([1:nBlock]', mean(N_corr)', 'k', 'linewidth', 3)
    shade_plot(1:nBlock, mean(N_corr), std(N_corr, [], 1)/sqrt(sum(sig_idx)), 'k', 0.4)
    title('A_hat')
    saveas(gca, [img_dir, '/block_A_hat.png'], 'png')
    %
    %     figure(3); clf
    %     plot(G_corr', 'linewidth', 2); hold on
    %     plot(1:nBlock, mean(N_corr), 'k', 'linewidth', 3)
    %     shade_plot(1:nBlock, mean(N_corr), std(N_corr, [], 1)/sqrt(nElec), 'k', 0.4)
    %     title('G')
    %     saveas(gca, [img_dir, '/block_G.png'], 'png')
end

writetable(block_data,[r_dir, 'block_searchlight.csv'])