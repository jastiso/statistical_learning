%% Compare theta power at within vs between cluster trantisions

clear
clc

addpath(genpath('/Users/stiso/Documents/MATLAB/IRASA/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
addpath(genpath('/Users/stiso/Documents/Code/graph_learning/functions/'))

% define variables
subjs = [{'1'}, {'2'}, {'3'}, {'4'}, {'6'}, {'8'}, {'10'}];
feature = 'lfp'; % pow or lfp
freqs = logspace(log10(3), log10(150), 50);

trans1 = [253,224,239]./255;
trans2 = [230,245,208]./255;
within1 = [233,163,201]./255;
within2 = [161,215,106]./255;
center1 = [197,27,125]./255;
center2 = [77,146,33]./255;
colors = [trans1; within1; center1; within1; trans1; trans2; within2; center2; within2; trans2];


for subj_idx = 1:numel(subjs)
    clear D_perm
    subj = subjs{subj_idx};
    save_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_raw/';
    r_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_analysis/';
    img_dir = ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/subj', subj];
    
    % make diractories
    if ~exist(img_dir, 'dir')
        mkdir(img_dir);
    end
    
    % load stuff
    load([save_dir, subj, '/ft_data.mat'])
    load([save_dir, subj, '/header_clean.mat'], 'elec_labels', 'srate', 'HUP_ID', 'subj')
    load([save_dir, subj, '/good_events.mat']) % in samples
    
    load([save_dir, subj, '/task_data.mat'])
    
    % analysis varaibles
    nElec = size(ft_data.trial{1},1);
    nNode = numel(unique(walk));
    tri_mask = logical(triu(ones(nNode),1));
    
    % timing varaibles
    exp_st = good_events(1);
    exp_en = good_events(end);
    
    % benchmarks
    
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
    if mod(str2double(subjs{subj_idx}),2) == 0
        A = M;
    else
        A = L;
    end
    
    G = expm(A);
    % load A_hat
    load('/Users/stiso/Documents/Code/graph_learning/ECoG_data/behavior_preprocessed/max_ent.mat', 'A_hat')
    A_hat = squeeze(A_hat(subj_idx,:,:));
    
    
    %% Field trip format
    
    % some useful variables
    % need at least 200 ms for power analyses, also anything short doesnt
    % really make sense
    min_dur = min(good_events(:,2) - good_events(:,1));
    good_walk = walk(good_trials) + 1; % sswitch to matlab indexing
    keep_idx = true(numel(good_walk),1);
    if min_dur < .2*srate
        keep_idx = keep_idx & (good_events(:,2) - good_events(:,1) > .2*srate);
    end
    good_walk = good_walk(keep_idx);
    good_events = good_events(keep_idx,:);
    
    % prewhiten
    cfg = [];
    cfg.derivative = 'yes';
    ft_data = ft_preprocessing(cfg, ft_data);
    nTrial = sum(keep_idx);
    
    if strcmp(feature, 'pow')
        % get frequency band feats
        cfg = [];
        cfg.method = 'mtmfft';
        cfg.channel = elec_labels;
        cfg.taper = 'dpss';
        cfg.ouput = 'pow';
        cfg.pad = 'nextpow2';
        cfg.trials = keep_idx;
        cfg.foi = freqs; % this includes notch filtered freqs within the band!!
        cfg.keeptrials = 'yes';
        cfg.tapsmofrq = 4; %smoothing index - check if same effect is present for others
        
        pow = ft_freqanalysis(cfg,ft_data);
        feats = log10(pow.powspctrm);
        
    elseif strcmp(feature, 'lfp')
        % cut to same number of timepoints
        cfg = [];
        cfg.trl = [good_events(:,1), good_events(:,1) + min_dur, zeros(size(good_events,1),1)];
        % if multiple sessions, add that
        if isfield(ft_data, 'trialinfo')
            cfg.trl = [cfg.trl, ft_data.trialinfo];
        end
        ft_data = ft_redefinetrial(cfg,ft_data);
        % reshape into Trial x (timepoint*elec)
        feats = zeros(nTrial, nElec, size(ft_data.trial{1},2));
        for i = 1:nTrial
            feats(i,:,:) = ft_data.trial{i};
        end
    end
    
    %% get CV normalized euclidean distance
    % same as Mahalanobis Distance
    
    G_corr = zeros(nElec,1);
    A_corr = zeros(nElec,1);
    A_hat_corr = zeros(nElec,1);
    for e = 1:nElec
        curr_feats = squeeze(feats(:,e,:));
        elec = ft_data.label{e};
        
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
            [d,m] = get_rdm(curr_feats, train, test, good_walk, nNode);
            D = D + d;
            N = N + m;
        end
        D = D./N;
        D(logical(eye(nNode))) = NaN;
        
        % average dist
        figure(1); clf
        imagesc(D); colorbar
        D(logical(eye(nNode))) = 0;
        
        %% MDS for each channel
        
        [Y, stress] = cmdscale(D,2);
        
        figure(2); clf
        scatter(Y(:,1), Y(:,2), 10000, colors, '.', 'MarkerFaceAlpha', 0.4)
        title(elec)
        saveas(gca, [img_dir, '/MDS_', feature, '_', elec, '.png'], 'png')
        
        % correlations
        G_corr = corr(reshape(G(tri_mask),[],1), reshape(D(tri_mask),[],1));
        A_corr = corr(reshape(A(tri_mask),[],1), reshape(D(tri_mask),[],1));
        A_hat_corr = corr(reshape(A_hat(tri_mask),[],1), reshape(D(tri_mask),[],1));
        
    end
end



%     %% Null Models
%
%     % permutation
%     nPerm = 100;
%     D_perm = zeros(nNode, nNode, nPerm);
%     for n = 1:nPerm
%         fprintf('\nPerm %d', n);
%         % count the number of times you see each transition
%         N = zeros(nNode);
%         r = randi([2,nTrial-1]);
%         perm_walk = [good_walk(r:end), good_walk(1:(r-1))];
%         % leave one out cv
%         k = nTrial;
%         for i = 1:k
%             % split
%             train = true(nTrial,1);
%             train(i) = false;
%             test = ~train;
%
%             % get dist
%             [d,m] = get_rdm(feats, train, test, perm_walk, nNode);
%             D_perm(:,:,n) = D_perm(:,:,n) + d;
%             N = N + m;
%         end
%         D_perm(:,:,n) = D_perm(:,:,n)./N;
%         figure(1); clf;
%         imagesc(D_perm(:,:,n))
%     end
%     % D_perm(logical(eye(nNode))) = 0;
%
%     % correlations
%     G_corr_perm = zeros(nPerm,1);
%     A_corr_perm = zeros(nPerm,1);
%     hat_corr_perm = zeros(nPerm,1);
%
%     for n = 1:nPerm
%         curr = D_perm(:,:,n);
%         curr(logical(eye(nNode))) = 0;
%         G_corr_perm(n) = corr(reshape(G(tri_mask),[],1), reshape(curr(tri_mask),[],1));
%         A_corr_perm(n) = corr(reshape(A(tri_mask),[],1), reshape(curr(tri_mask),[],1));
%         hat_corr_perm(n) = corr(reshape(A_hat(tri_mask),[],1), reshape(curr(tri_mask),[],1));% space
%     end
%
%     figure(1); clf
%     histogram(G_corr_perm, 'Normalization', 'probability', 'FaceColor',rgb('steelblue'),'EdgeColor','white','facealpha', 0.8); hold on
%     plot([G_corr, G_corr], [0, .3], 'red', 'linewidth', 3)
%     title('Communicability')
%     saveas(gca, [img_dir, '/G_rsa_', feature, '.png'], 'png')
%
%     figure(2); clf
%     histogram(A_corr_perm, 'Normalization', 'probability', 'FaceColor',rgb('steelblue'),'EdgeColor','white','facealpha', 0.8); hold on
%     plot([A_corr, A_corr], [0, .3], 'red', 'linewidth', 3)
%     title('A')
%     saveas(gca, [img_dir, '/A_rsa_', feature, '.png'], 'png')
%
%     figure(3); clf
%     histogram(hat_corr_perm, 'Normalization', 'probability', 'FaceColor',rgb('steelblue'),'EdgeColor','white','facealpha', 0.8); hold on
%     plot([A_hat_corr, A_hat_corr], [0, .3], 'red', 'linewidth', 3)
%     title('A_hat')
%     saveas(gca, [img_dir, '/A_hat_rsa_', feature, '.png'], 'png')
%
% end