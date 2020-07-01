%% Compare theta power at within vs between cluster trantisions

clear
clc

addpath(genpath('/Users/stiso/Documents/MATLAB/IRASA/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
addpath(genpath('/Users/stiso/Documents/MATLAB/BrainNetViewer_20171031/'))
addpath(genpath('/Users/stiso/Documents/Code/graph_learning/functions/'))

% define variables
subjs = [{'1'}, {'2'}, {'3'}, {'4'}, {'6'}, {'8'}, {'10'}];
feature = 'lfp'; % pow or lfp
freqs = logspace(log10(3), log10(150), 50);
nBlock = 4;
block_indices = [1, 250; 251, 500; 501, 750; 751, 1000];
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


trans1 = [253,224,239]./255;
trans2 = [230,245,208]./255;
within1 = [233,163,201]./255;
within2 = [161,215,106]./255;
center1 = [197,27,125]./255;
center2 = [77,146,33]./255;
cluster_colors = [trans1; within1; center1; within1; trans1; trans2; within2; center2; within2; trans2];


for subj_idx = 1:numel(subjs)
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
    load([save_dir, subj, '/header_clean.mat'])
    load([save_dir, subj, '/good_events.mat']) % in samples
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
    A_hat = squeeze(A_hat(subj_idx,:,:));
    
    % initialize data structure
    G_corr = nan(nElec,nBlock);
    A_corr = nan(nElec,nBlock);
    A_hat_corr = nan(nElec,nBlock);
            
    curr_trials = trial(good_trials);
    for b = 1:nBlock
        
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
        
        % get trial indicies for current block
        ind = block_indices(b,:);
        curr_idx = ((curr_trials > ind(1)) & (curr_trials <= ind(2)))' & keep_idx;
        
        if sum(curr_idx) == 0
           continue 
        end
        %% get features
        
        % prewhiten
        cfg = [];
        cfg.derivative = 'yes';
        curr_data = ft_preprocessing(cfg, ft_data);
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
            
            pow = ft_freqanalysis(cfg,curr_data);
            feats = log10(pow.powspctrm);
            
        elseif strcmp(feature, 'lfp')
            % cut to same number of timepoints
            cfg = [];
            cfg.trl = [good_events(:,1), good_events(:,1) + min_dur, zeros(size(good_events,1),1)];
            % if multiple sessions, add that
            if isfield(curr_data, 'trialinfo')
                cfg.trl = [cfg.trl, curr_data.trialinfo];
            end
            cfg.trl = cfg.trl(curr_idx,:);
            curr_data = ft_redefinetrial(cfg,curr_data);
            % reshape into Trial x (timepoint*elec)
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
            
            %% correlations
            
            % correlations
            G_corr(e,b) = corr(reshape(G(tri_mask),[],1), reshape(D(tri_mask),[],1));
            A_corr(e,b) = corr(reshape(A(tri_mask),[],1), reshape(D(tri_mask),[],1));
            A_hat_corr(e,b) = corr(reshape(A_hat(tri_mask),[],1), reshape(D(tri_mask),[],1));
            
        end
                        
    end
    save([r_dir, 'subj' subj, '/searchlight_corrs_block.mat'], 'G_corr', 'A_corr', 'A_hat_corr')
    
    %plot
    figure(1); clf
    plot(A_corr')
    title('A')
    saveas(gca, [img_dir, '/block_A.png'], 'png')
    
    figure(2); clf
    plot(A_hat_corr');
    title('A_hat')
    saveas(gca, [img_dir, '/block_A_hat.png'], 'png')
    
    figure(3); clf
    plot(G_corr');
    title('G')
    saveas(gca, [img_dir, '/block_G.png'], 'png')
end