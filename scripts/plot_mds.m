%% Plot MDS for all sig elecs

clear
clc

addpath(genpath('/Users/stiso/Documents/MATLAB/Colormaps/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
addpath(genpath('/Users/stiso/Documents/Code/graph_learning/functions/'))
save_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_raw/';
r_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_analysis/';


subjs = [{'1'}, {'2'}, {'3'}, {'4'}, {'8'}, {'10'}, {'12'}];
feat_type = 'lfp';
A_hat_order = [{'1'}, {'2'}, {'3'}, {'4'}, {'5'}, {'6'}, {'8'}, {'10'},{'12'},{'18'}];
load('/Users/stiso/Documents/Code/graph_learning/ECoG_data/behavior_preprocessed/max_ent.mat', 'A_hat', 'beta')
stresses = zeros(size(subjs));

trans1 = [253,224,239]./255;
trans2 = [230,245,208]./255;
within1 = [233,163,201]./255;
within2 = [161,215,106]./255;
center1 = [197,27,125]./255;
center2 = [77,146,33]./255;
cluster_colors = [trans1; within1; center1; within1; trans1; trans2; within2; center2; within2; trans2];


for s = 1:numel(subjs)
    subj = subjs{s};
    load([r_dir, 'subj' subj, '/searchlight_corrs.mat'], 'sig_idx', 'perm_corr')
    
    % load stuff
    load([save_dir, subj, '/ft_data.mat'], 'ft_data')
    load([save_dir, subj, '/header_clean.mat'], 'elec_labels', 'srate')
    load([save_dir, subj, '/good_events.mat'], 'good_events', 'good_trials') % in samples
    load([save_dir, subj, '/task_data.mat'], 'walk')
    
    % analysis varaibles
    nElec = size(ft_data.trial{1},1);
    nNode = numel(unique(walk));
    if mod(str2double(subj),2) == 0
        colors = cluster_colors;
    else
        colors = viridis(nNode);
    end
    img_dir = ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/subj', subj];

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
    
    if strcmp(feat_type, 'pow')
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
        
    elseif strcmp(feat_type, 'lfp')
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
    feats = feats(:,sig_idx,:);
    feats = reshape(feats, nTrial, sum(sig_idx)*size(feats,3), []);
    
    
    %% RDM
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
        [d,m] = get_rdm(feats, train, test, good_walk, nNode);
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
    
    [Y, E] = cmdscale(D,2);
    
    figure(2); clf
    scatter(Y(:,1), Y(:,2), 10000, colors, '.', 'MarkerFaceAlpha', 0.4)
    title(subj)
    saveas(gca, [img_dir, '/MDS_', feat_type, '.png'], 'png')
    
    
    %% Get stress
    
    %subj_idx = find(cellfun(@(x) strcmp(subj,x),A_hat_order));
    %A_hat = squeeze(A_hat(subj_idx,:,:));
    
    d_vect = nonzeros(triu(D,1));
    d_hat = squareform(pdist(Y));
    d_hat = nonzeros(triu(d_hat,1));
    stress = sqrt((sum(d_vect - d_hat)^2)/sum(d_vect.^2));
    mod_dist = pdist([mean(Y(1:5,1)),mean(Y(1:5,2));mean(Y(6:end,1)),mean(Y(6:end,2))]);
    
    stresses(s) = mod_dist;
    
end

beta = beta(cellfun(@(x) any(strcmp(x,subjs)), A_hat_order));
mod_idx = cellfun(@(x) mod(str2double(x),2) == 0,subjs);
stresses = stresses(mod_idx);
beta = beta(mod_idx');
stresses = stresses(beta ~= 0);
beta = nonzeros(beta);
figure(3); clf
scatter(stresses,log(beta))

[r,p] = corr(stresses', log(beta))