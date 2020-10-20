function [] = searchlight_rsa(subj, subj_idx, cluster_colors, nSim, feat_type)
% performs a searchlight RSA analysis on each electrode for a given
% subject. 
%   this is mostly to help with parallelization

save_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_raw/';
r_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_analysis/';
img_dir = ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/subj', subj];

% make diractories
if ~exist(img_dir, 'dir')
    mkdir(img_dir);
end

% load stuff
load([save_dir, subj, '/ft_data.mat'], 'ft_data')
load([save_dir, subj, '/header_clean.mat'], 'elec_labels', 'srate')
load([save_dir, subj, '/good_events.mat'], 'good_events', 'good_trials') % in samples
try
    fname = dir([save_dir, subj, '/*/electrodenames_coordinates_mni.csv']);
    coords = readtable([fname.folder, '/', fname.name]);
    coords.Var1 = fix_names(coords.Var1);
    coord_idx = ismember(cell2mat(coords.Var1),cell2mat(elec_labels),'rows');
    coords = coords(coord_idx,:);
    flag = true;
catch
    fprintf('No region data')
    flag = false;
end
load([save_dir, subj, '/task_data.mat'], 'walk')

% analysis varaibles
nElec = size(ft_data.trial{1},1);
nNode = numel(unique(walk));
tri_mask = logical(triu(ones(nNode),1));

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
if mod(str2double(subj),2) == 0
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

%% get CV normalized euclidean distance
% same as Mahalanobis Distance

G_corr = zeros(nElec,1);
A_corr = zeros(nElec,1);
A_hat_corr = zeros(nElec,1);
sig_idx = zeros(nElec,1);
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
    
    [Y, ~] = cmdscale(D,2);
    
    figure(2); clf
    scatter(Y(:,1), Y(:,2), 10000, colors, '.', 'MarkerFaceAlpha', 0.4)
    title(elec)
    saveas(gca, [img_dir, '/MDS_', feat_type, '_', elec, '.png'], 'png')
    
    % correlations
    G_corr(e) = corr(reshape(G(tri_mask),[],1), reshape(D(tri_mask),[],1));
    A_corr(e) = corr(reshape(A(tri_mask),[],1), reshape(D(tri_mask),[],1));
    A_hat_corr(e) = corr(reshape(A_hat(tri_mask),[],1), reshape(D(tri_mask),[],1));
    
    % test against permutation model
    perm_corr = zeros(nSim,1);
    for n = 1:nSim
       % pick random split
       split = 1;
       while (split == 1) || (split == nTrial)
           split = randi(nTrial);
       end
       % split walk
       perm_walk = [good_walk(split:end), good_walk(1:(split-1))];
       
       % rerun cross validated distance measure
       D_perm = zeros(nNode);
       N_perm = zeros(nNode);
       
       for i = 1:k
           
           % split
           train = true(nTrial,1);
           train(i) = false;
           test = ~train;
           
           % get dist
           [d,m] = get_rdm(curr_feats, train, test, perm_walk, nNode);
           D_perm = D_perm + d;
           N_perm = N_perm + m;
       end
       D_perm = D_perm./N_perm;
       D_perm(logical(eye(nNode))) = NaN;
       perm_corr(n) = corr(reshape(A_hat(tri_mask),[],1), reshape(D_perm(tri_mask),[],1));
    end
    % test if empirical corr is greater than 95% percent of nulls
    fprintf('\nFor subject %s, contact %s is more correlated than %d null models', subj, elec, sum(A_hat_corr(e) < perm_corr))
    figure(1); clf
    histogram(perm_corr); hold on
    plot([A_hat_corr(e), A_hat_corr(e)], [0,20], 'r')
    % update index of elecs
    sig_idx(e) = sum(A_hat_corr(e) < perm_corr) >= (.95*nSim);
end

save([r_dir, 'subj' subj, '/searchlight_corrs.mat'], 'G_corr', 'A_corr', 'A_hat_corr', 'sig_idx')

% save node file
sig_idx = logical(sig_idx);
if flag
    subj_node_file = [r_dir, 'subj', subj, '/A_hat_', feat_type, '.node'];
    write_bv_node( subj_node_file, coords.Var2(sig_idx), coords.Var3(sig_idx), coords.Var4(sig_idx),...
        A_hat_corr(sig_idx), [], elec_labels(sig_idx));
    
    subj_node_file = [r_dir, 'subj', subj, '/A_', feat_type, '.node'];
    write_bv_node( subj_node_file, coords.Var2(sig_idx), coords.Var3(sig_idx), coords.Var4(sig_idx)...
        , A_corr(sig_idx), [], elec_labels(sig_idx));
    
    subj_node_file = [r_dir, 'subj', subj, '/G_', feat_type, '.node'];
    write_bv_node( subj_node_file, coords.Var2(sig_idx), coords.Var3(sig_idx), coords.Var4(sig_idx)...
        , G_corr(sig_idx), [], elec_labels(sig_idx));
    
    % plot
    BrainNet_MapCfg('/Users/stiso/Documents/MATLAB/BrainNetViewer_20171031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv',...
        [r_dir, 'subj', subj, '/A_hat_', feat_type, '.node'],[r_dir, 'rsa_corr.mat'], ...
        ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/subj', subj, '/A_hat_brain_', feat_type, '.jpg']);
    
    
    BrainNet_MapCfg('/Users/stiso/Documents/MATLAB/BrainNetViewer_20171031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv',...
        [r_dir, 'subj', subj, '/A_', feat_type, '.node'],[r_dir, 'rsa_corr.mat'], ...
        ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/subj', subj, '/A_brain_', feat_type, '.jpg']);
    
    
    BrainNet_MapCfg('/Users/stiso/Documents/MATLAB/BrainNetViewer_20171031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv',...
        [r_dir, 'subj', subj, '/G_', feat_type, '.node'],[r_dir, 'rsa_corr.mat'], ...
        ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/subj', subj, '/G_brain_', feat_type, '.jpg']);
else
    [~,I] = max(A_corr);
    fprintf('Max A correlation at %s\n', elec_labels{I})
    
    [~,I] = max(G_corr);
    fprintf('Max G correlation at %s\n', elec_labels{I})
    
    [~,I] = max(A_hat_corr);
    fprintf('Max A_hat correlation at %s\n', elec_labels{I})
end
end

