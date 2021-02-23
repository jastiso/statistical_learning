%% Plot MDS for all sig elecs

clear
clc

addpath(genpath('/Users/stiso/Documents/MATLAB/Colormaps/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
addpath(genpath('/Users/stiso/Documents/Code/graph_learning/functions/'))
save_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_raw/';
r_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_analysis/';


subjs = [{'1'}, {'2'}, {'3'}, {'4'}, {'5'}, {'6'}, {'7'}, {'8'}, {'10'}, {'12'}];
feat_type = 'lfp';
A_hat_order = load([r_dir, 'ahat_order.mat']);
A_hat_order = A_hat_order.subjs;
load('/Users/stiso/Documents/Code/graph_learning/ECoG_data/behavior_preprocessed/max_ent.mat', 'A_hat', 'beta')
sep = zeros(size(subjs));
K = 5; % number of items in each fake cluster
n = 10; %total number of nodes
null_sep = zeros(numel(subjs),nchoosek(n,K)/2);
compressibility = zeros(numel(subjs),1);

trans1 = [239,169,186]./255;
trans2 = [177,191,146]./255;
within1 = [174,116,133]./255;
within2 = [126,138,96]./255;
cluster_colors = [trans1; within1; within1; within1; trans1; trans2; within2; within2; within2; trans2];
mid_blue = [108,111,147]./255;
dark_blue = [87,89,117]./255;
light_blue = [196,198,212]./255;
ring_colors = [[linspace(0, dark_blue(1), 5)'; linspace(mid_blue(1), light_blue(1), 5)'],...
    [  linspace(0, dark_blue(2), 5)'; linspace(mid_blue(2), light_blue(2), 5)'],...
    [ linspace(0, dark_blue(3), 5)'; linspace(mid_blue(3), light_blue(3), 5)']];

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
        colors = ring_colors;
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
    
    
    %% Get compression/separation stats
    
    [rd_upper, rd_lower, Cs, ~] = rate_distortion(D, 1, 100);
    compressibility(s) = mean(rd_upper(end) - rd_upper);
    
    %subj_idx = find(cellfun(@(x) strcmp(subj,x),A_hat_order));
    %A_hat = squeeze(A_hat(subj_idx,:,:));
    
    d_vect = nonzeros(triu(D,1));
    d_hat = squareform(pdist(Y));
    d_hat = nonzeros(triu(d_hat,1));
    stress = sqrt((sum(d_vect - d_hat)^2)/sum(d_vect.^2));
    if mod(str2double(subj),2) == 0
        mod1 = 1:5;
        mod2 = 6:10;
        %num = mean([pdist(Y(mod1,:)), pdist(Y(mod2,:))]);
        denom = 0;
        num = 0;
        for n = 1:nNode
            mod_flag = any(n == mod1);
            
            if mod_flag
                num = num + mean(sqrt((Y(mod1,1) - Y(n,1)).^2 + (Y(mod1,2) - Y(n,2)).^2));
                denom = denom + mean(sqrt((Y(mod2,1) - Y(n,1)).^2 + (Y(mod2,2) - Y(n,2)).^2));
            else
                num = num + mean(sqrt((Y(mod2,1) - Y(n,1)).^2 + (Y(mod2,2) - Y(n,2)).^2));
                denom = denom + mean(sqrt((Y(mod1,1) - Y(n,1)).^2 + (Y(mod1,2) - Y(n,2)).^2));
            end
        end
        denom = denom/nNode;
        sep(s) = num - denom;
    else
        C = nchoosek(1:10,K);
        for u = 1:(size(C,1)/2)
            mod1 = C(u,:);
            mod2 = C((end-(u-1)),:);
            
            % check if these overlap
            if any(mod1 == mod2)
                warning('Your set election process is wrong')
            else
                %num = mean([pdist(Y(mod1,:)), pdist(Y(mod2,:))]);
                denom = 0;
                num = 0;
                for n = 1:nNode
                    mod_flag = any(n == mod1);
                    if mod_flag
                        num = num + mean(sqrt((Y(mod1,1) - Y(n,1)).^2 + (Y(mod1,2) - Y(n,2)).^2));
                        denom = denom + mean(sqrt((Y(mod2,1) - Y(n,1)).^2 + (Y(mod2,2) - Y(n,2)).^2));
                    else
                        num = num + mean(sqrt((Y(mod2,1) - Y(n,1)).^2 + (Y(mod2,2) - Y(n,2)).^2));
                        denom = denom + mean(sqrt((Y(mod1,1) - Y(n,1)).^2 + (Y(mod1,2) - Y(n,2)).^2));
                    end
                end
                denom = denom/nNode;
                num = num/nNode;
                null_sep(s,u) = num - denom;
            end
        end
    end
    
end

%% Plot

beta = beta(cellfun(@(x) any(strcmp(x,subjs)), A_hat_order));
%subjs = subjs(mod_idx);
%sep = sep(mod_idx);
sep = sep(beta ~= 0);
subjs = subjs(beta ~= 0);
compressiblity = compressibility(beta~=0);
beta = beta(beta~=0);
mod_idx = cellfun(@(x) mod(str2double(x),2) == 0,subjs);

beta_mod = beta(mod_idx');
sep = sep(mod_idx');
beta_lat = beta(~mod_idx);
null_sep = null_sep(~mod_idx,:);
figure(3); clf
scatter(log10(beta_mod), sep, '.'); hold on;
for j = 1:size(null_sep,2)
    plot(log10(beta_lat), null_sep(:,j), 'k');
end
%saveas(gca, ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/mod_sep.pdf'], 'pdf')

[r,p] = corr(log(beta_mod), sep', 'type','pearson')

data = table(subjs(mod_idx)', sep', beta_mod, ...
    'VariableNames', [{'subj'}, {'module_dist'}, {'beta'}]);
writetable(data, [r_dir, 'mod_dist.csv']);

null_data = table(repmat(subjs(~mod_idx),1,nchoosek(n,K)/2)', reshape(null_sep,[],1), repmat(beta_lat',1,nchoosek(n,K)/2)',...
    repelem(1:nchoosek(nNode,K)/2,sum(~mod_idx))',...
    'VariableNames', [{'subj'}, {'module_dist'}, {'beta'}, {'set'}]);
writetable(null_data, [r_dir, 'null_mod_dist.csv']);

comp_data = table(subjs', compressibility, beta, ...
    'VariableNames', [{'subj'}, {'compress'}, {'beta'}]);
writetable(comp_data, [r_dir, 'comp_data.csv']);
