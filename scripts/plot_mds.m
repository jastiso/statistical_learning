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
sep_pca = zeros(size(subjs));
K = 500; % number random cluster selections
n = 10; %total number of nodes
null_sep = zeros(numel(subjs),K);
null_sep_pca = zeros(numel(subjs),K);
compressibility = zeros(numel(subjs),1);
losses = zeros(numel(subjs,1));

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
    ahat_idx = find(cellfun(@(x) strcmp(subj,x),A_hat_order));

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
    
    %% MDS 
    % neural data
    
    [Y, E] = cmdscale(D,2);
    
    figure(2); clf
    scatter(Y(:,1), Y(:,2), 10000, colors, '.', 'MarkerFaceAlpha', 0.4)
    title(subj)
    saveas(gca, [img_dir, '/MDS_', feat_type, '.png'], 'png')
    
    % now with pca
    coeff = pca(D);
    x = D - mean(D);
    y_pca = x*coeff;
    y_pca = y_pca(:,1:2);
    
    % behavioral data
    A_hat_dist = squeeze(A_hat(ahat_idx,:,:)) - diag(diag(squeeze(A_hat(ahat_idx,:,:))));
    [Yb, ~] = cmdscale(A_hat_dist,2);
    
    figure(3); clf
    scatter(Yb(:,1), Yb(:,2), 10000, colors, '.', 'MarkerFaceAlpha', 0.4)
    title(subj)
    saveas(gca, [img_dir, '/MDS_behave.png'], 'png')
    
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
        % discriminability
        ytab = table(Y(:,1), Y(:,2), 'VariableNames', [{'x1'}, {'x2'}]);
        fit = fitcdiscr(ytab, [{'1'},{'1'},{'1'},{'1'},{'1'},{'2'},{'2'},{'2'}, {'2'},{'2'}]);
        losses(s) = loss(fit, ytab, [{'1'},{'1'},{'1'},{'1'},{'1'},{'2'},{'2'},{'2'}, {'2'},{'2'}]);
        
        % discriminability
        ytab = table(Yb(:,1), Yb(:,2), 'VariableNames', [{'x1'}, {'x2'}]);
        fit = fitcdiscr(ytab, [{'1'},{'1'},{'1'},{'1'},{'1'},{'2'},{'2'},{'2'}, {'2'},{'2'}]);
        losses_b(s) = loss(fit, ytab, [{'1'},{'1'},{'1'},{'1'},{'1'},{'2'},{'2'},{'2'}, {'2'},{'2'}]);
        
        mod1 = 2:4;
        mod2 = 7:9;
        curr_dist = 0;
        curr_dist_pca = 0;
        for n = 1:nNode
            mod_flag = any(n == mod1);
            
            if mod_flag
                tmp_dist = module_sep(Y,mod1,n);
                curr_dist = curr_dist + tmp_dist;
                % pca
                tmp_dist = module_sep(y_pca,mod1,n);
                curr_dist_pca = curr_dist_pca + tmp_dist;
            elseif any(n == mod2)
                tmp_dist = module_sep(Y,mod2,n);
                curr_dist = curr_dist + tmp_dist;
                % pca
                tmp_dist = module_sep(y_pca,mod2,n);
                curr_dist_pca = curr_dist_pca + tmp_dist;
            end
        end
        sep(s) = curr_dist/numel([mod1,mod2]);
        sep_pca(s) = curr_dist_pca/numel([mod1,mod2]);
    else
        for u = 1:K
            nodes = datasample(1:nNode, 6, 'Replace', false);
            mod1 = nodes(1:3);
            mod2 = nodes(4:6);
            
            % check if these overlap
            if any(mod1 == mod2)
                warning('Your set selection process is wrong')
            else
                curr_dist = 0;
                curr_dist_pca = 0;
                for n = 1:nNode
                    mod_flag = any(n == mod1);
                    
                    if mod_flag
                        tmp_dist = module_sep(Y,mod1,n);
                        curr_dist = curr_dist + tmp_dist;
                        % pca
                        tmp_dist = module_sep(y_pca,mod1,n);
                        curr_dist_pca = curr_dist_pca + tmp_dist;
                    elseif any(n == mod2)
                        tmp_dist = module_sep(Y,mod2,n);
                        curr_dist = curr_dist + tmp_dist;
                        % pca
                        tmp_dist = module_sep(y_pca,mod2,n);
                        curr_dist_pca = curr_dist_pca + tmp_dist;
                    end
                end
                null_sep(s,u) = curr_dist/numel(nodes);
                null_sep_pca(s,u) = curr_dist_pca/numel(nodes);
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

figure(3); clf
scatter(log10(beta_mod), losses(mod_idx), '.');

[r,p] = corr(log(beta_mod), sep', 'type','spearman')

data = table(subjs(mod_idx)', sep', beta_mod, ...
    'VariableNames', [{'subj'}, {'module_dist'}, {'beta'}]);
writetable(data, [r_dir, 'mod_dist.csv']);

null_data = table(repmat(subjs(~mod_idx),1,K)', reshape(null_sep,[],1), repmat(beta_lat',1,K)',...
    repelem(1:K,sum(~mod_idx))',...
    'VariableNames', [{'subj'}, {'module_dist'}, {'beta'}, {'set'}]);
writetable(null_data, [r_dir, 'null_mod_dist.csv']);

comp_data = table(subjs', compressibility, beta, ...
    'VariableNames', [{'subj'}, {'compress'}, {'beta'}]);
writetable(comp_data, [r_dir, 'comp_data.csv']);

comp_data = table(subjs', losses, losses_b, beta, ...
    'VariableNames', [{'subj'}, {'losses'}, {'loss_beh'}, {'beta'}]);
writetable(comp_data, [r_dir, 'comp_data.csv']);
