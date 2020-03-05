%% Calculate theta power for every trial
clear

addpath(genpath('/Users/stiso/Documents/MATLAB/IRASA/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
% define variables
subj = '6';
save_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_raw/';
r_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_analysis/';
img_dir = ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/subj', subj];

% make diractories
if ~exist(img_dir, 'dir')
    mkdir(img_dir);
end

% load stuff
load([save_dir, subj, '/ft_data.mat'])
load([save_dir, subj, '/header_clean.mat'], 'elec_labels', 'srate', 'HUP_ID', 'subj', 'AAL')
load([save_dir, subj, '/events.mat']) % in samples
load([save_dir, subj, '/task_data.mat'])
load([save_dir, subj, '/good_events.mat'])
load([r_dir, 'subj', subj, '/theta_peaks.mat'])

% analysis varaibles
nTrial = size(events,1);
nGoodTrial = size(good_events,1);
nElec = numel(elec_labels);

% freq variables
hg = [70,150];
lf_width = 2;

% get shortest rt
min_dur = min(cell2mat(cellfun(@(x) size(x,2), ft_data.trial, 'UniformOutput', false)));


%% Get HG

hg_ts = zeros(min_dur,nElec,nGoodTrial);

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = hg; % needs to be have a width around center frquency greater than twice the modulating low frequency (Aru et al 2019)
cfg.bpfiltdf = 0.5; % bandpass transition width
cfg.bsfiltdf = 0.5; % bandstop transition width
cfg.bpfiltdev = 0.01; % bandpass max passband deviation
cfg.bsfiltdev = 0.05; % bandstp max passband deviation
cfg.bpfilttype = 'firws'; % or 'firls' (slower), but avoid the default 'but' (= not well-suited for hilbert phase estimate)
cfg.hilbert = 'complex';

pow = ft_preprocessing(cfg, ft_data);

% get amp envelope
for t = 1:nGoodTrial
    hg_ts(:,:,t) = abs(pow.trial{t}(:,1:min_dur))';
end


%% Get theta for elecs with oscillations
% should this be phase?

%get one for each elec
theta_ts = zeros(min_dur,nElec,nGoodTrial);
for i = 1:nElec
    % get peak oscillation for this elec
    curr_peak = peaks(i);
    
    
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.channel = elec_labels{i};
    cfg.bpfreq = [curr_peak-lf_width curr_peak+lf_width];
    cfg.bpfiltdf = 0.5; % bandpass transition width
    cfg.bsfiltdf = 0.5; % bandstop transition width
    cfg.bpfiltdev = 0.01; % bandpass max passband deviation
    cfg.bsfiltdev = 0.05; % bandstp max passband deviation
    cfg.bpfilttype = 'firws'; % or 'firls' (slower), but avoid the default 'but' (= not well-suited for hilbert phase estimate)
    cfg.hilbert = 'complex';
    
    curr_theta = ft_preprocessing(cfg, ft_data);
    
    % get inst amp
    for t = 1:nGoodTrial
        theta_ts(:,i,t) = abs(curr_theta.trial{t}(:,1:min_dur));
    end
    
end

%% Get erp

raw_feats = zeros(size(hg_ts));
for t = 1:nGoodTrial
    raw_feats(:,:,t) = ft_data.trial{t}(:,1:min_dur)';
end

%% Prewhiten all trials - 
%only do this when you are using multiple elctrodes
% 
% raw_feats_white = zeros(size(raw_feats));
% 
% for i = 1:nElec
%     % row is obs, column is var
%     curr = squeeze(raw_feats(:,i,:))';
%     % ZCA whitening
%     mu = mean(curr);
%     X = bsxfun(@minus, curr, mu);
%     A = curr'*curr;
%     [V,D,~] = svd(A);
%     whMat = sqrt(size(curr,1)-1)*V*sqrtm(inv(D + eye(size(D))*0.0001))*V';
%     Xwh = X*whMat;
%     raw_feats_white(:,i,:) = Xwh';
% end


%% Normalize

pow_feats = [hg_ts; theta_ts];
%theta_feats = [theta_ts];

raw_feats = raw_feats - mean(raw_feats,3);
pow_feats = pow_feats - mean(pow_feats,3);
%theta_feats = theta_feats - mean(theta_feats,3);

%% Get cross-validated distance metric

% get each condition
good_walk = walk(good_trials);
nStim = 10;

% get corr distance
raw_RDA = zeros(nStim, nStim, nElec);
pow_RDA = zeros(nStim, nStim, nElec);

for k = 1:nElec
    curr_raw = squeeze(raw_feats(:,k,:));
    %curr_pow = squeeze(pow_feats(:,k,:));
    for i = 1:nStim
        stim_i = mean(curr_raw(:,(good_walk == (i-1))),2);
        stim_i = stim_i - mean(stim_i);
        for j = i:nStim
            stim_j = mean(curr_raw(:,(good_walk == (j-1))),2);
            stim_j = stim_j - mean(stim_j);
            % corr distance
            raw_RDA(i,j,k) = (stim_i)'*(stim_j)...
                /(norm(stim_i)*norm(stim_j));
        end
    end
    raw_RDA(:,:,k) = raw_RDA(:,:,k) + raw_RDA(:,:,k)';
end




%% PCA

raw_coeff = zeros(nStim, 2, nElec);
%pow_coeff = zeros(nStim, 2, nElec);
%theta_coeff = zeros(nStim, 2, nElec);

trans1 = [253,224,239]./255;
trans2 = [230,245,208]./255;
within1 = [233,163,201]./255;
within2 = [161,215,106]./255;
center1 = [197,27,125]./255;
center2 = [77,146,33]./255;
colors = [trans1; within1; center1; within1; trans1; trans2; within2; center2; within2; trans2];

for i = 1:nElec
    raw_coeff(:,:,i) = pca(raw_RDA(:,:,i), 'Centered', true, 'NumComponents', 2);
    %pow_coeff(:,:,i) = pca(pow_RDA(:,:,i), 'Centered', true, 'NumComponents', 2);
    %theta_coeff(:,:,i) = pca(theta_RDA(:,:,i), 'Centered', true, 'NumComponents', 2);
    
    figure(1); clf
    scatter(raw_coeff(:,1,i), raw_coeff(:,2,i), 10000, colors, '.', 'MarkerFaceAlpha', 0.4)
    title(AAL{i,1})
    pause(0.5)
    saveas(gca, [img_dir, '/PCA_raw_', elec_labels{i}, '.png'], 'png')
    
%     figure(2); clf
%     scatter(pow_coeff(:,1,i), pow_coeff(:,2,i), 10000, colors, '.', 'MarkerFaceAlpha', 0.4)
%     title(AAL{i,1})
%     saveas(gca, [img_dir, '/PCA_pow_', elec_labels{i}, '.png'], 'png')
%     
%     figure(3); clf
%     scatter(theta_coeff(:,1,i), theta_coeff(:,2,i), 10000, colors, '.', 'MarkerFaceAlpha', 0.4)
%     title(AAL{i,1})
%     saveas(gca, [img_dir, '/PCA_theta_', elec_labels{i}, '.png'], 'png')
end

