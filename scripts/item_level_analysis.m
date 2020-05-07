%% Compare theta power at within vs between cluster trantisions

clear 
clc

addpath(genpath('/Users/stiso/Documents/MATLAB/IRASA/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
addpath(genpath('/Users/stiso/Documents/Code/graph_learning/functions/'))

% define variables
subj = '4';
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
load([save_dir, subj, '/events.mat']) % in samples
load([save_dir, subj, '/task_data.mat'])
load([r_dir, 'subj', subj, '/theta_peaks.mat'])

% analysis varaibles
nElec = size(ft_data.trial{1},1);
hg = [70:150];
nNode = numel(unique(walk));

% timing varaibles
exp_st = events(1);
exp_en = events(end);


%% get extra event variables

%get only good trials
good_trials = logical(correct) & cutoff;
good_events = events(logical(correct) & cutoff,:);


%% Field trip format

% cut out bad trials
min_dur = min(good_events(:,2) - good_events(:,1));
cfg = [];
cfg.trl = [good_events, zeros(size(good_events,1),1)];

ft_data = ft_redefinetrial(cfg,ft_data);
good_walk = walk(good_trials) + 1; % sswitch to matlab indexing

% prewhiten
cfg = [];
cfg.derivative = 'yes';
ft_data = ft_preprocessing(cfg, ft_data);
nTrial = numel(ft_data.trial);

% get frequency bany feats
cfg = [];
cfg.method = 'mtmfft'; 
cfg.channel = elec_labels;
cfg.taper = 'dpss';
cfg.ouput = 'pow';
cfg.pad = 'nextpow2';
cfg.foi = hg; % this includes two notch filtered freqs within the band!!
cfg.keeptrials = 'yes';
cfg.tapsmofrq = 4; %smoothing index - check if same effect is present for others

pow = ft_freqanalysis(cfg,ft_data);
ft_hg = squeeze(mean(log10(pow.powspctrm),3));

cfg = [];
cfg.method = 'mtmfft'; % only single taper for low freqs
cfg.channel = elec_labels;
cfg.taper = 'hanning';
cfg.ouput = 'pow';
cfg.pad = 'nextpow2';
cfg.foi = unique(round(peaks));
cfg.keeptrials = 'yes';

pow = ft_freqanalysis(cfg,ft_data);
ft_theta = zeros(size(ft_hg));
for i = 1:nElec
    curr_peak = round(peaks(i));
    peak_idx = find(curr_peak == pow.freq');
    ft_theta(:,i) = squeeze(log10(pow.powspctrm(:,i,peak_idx)));
end
ft_spec = [ft_theta, ft_hg];

% cut to same number of timepoints
cfg = [];
cfg.trl = [good_events(:,1), good_events(:,1) + min_dur, zeros(size(good_events,1),1)];

ft_data = ft_redefinetrial(cfg,ft_data);
% reshape into Trial x (timepoint*elec)
feats = zeros(nTrial, numel(ft_data.trial{1}));
for i = 1:nTrial
   feats(i,:) = reshape(ft_data.trial{i},1,[]);
end

%% get CV Mahalanobis Distance

D = zeros(nNode);
% leave one out cv
k = nTrial;

for i = 1:k
    N = zeros(nNode);
    % split
    train = true(nTrial,1);
    train(i) = false;
    test = ~train;
    
  % get dist
  curr_node = good_walk(test);
  for j = 1:nNode
      if j ~= curr_node
          X = ft_theta(train,:);
          X = X(good_walk(train) == j,:);
          Y = ft_theta(test,:);
          d = mahal(Y,X);
          D(curr_node,j) = D(curr_node,j) + d;
          D(j,curr_node) = D(j,curr_node) + d;
          N(curr_node,j) = N(curr_node,j) + 1;
          N(j,curr_node) = N(j,curr_node) + 1;
      end
  end
  D./N;
end
D = D./k;
D(logical(eye(nNode))) = 0;

% average dist
figure(1); clf
imagesc(D); colorbar

%% MDS

[Y, stress] = mdscale(D,2);

trans1 = [253,224,239]./255;
trans2 = [230,245,208]./255;
within1 = [233,163,201]./255;
within2 = [161,215,106]./255;
center1 = [197,27,125]./255;
center2 = [77,146,33]./255;
colors = [trans1; within1; center1; within1; trans1; trans2; within2; center2; within2; trans2];

figure(1); clf
scatter(Y(:,1), Y(:,2), 10000, colors, '.', 'MarkerFaceAlpha', 0.4)
    

%% Communcability

A = [0 1 1 1 0 0 0 0 0 1;
    1 0 1 1 1 0 0 0 0 0;
    1 1 0 1 1 0 0 0 0 0;
    1 1 1 0 1 0 0 0 0 0;
    0 1 1 1 0 1 0 0 0 0;
    0 0 0 0 1 0 1 1 1 0;
    0 0 0 0 0 1 0 1 1 1;
    0 0 0 0 0 1 1 0 1 1;
    0 0 0 0 0 1 1 1 0 1;
    1 0 0 0 0 0 1 1 1 0];
G = expm(A);
