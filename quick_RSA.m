%% Compare theta power at within vs between cluster trantisions

addpath(genpath('/Users/stiso/Documents/MATLAB/IRASA/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
% define variables
HUP_ID = 'HUP187';
subj = '2';
save_dir = '/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_raw/';
r_dir = '/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_analysis/';
img_dir = ['/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_img/subj', subj];

% make diractories
if ~exist(img_dir, 'dir')
    mkdir(img_dir);
end

% load stuff
load([save_dir, subj, '/data_clean.mat'])
load([save_dir, subj, '/header_clean.mat'], 'elec_labels', 'srate', 'HUP_ID', 'subj', 'AAL')
load([save_dir, subj, '/events.mat']) % in samples
load([save_dir, subj, '/task_data.mat'])
load([r_dir, 'subj', subj, '/theta_peaks.mat'])

% analysis varaibles
nTrial = size(events,1);
nElec = size(data,1);

% timing varaibles
exp_st = events(1);
exp_en = events(end);


%% get extra event variables

%get only good trials
good_trials = logical(correct) & cutoff;
good_events = events(logical(correct) & cutoff,:);


%% Field trip format

% get ft data
ft_data = fieldtrip_format(data, srate, elec_labels, []);


%% Get bandpasses high_frequency activity

cfg = [];
cfg.hpfilter = 'yes';
cfg.hpfreq = 70;
cfg.bsfilter = 'yes';
cfg.continuous = 'yes';
cfg.bsfreq = [118, 122; 178, 182; 238, 242];

ft_data = ft_preprocessing(cfg,ft_data);

cfg = [];
cfg.trl = [good_events, zeros(size(good_events,1),1)];

ft_data = ft_redefinetrial(cfg,ft_data);

%% get avg envelope

good_walk = walk(good_trials);
len = min(good_events(:,2) - good_events(:,1));
envelopes = zeros(numel(unique(good_walk)), len, nElec);


for i = 1:numel(good_walk)
    for j = 1:nElec
        curr = good_walk(i) + 1;
        envelopes(curr, :, j) = envelopes(curr, :, j) + ft_data.trial{i}(j,1:len).^2;
    end
end

% get avg
for i = 1:numel(unique(good_walk))
    N = numel(find(good_walk == i-1));
    envelopes(i,:,:) = envelopes(i,:,:)./N;
end

% corr plot
mat = zeros(numel(unique(good_walk)), numel(unique(good_walk)), nElec);
for i = 1:nElec
    curr = envelopes(:,:,i);
    mat(:,:,i) = corr(curr');
    mat(:,:,i) = mat(:,:,i) - eye(numel(unique(good_walk)));
    
    figure(1); clf;
    imagesc(mat(:,:,i)); c = colorbar; %caxis([-1,1])
    title(AAL{i})
    xlabel('Node'); ylabel('Node')
    c.Label.String = 'Pearsons R';
    saveas(gca, [img_dir, '/HG_rsa_', elec_labels{i}, '.png'], 'png')
end

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
nNode = size(G,1);

% get corrs
rsa_p = zeros(nElec,1);

for i = 1:nElec
   G_vect = G(triu(true(nNode),1));
   A = mat(:,:,i);
   vect = A(triu(true(nNode),1));
   [r,p] = corr(G_vect,vect);
   rsa_p(i) = p;
end

sum(rsa_p < 0.05)