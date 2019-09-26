%% Phase locking

% see if electrodes are showing phase locking to stim onset across all
% trials
clear

addpath(genpath('/Users/stiso/Documents/MATLAB/IRASA/'))
addpath(genpath('/Users/stiso/Documents/Code/graph_learning/'))
addpath(genpath('/Users/stiso/Documents/MATLAB/CircStat2012a/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
% define variables
subj = '4';
save_dir = '/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_raw/';
r_dir = '/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_analysis/';
img_dir = ['/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_img/subj', subj];

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

% freq variables
theta = [4, 12];

%% Channel Selection

% only elecs with signigifant peak
elec_labels = elec_labels(logical(ps),:);
peaks = peaks((logical(ps)));
AAL = AAL(logical(ps),:);
nElec = size(elec_labels,1);

%% Get fourier spectrum via wavelet

% wavelet centered on every time point, from the shortest RT
min_rt = min(good_events(:,2) - good_events(:,1));

cfg = [];
cfg.bpfilter = 'yes';
cfg.channel = elec_labels;
cfg.bpfreq = theta;
cfg.bpfiltdf = 0.5; % bandpass transition width
cfg.bsfiltdf = 0.5; % bandstop transition width
cfg.bpfiltdev = 0.01; % bandpass max passband deviation
cfg.bsfiltdev = 0.05; % bandstp max passband deviation
cfg.bpfilttype = 'firws'; % or 'firls' (slower), but avoid the default 'but' (= not well-suited for hilbert phase estimate)
cfg.hilbert = 'complex';

wave = ft_preprocessing(cfg, ft_data);


%% Get inst phase

phi = zeros(nElec, min_rt, nGoodTrial);
for t = 1:nGoodTrial
    phi(:,:,t) = (atan2(imag(wave.trial{t}(:,1:min_rt)),real(wave.trial{t}(:,1:min_rt))));
end


%% test uniformity for duration of trial
% at some point...we want to say they are less uniform than thier RTs
% get index of rayleight stat for every point in time, for every elec

ray_stat = zeros(nElec,min_rt);
ray_pval = zeros(nElec,min_rt);

for i = 1:nElec
    
    curr = squeeze(phi(:,i,:));
    for t = 1:min_rt
       [ray_stat(i,t), ray_pval(i,t)] = rayleigh_stat(curr(:,t));
    end
    
    figure(1); clf
    subplot(1,2,1)
    plot(ray_stat(i,:), 'linewidth', 2); hold on
    title('rayleigh stat')
    subplot(1,2,2)
    plot(ray_pval(i,:), 'linewidth', 2); hold on
    title('p-val')
    pause(0.1)
end


%% Plot raster-like heatmap for each electrode



%% Get plot of uniformity for each transition

plot_data = ft_data;
% get step and impulse
step = [ones(1,2000), ones(1,2000)*2];
impulse = ones(1,4000);
impulse(2000) = 2;
% get new data strucutre
plot_data.trial = {[impulse; step]};
plot_data.label = [{'impulse'}, {'step'}];
plot_data.sampleinfo = [1*srate, 3000*srate];
plot_data.time = {0.001:0.001:4};
% change some things for plotting
cfg.plotfiltresp = 'yes';
cfg.hilbert = 'no';
cfg.channel = plot_data.label;
filter_resp = ft_preprocessing(cfg, plot_data);
saveas(gca, [img_dir, '/ft_filter_resp.png'], 'png')


figure(1); clf
subplot(2,2,1)
plot(filter_resp.trial{1}(1,:)', 'linewidth', 2); hold on
title('Impulse Response')
xlabel('Time')
ylabel('Amplitude')
subplot(2,2,2)
plot(plot_data.trial{1}(1,:)', 'k', 'linewidth', 2);
title('Impulse')
xlabel('Time')
ylabel('Amplitude')
subplot(2,2,3)
plot(filter_resp.trial{1}(2,:)', 'linewidth', 2); hold on
title('Step Response')
xlabel('Time')
ylabel('Amplitude')
subplot(2,2,4)
plot(plot_data.trial{1}(2,:)', 'k', 'linewidth', 2)
title('Step')
xlabel('Time')
ylabel('Amplitude')
saveas(gca, [img_dir, '/filter_resp.png'], 'png')