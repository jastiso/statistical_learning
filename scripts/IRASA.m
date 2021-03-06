%% IRASA pipeline
clear

% Identify electrodes with theta range electrodes
% repeat for different windows, with and without lowpass filtering

addpath(genpath('/Users/stiso/Documents/MATLAB/IRASA/'))
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current/'))
removepath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')) % this hides the built in "hann" function

% define variables
subj = '8';
save_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_raw/';
r_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_analysis/';
img_dir = ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/subj', subj];

% make diractories
if ~exist(img_dir, 'dir')
    mkdir(img_dir);
end

% load stuff
load([save_dir, subj, '/data_clean.mat'])
load([save_dir, subj, '/header_clean.mat'], 'elec_labels', 'srate', 'HUP_ID', 'subj', 'regions', 'sessions')
load([save_dir, subj, '/events.mat'])
%load([save_dir, subj, '/trans_idx.mat'])
if exist('data_all','var')
    data = data_all;
else
    tmp = data;
    data = [];
    data(1).sess = tmp;
end

% analysis varaibles
freqs = 0:2:srate/2;
nTrial = size(events,1);
nElec = numel(elec_labels);
nFreq = numel(freqs);
nSess = numel(sessions);

% timing varaibles
exp_st = zeros(nSess,1);
exp_en = zeros(nSess,1);

for i = 1:nSess
    try
        curr = events(events(:,3) == i,:);
        exp_st(i) = curr(1,1);
        exp_en(i) = curr(end,2);
    catch % some subjects don't have session in the events
        exp_st(i) = events(1,1);
        exp_en(i) = events(end,2);
    end
end

%freqs variables
theta = [3,12];


%% IRASA
% separates scale free from oscillatory component of the power spectra.

win_length = 3000; % in ms
step = 1000; %in ms
filter = 1; % with ot without lowpass alaising filter

spec = cell(nElec,nSess); % initialize
for s = 1:numel(sessions)
    for j = 1:nElec
        
        
        spec{j,s} = get_IRASA_spec(data(s).sess(j,:), exp_st(s), exp_en(s), srate, win_length, step, filter);
        
        % plot
        % show averaged fractal and oscillatory power spectrum
        figure(1); clf;
        subplot(2,1,1);
        loglog(spec{j}.freq,mean(spec{j}.mixd,2),'b',  'linewidth', 2); hold on
        loglog(spec{j}.freq,mean(spec{j}.frac,2),'r', 'linewidth', 2);
        title(regions{j});
        subplot(2,1,2);
        plot(spec{j}.freq, mean(spec{j}.osci,2), 'linewidth', 2); hold on
        shade_plot(spec{j}.freq', mean(spec{j}.osci,2)', (std(spec{j}.osci,[],2))', rgb("slategrey"), 0.4);
        saveas(gca, [img_dir, '/IRASA_sess', num2str(s), '_', elec_labels{j}, '.png'], 'png')
        
    end
end

save([r_dir, 'subj', subj, '/IRASA.mat'], 'spec')


%% IRASA- #no filter

win_length = 3000; % in ms
step = 1000; %in ms
filter = 0; % with ot without lowpass alaising filter

spec_nf = cell(nElec,nSess); % initialize
for s = 1:nSess
    for j = 1:nElec
        fprintf('Elec %s...\n', elec_labels{j})
        spec_nf{j,s} = get_IRASA_spec(data(s).sess(j,:), exp_st(s), exp_en(s), srate, win_length, step, filter);
        
    end
end
save([r_dir, 'subj', subj, '/IRASA_no_filter.mat'], 'spec_nf')


%% Stats
% want to test if the peak in the theta range is greater than 1/f,
% how to correct across electrodes?
% time windows are not independent


alpha = 0.05/nElec;

peak_power = zeros(nElec,1);
peaks = zeros(nElec,1);
ps = zeros(nElec,1);
for j = 1:nElec
    curr = spec{j};
    % find peak in theta range
    theta_idx = curr.freq >= theta(1) & curr.freq <= theta(2);
    osci = mean(curr.osci,2);
    [~, peak_idx] = max(osci(theta_idx));
    tmp = curr.freq(theta_idx);
    peaks(j) = tmp(peak_idx);
    peak_log = curr.freq == peaks(j);
    
    %perm test with 1/f and psd
    scale_free = curr.frac(peak_log,:);
    psd = curr.mixd(peak_log,:);
    % is the peak greater than 3stds in scale free component?
    ps(j) = (mean(scale_free) + 3*std(scale_free)) < mean(curr.mixd(peak_log,:));
    % how do I correct this across elecs?
    
    % get distribution of values at peak, and see if mean or median makes
    % more sense - looks like median
    hist(curr.osci(peak_log,:));
    pause(0.1)
    % get the peak of the oscillation
    peak_power(j) = median(curr.osci(peak_log,:));
end
sum(ps)

hist(peak_power,30)

save([r_dir, 'subj', subj, '/theta_peaks.mat'], 'peaks', 'ps')

%% Stats - no filter
% want to test if the peak in the theta range is greater than 1/f,
% how to correct across electrodes?
% time windows are not independent


ps_nf = zeros(nElec,1);
for j = 1:nElec
    curr = spec_nf{j};
    % find peak in theta range
    theta_idx = curr.freq >= theta(1) & curr.freq <= theta(2);
    osci = mean(curr.osci,2);
    [~, peak_idx] = max(osci(theta_idx));
    tmp = curr.freq(theta_idx);
    peaks(j) = tmp(peak_idx);
    peak_log = curr.freq == peaks(j);
    
    % test with 1/f and psd
    scale_free = curr.frac(peak_log,:);
    psd = curr.mixd(peak_log,:);
    % is the peak greater than 3stds in scale free component?
    ps_nf(j) = (mean(scale_free) + 3*std(scale_free)) < mean(curr.mixd(peak_log,:));
    % how do I correct this across elecs?
    
    
end

sum(ps_nf)
sum(ps == ps_nf)/nElec
%ps(ps > alpha) = 1;


%% Explore window length

s = 1; % just the first session, the important variable should be subject, not day
wins = [500, 1000, 2000, 3000, 4000, 5000, 7000]; % in ms
filter = 1; % with or without lowpass alaising filter

mean_scale_free = zeros(nElec, numel(wins));
std_scale_free = zeros(nElec, numel(wins));
mean_psd = zeros(nElec, numel(wins));

for i = 1:numel(wins)
    win_length = wins(i);
    step = win_length/2;
    for j = 1:nElec
        fprintf('Elec %s...\n', elec_labels{j})
        curr = get_IRASA_spec(data(s).sess(j,:), exp_st(s), exp_en(s), srate, win_length, step, filter);
        
        % get relevant stats, to see which is changing
        theta_idx = curr.freq >= theta(1) & curr.freq <= theta(2);
        osci = mean(curr.osci,2);
        [~, peak_idx] = max(osci(theta_idx));
        tmp = curr.freq(theta_idx);
        peaks(j) = tmp(peak_idx);
        peak_log = curr.freq == peaks(j);
        
        mean_scale_free(j,i) = mean(curr.frac(peak_log,:));
        std_scale_free(j,i) = std(curr.frac(peak_log,:));
        mean_psd(j,i) = mean(curr.mixd(peak_log,:));
    end
end

figure(1); clf
subplot(3,1,1)
plot(wins,mean_scale_free', 'linewidth', 2); title('Mean 1/f')
xlim([500,7000])
subplot(3,1,2)
plot(wins,std_scale_free', 'linewidth', 2); title('Std 1/f')
xlim([500,7000])
subplot(3,1,3)
plot(wins,mean_psd', 'linewidth', 2); title('Mean PSD')
xlim([500,7000])
saveas(gca, [img_dir, 'win_size.png'], 'png')
