%% IRASA pipeline

% Identify electrodes with theta range electrodes
% repeat for different windows, with and without lowpass filtering

addpath(genpath('/Users/stiso/Documents/MATLAB/IRASA/'))
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current/'))
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
load([save_dir, subj, '/events.mat'])
load([save_dir, subj, '/trans_idx.mat'])

% analysis varaibles
freqs = 0:2:srate/2;
nTrial = size(events,1);
nElec = size(data,1);
nFreq = numel(freqs);

% timing varaibles
exp_st = events(1);
exp_en = events(end);

%% IRASA
% separates scale free from oscillatory component of the power spectra.

win_length = 3000; % in ms
step = 1000; %in ms
filter = 1; % with ot without lowpass alaising filter

spec = cell(nElec,1); % initialize
for j = 1:nElec
    
    
    spec{j} = get_IRASA_spec(data(j,:), exp_st, exp_en, srate, win_length, step, filter);
    
    % plot
    % show averaged fractal and oscillatory power spectrum
    figure(1); clf;
    subplot(2,1,1);
    loglog(spec{j}.freq,mean(spec{j}.mixd,2),'b',  'linewidth', 2); hold on
    loglog(spec{j}.freq,mean(spec{j}.frac,2),'r', 'linewidth', 2);
    title(AAL{j});
    subplot(2,1,2);
    plot(spec{j}.freq, mean(spec{j}.osci,2), 'linewidth', 2); hold on
    shade_plot(spec{j}.freq', mean(spec{j}.osci,2)', (std(spec{j}.osci,[],2))', rgb("slategrey"), 0.4);
    saveas(gca, [img_dir, '/IRASA_', elec_labels{j}, '.png'], 'png')
end

save([r_dir, 'IRASA.mat'], 'spec')


%% IRASA- #no filter

win_length = 3000; % in ms
step = 1000; %in ms
filter = 0; % with ot without lowpass alaising filter

spec_nf = cell(nElec,1); % initialize
for j = 1:nElec
    fprintf('Elec %s...\n', elec_labels{j})
    spec_nf{j} = get_IRASA_spec(data(j,:), exp_st, exp_en, srate, win_length, step, filter);
    
end

save([r_dir, 'IRASA_no_filter.mat'], 'spec_nf')

%% IRASA- short window

win_length = 1000; % in ms
step = 100; %in ms
filter = 1; % with ot without lowpass alaising filter

spec_short = cell(nElec,1); % initialize
for j = 1:nElec
    fprintf('Elec %s...\n', elec_labels{j})
    spec_short{j} = get_IRASA_spec(data(j,:), exp_st, exp_en, srate, win_length, step, filter);
    
end

save([r_dir, 'IRASA_short.mat'], 'spec_short')

%% Stats
% want to test if the peak in the theta range is greater than 1/f, 
% how to correct across electrodes?
% time windows are not independent


theta = [4,10];
alpha = 0.05/nElec;

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
    
    
end
sum(ps)
%ps(ps > alpha) = 1;
save([r_dir, 'theta_peaks.mat'], 'peaks', 'ps')

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
    
    %perm test with 1/f and psd
    scale_free = curr.frac(peak_log,:);
    psd = curr.mixd(peak_log,:);
    % is the peak greater than 3stds in scale free component?
    ps_nf(j) = (mean(scale_free) + 3*std(scale_free)) < mean(curr.mixd(peak_log,:));
    % how do I correct this across elecs?
    
    
end

sum(ps_nf)
sum(ps == ps_nf)/nElec
%ps(ps > alpha) = 1;

%% Stats - small window
% want to test if the peak in the theta range is greater than 1/f, 
% how to correct across electrodes?
% time windows are not independent


ps_short = zeros(nElec,1);
for j = 1:nElec
    curr = spec_short{j};
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
    ps_short(j) = (mean(scale_free) + 3*std(scale_free)) < mean(curr.mixd(peak_log,:));
    % how do I correct this across elecs?
    
    
end

sum(ps_short)
sum(ps == ps_short)/nElec
%ps(ps > alpha) = 1;
sum(ps & ps_short)