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

%win_length = floor((min(rt) + 50)*(srate/1000)) - 1;
win_length = floor((3000)*(srate/1000)) - 1; % in samples
step = floor(1000*(srate/1000));
nWin = floor(((exp_en-exp_st)*(srate/1000) - win_length)/step); % (length - win)/step

spec = cell(nElec,1); % initialize
for j = 1:nElec
sig = zeros(win_length + 1, nWin);
for i = 1:nWin
    onset = events(i,1);
    %sig(:,i) = data(j,onset:(onset + win_length));
    st = ceil((exp_st*(srate/1000)) + (i-1)*step + 1);
    en = ceil(st + win_length);
    
    sig(:,i) = data(j, st:en);
end

% spec = amri_sig_fractal(sig,srate,...); optional arguments frange,
% detrend, filter
spec{j} = amri_sig_fractal(sig,srate, 'frange', [2, 70], 'detrend', 1, 'filter', 1);

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
    