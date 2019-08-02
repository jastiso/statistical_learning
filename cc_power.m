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

trans_nodes = [4,5,9,0];
alt_nodes = [0,2,7,5];

%get only good trials
good_events = events(logical(correct) & cutoff,:);

trans_idx = false(nTrial,1);
alt_idx = false(nTrial,1);
for i = 2:nTrial
    if any(trans_nodes == walk(i-1)) && any(trans_nodes == walk(i))
        trans_idx(i) = true;
    end
    if any(alt_nodes == walk(i-1)) && any(alt_nodes == walk(i))
        alt_idx(i) = true;
    end
end

% match to events
trans_idx = trans_idx(logical(correct) & cutoff);
alt_idx = alt_idx(logical(correct) & cutoff);
cc_events = good_events(trans_idx,:);
alt_events = good_events(alt_idx,:);

%% Field trip format

% only elecs with signigifant peak
data = data(logical(ps),:);
elec_labels = elec_labels(logical(ps),:);
peaks = peaks((logical(ps)));
AAL = AAL(logical(ps),:);
nElec = size(data,1);

% cross cluster data
ft_data = [];
ft_data.trial{1,1} = data;
ft_data.time{1,1} = (1/srate):(1/srate):(size(data,2)/srate);
ft_data.label = elec_labels;
ft_data.fsample = srate;
ft_data.nSamples = size(data,2);

cfg = [];
cfg.trl = [good_events, zeros(size(good_events,1),1)];
ft_data = ft_redefinetrial(cfg, ft_data);

nGoodTrial = size(good_events,1);


%% Get power at peak index

cfg = [];
cfg.method = 'mtmfft';
cfg.ouput = 'pow';
cfg.pad = 'nextpow2';
cfg.foi = unique(round(peaks));
cfg.keeptrials = 'yes';
cfg.tapsmofrq = 4; %smoothing index - check if same effect is present for others

pow = ft_freqanalysis(cfg,ft_data);


%% Compare

cc_ps = zeros(nElec,1);
corr_ps = zeros(nElec,1);
means =zeros(nElec,1);
peak_power = zeros(nElec, nGoodTrial);
for i = 1:nElec
    curr_peak = round(peaks(i));
    peak_idx = find(curr_peak == pow.freq');
    
    figure(1); clf
    [N,X] = hist(log10(pow.powspctrm(trans_idx,i,peak_idx)),30);
    bar(X./sum(trans_idx), N, 0.8, 'facecolor', rgb('blue'), 'facealpha', .5); hold on
    [N2,X2] = hist(log10(pow.powspctrm(alt_idx,i,peak_idx)),30);
    bar(X2./sum(alt_idx), N2, 0.8, 'facecolor', rgb('orange'), 'facealpha', .5); hold on
    legend([{'Between Cluster'},{'Within Cluster'}])
    xlabel('log(power)')
    %pause(0.5)
    
    [p] = permtest(log10(pow.powspctrm(alt_idx,i,peak_idx)), log10(pow.powspctrm(trans_idx,i,peak_idx)), 10000, 'conservative');
    means(i) = mean(log10(pow.powspctrm(alt_idx,i,peak_idx))) - log10(mean(pow.powspctrm(trans_idx,i,peak_idx)));
    cc_ps(i) = p;
    
    
    [r,pcorr] = corr(log10(pow.powspctrm(:,i,peak_idx)),good_events(:,2) - good_events(:,1));
    corr_ps(i) = pcorr;
    figure(2); clf
    scatter(log10(pow.powspctrm(:,i,peak_idx)),good_events(:,2) - good_events(:,1));
    ylabel('log(power)')
    xlabel('RT')
    pause(0.01);
    
    %add to peak data
    peak_power(i,:) = pow.powspctrm(:,i,peak_idx);
end

sig = cc_ps < 0.05/nElec;
sum(cc_ps < 0.05/nElec)

cc_ps(sig)
means(sig)
AAL(sig,:)

% get data for python anaylsis
save([r_dir, 'peak_power.mat'], 'peak_power', 'elec_labels', 'AAL')
