%% Compare theta power at within vs between cluster trantisions

addpath(genpath('/Users/stiso/Documents/MATLAB/IRASA/'))
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
alt_nodes = [0,1,6,5];
module0 = [0 1 2 3 4];
module1 = [5 6 7 8 9];

%get only good trials
good_trials = logical(correct) & cutoff;
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

% constrast for increasign activity within a module
module_idx = zeros(nTrial,1);
prev_mod = any(find(walk(1) == module1)); % flag for which module we're in
for i = 2:nTrial
    curr_mod = any(find(walk(i) == module1));
    if prev_mod == curr_mod
       module_idx(i) = module_idx(i-1) + 1; 
    else
       module_idx(i) = 0; 
    end
    prev_mod = curr_mod;
end

% match to events
trans_idx = trans_idx(logical(correct) & cutoff);
module_idx = module_idx(logical(correct) & cutoff);
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

% firldtrip data
ft_data = fieldtrip_format(data, srate, elec_labels, [good_events, zeros(size(good_events,1),1)]);

nGoodTrial = size(good_events,1);

% demean and detrend data
for i = 1:numel(ft_data.trial)
    ft_data.trial{i} = ft_data.trial{i} - mean(ft_data.trial{i},2);
    ft_data.trial{i} = detrend(ft_data.trial{i});
end


%% Get power at peak index

cfg = [];
cfg.method = 'mtmfft'; % only single taper for low freqs
cfg.taper = 'hanning';
cfg.ouput = 'pow';
cfg.pad = 'nextpow2';
cfg.foi = unique(round(peaks));
cfg.keeptrials = 'yes';
%cfg.tapsmofrq = 4; %smoothing index - check if same effect is present for others

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
    
    [p] = permtest(log10(pow.powspctrm(trans_idx,i,peak_idx)), log10(pow.powspctrm(alt_idx,i,peak_idx)), 10000, 'conservative');
    means(i) = mean(log10(pow.powspctrm(trans_idx,i,peak_idx))) - log10(mean(pow.powspctrm(alt_idx,i,peak_idx)));
    cc_ps(i) = p;
    
    
    [r,pcorr] = corr(log10(pow.powspctrm(:,i,peak_idx)),good_events(:,2) - good_events(:,1));
    corr_ps(i) = pcorr;
    figure(2); clf
    scatter(log10(pow.powspctrm(:,i,peak_idx)),good_events(:,2) - good_events(:,1));
    ylabel('log(power)')
    xlabel('RT')
    pause(0.01);
    
    %add to peak data
    peak_power(i,:) = log10(pow.powspctrm(:,i,peak_idx));
end

sig = cc_ps < 0.05;
sum(cc_ps < 0.05/nElec)

cc_ps(sig)
means(sig)
AAL(sig,:)
find(sig)
%things for r
region = AAL(1,:);
good_trial_idx = trial(good_trials);

% get data for python anaylsis
save([r_dir, 'subj', subj, '/peak_power.mat'], 'peak_power', 'elec_labels', 'region', 'good_trial_idx', 'module_idx')

%% Plots

nNode = numel(unique(walk));
RT = zeros(nNode);
power_node = zeros(nNode, nNode, nElec);

for i = 1:nElec
    
    count_node = zeros(nNode);
    
    cnt = 0; % not zero because j counter starts at 2
    for j = 2:numel(walk)
        prev = walk(j-1) + 1;
        curr = walk(j) + 1;
        if good_trials(j)
            cnt = cnt + 1;
            power_node(prev,curr,i) = power_node(prev,curr,i) + peak_power(i,cnt);
            count_node(prev,curr) = count_node(prev,curr) + 1;
            
            if i == nElec % only need to do this once
                RT(prev,curr) = RT(prev,curr) + events(cnt,2) - events(cnt,1);
            end
        end
    end
    power_node(:,:,i) = power_node(:,:,i)./count_node;
    
    %plot
    figure(3); clf
    imagesc(power_node(:,:,i)); c = colorbar;
    title(AAL{i})
    xlabel('current node')
    ylabel('previous node')
    c.Label.String = 'power';
    saveas(gca, [img_dir, '/pow_by_trans_', elec_labels{i}, '.png'], 'png')
    
end

RT = RT./count_node;
figure(4); clf
imagesc(RT); c = colorbar;
xlabel('current node')
ylabel('previous node')
c.Label.String = 'RT';
saveas(gca, [img_dir, '/rt_by_trans.png'], 'png')


