%% Calculate theta power for every trial
clear

addpath(genpath('/Users/stiso/Documents/MATLAB/IRASA/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
% define variables
subj = '10';
save_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_raw/';
r_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_analysis/';
img_dir = ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/subj', subj];

% make diractories
if ~exist(img_dir, 'dir')
    mkdir(img_dir);
end

% load stuff
load([save_dir, subj, '/ft_data.mat'])
load([save_dir, subj, '/header_clean.mat'], 'elec_labels', 'srate', 'HUP_ID', 'subj', 'regions')
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
regions = regions(logical(ps),:);
nElec = size(elec_labels,1);


%% Get power at peak index

cfg = [];
cfg.method = 'mtmfft'; % only single taper for low freqs
cfg.channel = elec_labels;
cfg.taper = 'hanning';
cfg.ouput = 'pow';
cfg.pad = 'nextpow2';
cfg.foi = unique(round(peaks));
cfg.keeptrials = 'yes';
%cfg.tapsmofrq = 4; %smoothing index - check if same effect is present for others

pow = ft_freqanalysis(cfg,ft_data);


%% Get power at peaks only


corr_rs = zeros(nElec,1);
peak_power = zeros(nElec, nGoodTrial);
for i = 1:nElec
    curr_peak = round(peaks(i));
    peak_idx = find(curr_peak == pow.freq');

    
    [r,pcorr] = corr(log10(pow.powspctrm(:,i,peak_idx)),good_events(:,2) - good_events(:,1));
    corr_rs(i) = r;
    figure(2); clf
    scatter(log10(pow.powspctrm(:,i,peak_idx)),good_events(:,2) - good_events(:,1));
    ylabel('log(power)')
    xlabel('RT')
    pause(0.01);
    
    %add to peak data
    peak_power(i,:) = log10(pow.powspctrm(:,i,peak_idx));
end


%things for r
region = regions(:,1);
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
    title(regions{i})
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


