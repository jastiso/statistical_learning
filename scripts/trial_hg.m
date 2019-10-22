%% Calculate theta power for every trial
clear

addpath(genpath('/Users/stiso/Documents/MATLAB/IRASA/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
% define variables
subj = '2';
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

% analysis varaibles
nTrial = size(events,1);
nGoodTrial = size(good_events,1);
nElec = numel(elec_labels);

% freq variables
hg = [70:250];

% task variables - recorded walk (only different if the photodiode messed
% up, in which case a separate variable is saved in preprocessing)
if nTrial == numel(walk)
    rec_walk = walk;
end

%% Get power at hg index

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


%% Log transform and average


corr_rs = zeros(nElec,1);
hg_power = zeros(nElec, nGoodTrial);
for i = 1:nElec
    curr = squeeze(mean(log10(pow.powspctrm(:,i,:)),3));
    
    [r,pcorr] = corr(curr,good_events(:,2) - good_events(:,1));
    corr_rs(i) = r;
    figure(2); clf
    scatter(curr,good_events(:,2) - good_events(:,1));
    ylabel('log(power)')
    xlabel('RT')
    pause(0.01);
    
    %add to peak data
    hg_power(i,:) = curr;
end


%things for r
region = AAL(:,1);
good_trial_idx = trial(good_trials);

% get data for python anaylsis
save([r_dir, 'subj', subj, '/hg_power.mat'], 'hg_power', 'elec_labels', 'region', 'good_trial_idx', 'module_idx')



%% Plots

nNode = numel(unique(walk));
power_node = zeros(nNode, nNode, nElec);

for i = 1:nElec
    
    count_node = zeros(nNode);
    
    cnt = 0; % not zero because j counter starts at 2
    for j = 2:numel(rec_walk)
        prev = rec_walk(j-1) + 1;
        curr = rec_walk(j) + 1;
        if good_trials(j)
            cnt = cnt + 1;
            power_node(prev,curr,i) = power_node(prev,curr,i) + hg_power(i,cnt);
            count_node(prev,curr) = count_node(prev,curr) + 1;
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
    saveas(gca, [img_dir, '/hg_by_trans_', elec_labels{i}, '.png'], 'png')
    
end


