%% Calculate theta power for every trial

addpath(genpath('/Users/stiso/Documents/MATLAB/IRASA/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
% define variables
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

theta = [4, 12];

%% get extra event variables

trans_nodes = [4,5,9,0];
module0 = [0 1 2 3 4];
module1 = [5 6 7 8 9];

%get only good trials
good_trials = logical(correct) & cutoff;
good_events = events(logical(correct) & cutoff,:);

trans_idx = false(nTrial,1);
for i = 2:nTrial
    if any(trans_nodes == walk(i-1)) && any(trans_nodes == walk(i))
        trans_idx(i) = true;
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
cc_events = good_events(trans_idx,:);


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


%% Get power via hilbert transform...
% filter from Wilson paper: ?band-pass filtered in the ? band (4?10 Hz) 
% using digital filters constructed via the Parks-McClellan 
% (Chebychev, similar to windowed sinc, or firws in fieldtrip) optimal equiripple FIR filter design. 
% Transition bands were 4 Hz?4.5 Hz and 10 Hz?10.5 Hz. 
% Maximal ripple was 0.05 in the stop bands and 0.01 in the pass band.

% Since this method is exact, you need to use a narrow band filter to
% remove as much noise as possible

% However, using a narrow band filter can distort your signal... it is
% important to plot your filter waveform

% This method does not take into account the different peak
% frequencies...but also maybe makes more sense for electrodes that have
% multiple peaks

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = theta;
cfg.bpfiltdf = 0.5; % bandpass transition width
cfg.bsfiltdf = 0.5; % bandstop transition width
cfg.bpfiltdev = 0.01; % bandpass max passband deviation
cfg.bsfiltdev = 0.05; % bandstp max passband deviation
cfg.bpfilttype = 'firws'; % or 'firls' (slower), but avoid the default 'but' (= not well-suited for hilbert phase estimate)
cfg.hilbert = 'complex';

inst_amp = ft_preprocessing(cfg, ft_data);


%% plot response (only have to do this once)
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

%% Compare


corr_ps = zeros(nElec,1);
peak_power = zeros(nElec, nGoodTrial);
for i = 1:nElec
    curr_peak = round(peaks(i));
    peak_idx = find(curr_peak == pow.freq');

    
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


%things for r
region = AAL(:,1);
good_trial_idx = trial(good_trials);

% get data for python anaylsis
save([r_dir, 'subj', subj, '/peak_power.mat'], 'peak_power', 'elec_labels', 'region', 'good_trial_idx', 'module_idx')

%% Get instantaneous amp from hilbert transformed data

ia_mat = zeros(nElec, nGoodTrial);
for i = 1:nGoodTrial
    %add to data structure
    % is the mean what I want here?
    ia_mat(:,i) = mean(abs(inst_amp.trial{i}),2);
    
    figure(1); clf
    plot(abs(inst_amp.trial{i})', 'linewidth', 2);
    ylabel('amp')
    xlabel('time')
    %pause(0.01);  
end

for i = 1:nElec
    figure(2); clf
    scatter(ia_mat(i,:),good_events(:,2) - good_events(:,1));
    ylabel('log(power)')
    xlabel('RT')
    %pause(0.01);
end

% get instantaneous amp for python anaylsis
save([r_dir, 'subj', subj, '/inst_amp.mat'], 'ia_mat', 'elec_labels', 'region', 'good_trial_idx', 'module_idx')

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


