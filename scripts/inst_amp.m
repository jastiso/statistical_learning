%% Instantaneous amplitude
clear

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
cfg.channel = elec_labels;
cfg.bpfreq = theta;
cfg.bpfiltdf = 0.5; % bandpass transition width
cfg.bsfiltdf = 0.5; % bandstop transition width
cfg.bpfiltdev = 0.01; % bandpass max passband deviation
cfg.bsfiltdev = 0.05; % bandstp max passband deviation
cfg.bpfilttype = 'firws'; % or 'firls' (slower), but avoid the default 'but' (= not well-suited for hilbert phase estimate)
cfg.hilbert = 'complex';

inst = ft_preprocessing(cfg, ft_data);


%% plot response (only have to do this once)

plot_data = filter_resp(ft_data, 4000, srate);
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

%% Get instantaneous amp from hilbert transformed data

ia_mat = zeros(nElec, nGoodTrial);
for i = 1:nGoodTrial
    %add to data structure
    % is the mean what I want here?
    ia_mat(:,i) = mean(abs(inst.trial{i}),2);
    
    figure(1); clf
    plot(abs(inst.trial{i})', 'linewidth', 2);
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

% get instantaneous amp for r anaylsis
region = AAL(:,1);
good_trial_idx = trial(good_trials);
save([r_dir, 'subj', subj, '/inst_amp.mat'], 'ia_mat', 'elec_labels', 'region', 'good_trial_idx', 'module_idx')

%% Plots

nNode = numel(unique(walk));
ia_node = zeros(nNode, nNode, nElec);

for i = 1:nElec
    
    count_node = zeros(nNode);
    
    cnt = 0; % not zero because j counter starts at 2
    for j = 2:numel(walk)
        prev = walk(j-1) + 1;
        curr = walk(j) + 1;
        if good_trials(j)
            cnt = cnt + 1;
            ia_node(prev,curr,i) = ia_node(prev,curr,i) + ia_mat(i,cnt);
            count_node(prev,curr) = count_node(prev,curr) + 1;
            
           
        end
    end
    ia_node(:,:,i) = ia_node(:,:,i)./count_node;
    
    %plot
    figure(3); clf
    imagesc(ia_node(:,:,i)); c = colorbar;
    title(AAL{i})
    xlabel('current node')
    ylabel('previous node')
    c.Label.String = 'power';
    saveas(gca, [img_dir, '/ia_by_trans_', elec_labels{i}, '.png'], 'png')
    
end
