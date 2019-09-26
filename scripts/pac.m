%% PAC
% get phase-amplitude coupling for oscillation in theta and broadband gamma
% get MI via KL distance method
% compare to surrogate
% control for power in theta

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
nElec = numel(elec_labels);

% freq variables
hg = [70:250];


%% Channel Selection

% only elecs with signigifant peak
elec_labels = elec_labels(logical(ps),:);
peaks = peaks((logical(ps)));
AAL = AAL(logical(ps),:);
nElec = size(elec_labels,1);

%% Get phase at peak index

% get one for each elec
phase = cell(nElec,1);
for i = 1:numel(elec_labels)
    % get peak oscillation for this elec
    curr_peak = peaks(i);
    
    
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.channel = elec_labels{i};
    cfg.bpfreq = [curr_peak-3 curr_peak+3];
    cfg.bpfiltdf = 0.5; % bandpass transition width
    cfg.bsfiltdf = 0.5; % bandstop transition width
    cfg.bpfiltdev = 0.01; % bandpass max passband deviation
    cfg.bsfiltdev = 0.05; % bandstp max passband deviation
    cfg.bpfilttype = 'firws'; % or 'firls' (slower), but avoid the default 'but' (= not well-suited for hilbert phase estimate)
    cfg.hilbert = 'complex';
    
    curr_phase = ft_preprocessing(cfg, ft_data);
    
    % get inst phase
    for t = 1:nGoodTrial
        curr_phase.trial{t} = (atan2(imag(curr_phase.trial{t}),real(curr_phase.trial{t})));
    end
    
    phase{i} = curr_phase;
end


%% get broadband power

% get one for each elec
amp = cell(nElec,1);
for i = 1:numel(elec_labels)

    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.channel = elec_labels{i};
    cfg.bpfreq = [70, 250];
    cfg.bpfiltdf = 0.5; % bandpass transition width
    cfg.bsfiltdf = 0.5; % bandstop transition width
    cfg.bpfiltdev = 0.01; % bandpass max passband deviation
    cfg.bsfiltdev = 0.05; % bandstp max passband deviation
    cfg.bpfilttype = 'firws'; % or 'firls' (slower), but avoid the default 'but' (= not well-suited for hilbert phase estimate)
    cfg.hilbert = 'complex';
    
    curr_amp = ft_preprocessing(cfg, ft_data);
    
    % get amp envelope
    for t = 1:nGoodTrial
        curr_amp.trial{t} = abs(curr_amp.trial{t});
    end
    
    amp{i} = curr_amp;
end

%% bin phases

nBin = 18;
binned_amp = zeros(nBin,nElec);
% save indicies for surrogate analysis later
indices = cell(nElec, 1);

% -pi/2 - pi/2
angles = linspace(-pi, pi, 100);
[~,edges] = discretize(angles, nBin);

for i = 1:nElec
    % vectorize
    phase_vect = cell2mat(phase{i}.trial);
    amp_vect = cell2mat(amp{i}.trial);
    
    % get indices for bins
    ind = discretize(phase_vect,edges);
    indices{i} = ind;
    
    % avg amp over bins
    binned_amp(:,i) = accumarray(ind', amp_vect, [nBin, 1], @mean);
    
end

% normalize over all bins
binned_amp = binned_amp./sum(binned_amp,1);

% plot histogram
for i = 1:nElec
    bar(edges(1:end-1), binned_amp(:,i), 'FaceColor',rgb('steelblue'),'EdgeColor','white','facealpha', 0.8)
    title(['PAC ', AAL{i,1}])
    saveas(gca, [img_dir, '/PAC_', elec_labels{i}, '.png'], 'png')
end


%% Get modulation index

uni_dist = repmat(unifpdf(linspace(-pi, pi, nBin),-pi,pi), nElec, 1);
KL_dist = KLDiv(binned_amp', uni_dist);
mod_idx = KL_dist./log10(nBin);


%% Make surrogate data

% load sig elecs from linear models
ramp = readtable([r_dir, 'subj' subj, '/ramp_stats_hg_osc.csv']);
mod = readtable([r_dir, 'subj' subj, '/mod_stats_hg_osc.csv']);
sig_contrast = mod.p < 0.05 | ramp.p < 0.05;
sig_elecs = ramp.elecs(sig_contrast);

nSim = 200;
uni_dist = repmat(unifpdf(linspace(-pi, pi, nBin),-pi,pi), nSim, 1);
surrogate_dist = zeros(sum(sig_contrast), nBin, nSim);
mi_surr = zeros(nSim, sum(sig_contrast));

for i = 1:sum(sig_contrast)
    elec_idx = find(strcmpi(elec_labels, sig_elecs{i}));
    
    % amplitudes
    amp_vect = cell2mat(amp{i}.trial);
    for n = 1:nSim
        % shuffle indices
        shuff = shuffle(indices{elec_idx});
        
        % get avg amp
        surrogate_dist(i, :, n) = accumarray(shuff', amp_vect, [nBin, 1], @mean);
        surrogate_dist(i, :, n) = surrogate_dist(i, :, n)./sum(surrogate_dist(i, :, n));
    end
    
    KL_surr = KLDiv(squeeze(surrogate_dist(i, :, :))', uni_dist);
    mi_surr(:,i) = KL_surr./log10(nBin);
    
    % plot
    figure(i); clf
    histogram(mi_surr(:,i), 'Normalization', 'probability', 'FaceColor',rgb('steelblue'),'EdgeColor','white','facealpha', 0.8); hold on
    plot([mod_idx(elec_idx), mod_idx(elec_idx)], [0, 1], 'red', 'linewidth', 3)
    title(['MI vs Surrogate ', AAL{elec_idx,1}])
    saveas(gca, [img_dir, '/PAC_surr_', elec_labels{elec_idx}, '.png'], 'png')
    
end

for i = 1:sum(sig_contrast)
    elec_idx = find(strcmpi(elec_labels, sig_elecs{i}));
    % z-score
    z = (mod_idx(elec_idx) - mean(mi_surr(:,i)))/std(mi_surr(:,i))
end


%% control for power in theta




