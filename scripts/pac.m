%% PAC
% get phase-amplitude coupling for oscillation in theta and broadband gamma
% get MI via KL distance method
% compare to surrogate
% control for power in theta

clear

addpath(genpath('/Users/stiso/Documents/MATLAB/IRASA/'))
addpath(genpath('/Users/stiso/Documents/Code/graph_learning/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
% define variables
subj = '4';
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
nElec = numel(elec_labels);

% freq variables
hg = [70:250];
lf_width = 2;

%% Channel Selection

% only elecs with signigifant peak
peak_elec_labels = elec_labels(logical(ps),:);
peaks = peaks((logical(ps)));
%regions = regions(logical(ps),:);
nPeakElec = size(peak_elec_labels,1);

%% Get phase at peak index

% get one for each elec
phase = cell(nPeakElec,1);
for i = 1:numel(peak_elec_labels)
    % get peak oscillation for this elec
    curr_peak = peaks(i);
    
    
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.channel = peak_elec_labels{i};
    cfg.bpfreq = [curr_peak-lf_width curr_peak+lf_width];
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

cfg = [];
cfg.bpfilter = 'yes';
cfg.channel = elec_labels;
cfg.bpfreq = [70, 150]; % needs to be have a width around center frquency greater than twice the modulating low frequency (Aru et al 2019)
cfg.bpfiltdf = 0.5; % bandpass transition width
cfg.bsfiltdf = 0.5; % bandstop transition width
cfg.bpfiltdev = 0.01; % bandpass max passband deviation
cfg.bsfiltdev = 0.05; % bandstp max passband deviation
cfg.bpfilttype = 'firws'; % or 'firls' (slower), but avoid the default 'but' (= not well-suited for hilbert phase estimate)
cfg.hilbert = 'complex';

pow = ft_preprocessing(cfg, ft_data);

% get amp envelope
for i = 1:nElec
    for t = 1:nGoodTrial
        amp{i}.trial{t} = abs(pow.trial{t}(i,:));
    end
end



%% bin phases

nBin = 18;
binned_amp = zeros(nBin,nPeakElec);
% save indicies for surrogate analysis later
indices = cell(nPeakElec, 1);

% -pi to pi
angles = linspace(-pi, pi, 100);
[~,edges] = discretize(angles, nBin);

cnt = 0;
for i = 1:nElec
    if ps(i) % if theres an oscillation
        cnt = cnt + 1;
        % vectorize
        phase_vect = cell2mat(phase{cnt}.trial);
        amp_vect = cell2mat(amp{i}.trial);
        
        % get indices for bins
        ind = discretize(phase_vect,edges);
        indices{cnt} = ind;
        
        % avg amp over bins
        binned_amp(:,cnt) = accumarray(ind', amp_vect, [nBin, 1], @mean);
    end
end

% normalize over all bins
binned_amp = binned_amp./sum(binned_amp,1);

% plot histogram
for i = 1:nPeakElec
    bar(edges(1:end-1), binned_amp(:,i), 'FaceColor',rgb('steelblue'),'EdgeColor','white','facealpha', 0.8)
    title(['PAC ', regions{i}])
    saveas(gca, [img_dir, '/PAC_', peak_elec_labels{i}, '.png'], 'png')
end

% Get modulation index

uni_dist = repmat(unifpdf(linspace(-pi, pi, nBin),-pi,pi), nPeakElec, 1);
mod_idx = mod_index(binned_amp', uni_dist);

%% Make surrogate data for within elec PAC
% only check if this elec has oscillations
% using surrogate approach that cuts the raw data at a single point and
% flips it (Aru et al)

% load sig elecs from linear models
osc = readtable([r_dir, 'subj' subj, '/max_ent_stats_hg_osc.csv']);
reg = readtable([r_dir, 'subj' subj, '/max_ent_stats.csv']);
sig_contrast = [osc.p < 0.05; reg.p < 0.05];
sig_elecs = [osc.elecs(osc.p < 0.05); reg.elecs(reg.p < 0.05)];

nSim = 500;
min_cut = 1;
uni_dist = repmat(unifpdf(linspace(-pi, pi, nBin),-pi,pi), nSim, 1);
surrogate_dist = zeros(sum(sig_contrast), nBin, nSim);
mi_surr = zeros(nSim, sum(sig_contrast));

for i = 1:sum(sig_contrast)
    elec_idx = find(strcmpi(peak_elec_labels, sig_elecs{i}));
    
    if ~isempty(elec_idx) % if this elec has oscillations
        % amplitudes
        amp_vect = cell2mat(amp{elec_idx}.trial);
        
        surrogate_dist(i,:,:) = pac_surr_cut(nSim, min_cut, phase{i}.trial, amp_vect, nBin, edges);
        
        mi_surr(:,i) = mod_index(squeeze(surrogate_dist(i, :, :))', uni_dist);
        
        % plot
        figure(1); clf
        histogram(mi_surr(:,i), 'Normalization', 'probability', 'FaceColor',rgb('steelblue'),'EdgeColor','white','facealpha', 0.8); hold on
        plot([mod_idx(elec_idx), mod_idx(elec_idx)], [0, .3], 'red', 'linewidth', 3)
        title(['MI vs Surrogate ', regions{elec_idx,1}])
        saveas(gca, [img_dir, '/PAC_surr_', peak_elec_labels{elec_idx}, '.png'], 'png')
    end
end

for i = 1:sum(sig_contrast)
    elec_idx = find(strcmpi(peak_elec_labels, sig_elecs{i}));
    % z-score
    fprintf('For electrode %s, z = %d\n', peak_elec_labels{elec_idx},...
        (mod_idx(elec_idx) - mean(mi_surr(:,i)))/std(mi_surr(:,i)));
    % show peak phase
    [~,ind] = max(binned_amp(:,elec_idx));
    fprintf('The peak of the bined amplitudes is at %d pi \n', edges(ind)./pi)
end


%% Test for coupling to the hippocampus from other regions

reg = readtable([r_dir, 'subj' subj, '/max_ent_stats_hg.csv']);
sig_contrast = reg.p < 0.05;
sig_elecs = reg.elecs(reg.p < 0.05);
sig_regions = reg.region(sig_contrast);

hpc_idx = cellfun(@(x) contains(x, 'hippocamp', 'IgnoreCase', true), regions);
hpc_phase = indices(hpc_idx);

binned_amp_hpc = zeros(nBin,sum(sig_contrast),sum(hpc_idx));

for i = 1:sum(sig_contrast)
    curr_elec = find(strcmpi(elec_labels, sig_elecs{i}));
    
    for j = 1:sum(hpc_idx)
        
        % vectorize
        amp_vect = cell2mat(amp{curr_elec}.trial);
        
        % get indices for bins
        ind = hpc_phase{j};
        
        % avg amp over bins
        binned_amp_hpc(:,i,j) = accumarray(ind', amp_vect, [nBin, 1], @mean);
        
        figure(1); clf;
        bar(edges(1:end-1), binned_amp_hpc(:,i,j), 'FaceColor',rgb('steelblue'),'EdgeColor','white','facealpha', 0.8)
        title(['PAC ', sig_regions{i}])
        saveas(gca, [img_dir, '/PAC_HPC_', elec_labels{curr_elec}, '_', num2str(j), '.png'], 'png')
    end
end

% normalize over all bins
binned_amp_hpc = binned_amp_hpc./sum(binned_amp_hpc,1);

% get MI
uni_dist = repmat(unifpdf(linspace(-pi, pi, nBin),-pi,pi), sum(sig_contrast), 1);
mi_hpc = zeros(sum(sig_contrast),sum(hpc_idx));
for i = 1:sum(hpc_idx)
    mi_hpc(:,i) = mod_index(binned_amp_hpc(:,:,i)', uni_dist);
end

% test surrogate
uni_dist = repmat(unifpdf(linspace(-pi, pi, nBin),-pi,pi), nSim, 1);
surrogate_dist_hpc = zeros(sum(sig_contrast), sum(hpc_idx), nBin, nSim);
mi_surr_hpc = zeros(nSim, sum(sig_contrast), sum(hpc_idx));
hpc = find(hpc_idx);

for i = 1:sum(sig_contrast)
    elec_idx = find(strcmpi(elec_labels, sig_elecs{i}));
    
    for j = 1:sum(hpc_idx)
        curr_hpc = hpc(j);
        
        % amplitudes
        amp_vect = cell2mat(amp{elec_idx}.trial);
        
        surrogate_dist_hpc(i, j, :, :) = pac_surr_cut(nSim, min_cut, phase{curr_hpc}.trial, amp_vect, nBin, edges);
        
        mi_surr_hpc(:,i,j) = mod_index(squeeze(surrogate_dist_hpc(i, j, :, :))', uni_dist);
        
        % plot
        figure(1); clf
        histogram(mi_surr_hpc(:,i,j), 'Normalization', 'probability', 'FaceColor',rgb('steelblue'),'EdgeColor','white','facealpha', 0.8); hold on
        plot([mi_hpc(i, j), mi_hpc(i, j)], [0, .3], 'red', 'linewidth', 3)
        title(['MI vs Surrogate ', sig_regions{i}])
        saveas(gca, [img_dir, '/PAC_surr_HPC_', elec_labels{elec_idx}, '_', num2str(j), '.png'], 'png')
    end
end

for i = 1:sum(sig_contrast)
    elec_idx = find(strcmpi(elec_labels, sig_elecs{i}));
    
    for j = 1:sum(hpc_idx)
        % z-score
        fprintf('\nFor electrode %s with %s, z = %d\n', elec_labels{elec_idx},...
            peak_elec_labels{hpc(j)}, round((mi_hpc(i, j) - mean(mi_surr_hpc(:,i,j)))/std(mi_surr_hpc(:,i,j)),2));
        % show peak phase
        [~,ind] = max(binned_amp_hpc(:,i,j));
        fprintf('The peak of the bined amplitudes is at %d pi \n', edges(ind)./pi)
    end
end

