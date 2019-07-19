%% Exploratory Theta

% I want to do this in python eventually...but for the purposes of getting
% quick preliminary data for grants I have a matlab one too. Also for
% checking my python results

addpath(genpath('/Users/stiso/Documents/MATLAB/ieeg-matlab-1.13.2/'))
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
load([save_dir, subj, '/raw_data.mat'])
load([save_dir, subj, '/header.mat'], 'labels', 'srate', 'HUP_ID', 'subj')
load([save_dir, subj, '/events.mat'])
load([save_dir, subj, '/trans_idx.mat'])
elec_labels = labels([1:3, 6:end]);

% analysis varaibles
freqs = 0:2:srate/2;
nTrial = size(events,1);
nElec = size(data,1);
nFreq = numel(freqs);

%% Get PSD for each trial
% using Welch's method for noe

psds = zeros(nFreq, nElec, nTrial);

for t = 1:nTrial
    curr = data(:,events(t,1):events(t,2))';
    % pwelch works over columns
    [psd, freq] = pwelch(curr, [],[],[], srate, 'power');
    psds(:,:,t) = log10(psd);
end

%% Format for r

psds_vect = reshape(psds, [], 1);
freq_id = repmat(freq', 1, nElec*nTrial)';
elec_id = repmat(reshape(repmat(1:nElec, nFreq, 1), [] , 1), nTrial, 1);
trans_exp = reshape(repmat(is_crosscluster, 1, nElec*nFreq)', [], 1);
alt_exp = reshape(repmat(alt_idx, 1, nElec*nFreq)', [], 1);
trial = reshape(repmat(1:nTrial, 1, nElec*nFreq)', [], 1);
save([r_dir, '/subj', subj '/psds.mat'], 'psds_vect', 'elec_labels', 'elec_id', 'freq_id', 'trans_exp', 'alt_exp', 'trial')

%% Exploratory plot

avg_psd = mean(psds,3);
plot(freqs,avg_psd);

for i = 1:nElec
    figure(1); clf;
    plot(freqs, squeeze(psds(:,i,:)));
    xlim([0,256])
    pause(0.1)
end



%% Average by cross cluster

alt_idx = circshift(is_crosscluster,1);

avg_trans = mean(psds(:,:,logical(is_crosscluster)),3);
avg_alt = mean(psds(:,:,logical(alt_idx)),3);
avg_within = mean(psds(:,:,~logical(is_crosscluster)),3);
se_trans = std(psds(:,:,logical(is_crosscluster)),[],3)./sqrt(sum(is_crosscluster));
se_alt = std(psds(:,:,logical(alt_idx)),[],3)./sqrt(sum(alt_idx));
se_within = std(psds(:,:,~logical(is_crosscluster)),[],3)./sqrt(sum(~is_crosscluster));

% remove PD
elec_labels = labels([1:2, 5:end]);

for i = 1:nElec
    figure(1); clf;
    plot(log10(freqs), avg_trans(:,i), "color", rgb("steelblue"), "linewidth", 2); hold on
    shade_plot(log10(freqs)', avg_trans(:,i),  se_trans(:,i), rgb("slategrey"));
    plot(log10(freqs), avg_within(:,i),  "color", rgb("chocolate"), "linewidth", 2);
    shade_plot(log10(freqs)', avg_within(:,i), se_within(:,i), rgb("sandybrown"));
    plot(log10(freqs), avg_alt(:,i),  "color", rgb("hotpink"), "linewidth", 2);
    shade_plot(log10(freqs)', avg_alt(:,i), se_alt(:,i), rgb("pink"));
    legend([{'Cross'}, {'Cross SE'}, {'Cross SE'},{'Within'}, {'Within SE'}, {'Within SE'}, {'Alt'}, {'Alt SE'}])
    xlim([0,log10(256)])
    title(elec_labels{i})
    saveas(gca, [img_dir, '/psd_', elec_labels{i}, '.png'], 'png')
end


%% Stats

theta = [4,10];
beta = [24, 30];
theta_idx = (freqs >= theta(1)) & (freqs <= theta(2));
beta_idx = (freqs >= beta(1)) & (freqs <= beta(2));
p = zeros(nElec,1);
p_beta = zeros(nElec,1);

for i = 1:nElec
    theta_alt = squeeze(mean(psds(theta_idx,i,logical(alt_idx)),1));
    theta_trans = squeeze(mean(psds(theta_idx,i,logical(is_crosscluster)),1));
    p(i) = permtest(theta_alt, theta_trans);
end

for i = 1:nElec
    beta_alt = squeeze(mean(psds(beta_idx,i,logical(alt_idx)),1));
    beta_trans = squeeze(mean(psds(beta_idx,i,logical(is_crosscluster)),1));
    p_beta(i) = permtest(beta_alt, beta_trans);
end

[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p(24:35),0.05,'dep','yes');
[hb, crit_pb, adj_ci_cvrgb, adj_pb]=fdr_bh(p_beta(24:35),0.05,'dep','yes');

e = 24:35;
save([r_dir, '/subj', subj '/stats.mat'], 'adj_p', 'adj_pb', 'e')


%% Plot raw data by type and elec

rt = events(:,2) - events(:,1);
%ERPs = nan(nElec, max(rt), nTrial);
%x_vect = (1:(max(rt)+1))./srate;
x_vect = (1:500)./srate;
ERPs = nan(nElec, 500, nTrial);

for t = 1:nTrial
    %curr = data(:,events(t,1):(events(t,2)));
    curr = data(:,(events(t,2)):events(t,2) + 499);
    %ERPs(:,1:rt(t)+1,t) = curr;
    ERPs(:,:,t) = curr;
end

avg_erp_trans = squeeze(nanmean(ERPs(:,:,logical(is_crosscluster)),3));
avg_erp_alt = squeeze(nanmean(ERPs(:,:,logical(alt_idx)),3));
avg_erp_within = squeeze(nanmean(ERPs(:,:,logical(~is_crosscluster)),3));
se_erp_trans = squeeze(nanstd(ERPs(:,:,logical(is_crosscluster)),[],3))./sqrt(sum(is_crosscluster));
se_erp_alt = squeeze(nanstd(ERPs(:,:,logical(alt_idx)),[],3))./sqrt(sum(alt_idx));
se_erp_within = squeeze(nanstd(ERPs(:,:,logical(is_crosscluster)),[],3))./sqrt(sum(~is_crosscluster));

for i = 1:nElec
    figure(1); clf;
    plot(x_vect, avg_erp_trans(i,:),  "color", rgb("steelblue"), 'linewidth', 2); hold on
    shade_plot(x_vect, avg_erp_trans(i,:),  se_erp_trans(i,:), rgb("slategrey"), 0.4);
    plot(x_vect, avg_erp_within(i,:),  "color", rgb("chocolate"), 'linewidth', 2); hold on
    shade_plot(x_vect, avg_erp_within(i,:),  se_erp_within(i,:), rgb("sandybrown"), 0.4);
    plot(x_vect, avg_erp_alt(i,:),  "color", rgb("hotpink"), 'linewidth', 2); hold on
    shade_plot(x_vect, avg_erp_alt(i,:),  se_erp_alt(i,:), rgb("pink"), 0.4);
    legend([{'Cross'}, {'Cross SE'}, {'Within'}, {'Within SE'}, {'Alt'}, {'Alt SE'}])
    title(elec_labels{i})
    saveas(gca, [img_dir, '/erp_', elec_labels{i}, '_stim.png'], 'png')
end




%% IRASA
% separates scale free from oscillatory component of the power spectra.

%win_length = floor((min(rt) + 50)*(srate/1000)) - 1;
win_length = floor((3000)*(srate/1000)) - 1;
step = floor((1000)*(srate/1000));
nWin = floor(size(data,2)/win_length);

for j = 1:nElec
sig = zeros(win_length + 1, nWin);
for i = 1:nWin
    onset = events(i,1);
    %sig(:,i) = data(j,onset:(onset + win_length));
    sig(:,i) = data(j,ceil((i-1)*step) + 1 : ceil((i-1)*step)+win_length+1);
end

% spec = amri_sig_fractal(sig,srate,...); optional arguments frange,
% detrend, filter
spec = amri_sig_fractal(sig,srate, 'detrend', 1, 'filter', 1);

% plot
% show averaged fractal and oscillatory power spectrum
figure(1); clf;
subplot(2,1,1);
loglog(spec.freq,mean(spec.mixd,2),'b',  'linewidth', 2); hold on
loglog(spec.freq,mean(spec.frac,2),'r', 'linewidth', 2);
subplot(2,1,2);
plot(spec.freq, mean(spec.osci,2), 'linewidth', 2);
saveas(gca, [img_dir, '/IRASA_', elec_labels{j}, '.png'], 'png')
end
    
