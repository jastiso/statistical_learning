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
img_dir = ['/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_img/subj', subj];

% make diractories
if ~exist(img_dir, 'dir')
    mkdir(img_dir);
end

% load stuff
load([save_dir, subj, '/raw_data.mat'])
load([save_dir, subj, '/header.mat'], 'labels', 'srate', 'HUP_ID', 'subj')
load([save_dir, subj, '/events.mat'])

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

avg_trans = mean(psds(:,:,logical(is_crosscluster)),3);
avg_within = mean(psds(:,:,~logical(is_crosscluster)),3);
se_trans = std(psds(:,:,logical(is_crosscluster)),[],3)./sqrt(sum(is_crosscluster));
se_within = std(psds(:,:,~logical(is_crosscluster)),[],3)./sqrt(sum(~is_crosscluster));
% remove PD
elec_labels = labels([1:2, 5:end]);

for i = 1:nElec
    figure(1); clf;
    plot(log10(freqs), avg_trans(:,i), "color", rgb("steelblue"), "linewidth", 2); hold on
    shade_plot(log10(freqs)', avg_trans(:,i),  se_trans(:,i), rgb("slategrey"));
    plot(log10(freqs), avg_within(:,i),  "color", rgb("chocolate"), "linewidth", 2);
    shade_plot(log10(freqs)', avg_within(:,i), se_within(:,i), rgb("sandybrown"));
    legend([{'Cross'}, {'Cross SE'}, {'Cross SE'},{'Within'}, {'Within SE'}])
    xlim([0,log10(256)])
    title(elec_labels{i})
    saveas(gca, [img_dir, '/psd_', elec_labels{i}, '.png'], 'png')
end

%% Plot raw data by type and elec

rt = events(:,2) - events(:,1);
%ERPs = nan(nElec, max(rt), nTrial);
%x_vect = (1:(max(rt)+1))./srate;
x_vect = (1:250)./srate;
ERPs = nan(nElec, 250, nTrial);

for t = 1:nTrial
   %curr = data(:,events(t,1):(events(t,2)));
   curr = data(:,(events(t,2)-249):events(t,2));
   %ERPs(:,1:rt(t)+1,t) = curr;
   ERPs(:,:,t) = curr;
end

avg_erp_trans = squeeze(nanmean(ERPs(:,:,logical(is_crosscluster)),3));
avg_erp_within = squeeze(nanmean(ERPs(:,:,logical(~is_crosscluster)),3));
se_erp_trans = squeeze(nanstd(ERPs(:,:,logical(is_crosscluster)),[],3))./sqrt(sum(is_crosscluster));
se_erp_within = squeeze(nanstd(ERPs(:,:,logical(is_crosscluster)),[],3))./sqrt(sum(~is_crosscluster));

for i = 1:nElec
    figure(1); clf;
    plot(x_vect, avg_erp_trans(i,:),  "color", rgb("steelblue"), 'linewidth', 2); hold on
    shade_plot(x_vect, avg_erp_trans(i,:),  se_erp_trans(i,:), rgb("slategrey"), 0.4);
    plot(x_vect, avg_erp_within(i,:),  "color", rgb("chocolate"), 'linewidth', 2); hold on
    shade_plot(x_vect, avg_erp_within(i,:),  se_erp_within(i,:), rgb("sandybrown"), 0.4);
    legend([{'Cross'}, {'Cross SE'}, {'Within'}, {'Within SE'}])
    title(elec_labels{i})
    saveas(gca, [img_dir, '/erp_', elec_labels{i}, '_resp.png'], 'png')
end
