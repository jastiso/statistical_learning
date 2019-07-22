%% Preprocessing for subj 2, HUP187

addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current/'))

% define variables
HUP_ID = 'HUP187';
subj = '2';
save_dir = '/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_raw/';

load([save_dir, subj, '/raw_data.mat'], 'data')
load([save_dir, subj, '/header.mat'], 'elec_labels', 'srate', 'HUP_ID', 'subj')

%% Remove noisy elecs

%data-driven: find channels with large kurtosis
rmv = reject_elecs(data, 2); % thr of 2 stds

eegplot(data(rmv,:), 'srate', srate)

eegplot(data(~rmv,:), 'srate', srate)

data = data(~rmv,:);
elec_labels = elec_labels(~rmv,:);

%% Look at line noise
% check for things outside of 60Hz...especially in theta/alpha range

spectopo(data, 0, srate)


%% CAR
% if data looks to have different levels of noise based on elec, might want
% to do this in groups

data = data - mean(data,2); % demean
data = detrend(data')'; % as opposed to low pass filtering
data = get_CAR(data); % CAR by group

%% Get ictal spikes