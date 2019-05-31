%% Get trial indices, and ictal periods

addpath(genpath('/Users/stiso/Documents/MATLAB/ieeg-matlab-1.13.2/'))
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current/'))
% define variables
HUP_ID = 'HUP187';
subj = '2';
save_dir = '/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_raw/';

% load stuff
load([save_dir, subj, '/raw_pd.mat'], 'pd')
load([save_dir, subj, '/header.mat'], 'labels', 'srate', 'HUP_ID', 'subj')

%% Interactive plot for getting photodiode events

% events = find_events(data, srate, minITI, minDuration, event_thresh, show)

% for this task, ITI is always 50, event duration can range as long as RT
% can. All times are in ms, not samples
%
% you might have to play around with the threshold to get something that
% works

pd(isnan(pd)) = 0;

events = find_events(pd(1,:), srate, 40, 50, 50000, "show") ;
events = [events.onsets events.offsets];

%% Remove events that weren't marked well, or are demo
% here, only one bad trial at the end when the PD is unplugged

% if all trials are present, there should be 1000 at the end
events = events(11:end-1,:);

save([save_dir, subj, '/events.mat'], 'events')
