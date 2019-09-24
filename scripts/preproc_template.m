%% Preprocessing for subj X, HUPXXX

addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current/'))

% define variables
HUP_ID = 'HUP187';
RID = '522';
subj = '2';
save_dir = '/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_raw/';

% get AAL regions
load('/Users/stiso/Documents/MATLAB/AAL_regions_by_subject.mat')
AAL_all = data;
eval(['AAL = AAL_all.RID', RID, ';']);

load([save_dir, subj, '/raw_data.mat'], 'data')
load([save_dir, subj, '/header.mat'], 'elec_labels', 'srate', 'HUP_ID', 'subj')
AAL = AAL(~rmv, :);

%% Remove noisy elecs

% sidebar: align the AAL names
[unaligned_elecs] = AAL_match(elec_labels)

%data-driven: find channels with large kurtosis
rmv = reject_elecs(data, 2); % thr of 2 stds

eegplot(data(rmv,:), 'srate', srate)

eegplot(data(~rmv,:), 'srate', srate)

data = data(~rmv,:);
elec_labels = elec_labels(~rmv,:);

%% Look at line noise
% check for things outside of 60Hz...especially in theta/alpha range

spectopo(data, 0, srate)

% filter out 60 Hz harmonics
[b,a] = butter(4, [59/(srate/2), 61/(srate/2)], 'stop');
data = filtfilt(b,a,data')';

[b,a] = butter(4, [119/(srate/2), 121/(srate/2)], 'stop');
data = filtfilt(b,a,data')';

[b,a] = butter(4, [179/(srate/2), 181/(srate/2)], 'stop');
data = filtfilt(b,a,data')';

%% CAR
% if data looks to have different levels of noise based on elec, might want
% to do this in groups

data = data - mean(data,2); % demean
data = detrend(data')'; % as opposed to low pass filtering
data = get_CAR(data); % CAR by group


%% Get ictal spikes


save([save_dir, subj, '/data_clean.mat'], 'data')
save([save_dir, subj, '/header_clean.mat'], 'elec_labels', 'srate', 'HUP_ID', 'subj', 'AAL')

%% Get extra event fields

load([save_dir, subj, '/events.mat']) % in samples
load([save_dir, subj, '/task_data.mat'])

nTrial = size(events,1);
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

save([save_dir, subj, '/good_events.mat'], 'good_events', 'trans_idx', 'module_idx');

%% Get into fieldtrip format

ft_data = fieldtrip_format(data, srate, elec_labels, [good_events, zeros(size(good_events,1),1)]);

% demean and detrend data
for i = 1:numel(ft_data.trial)
    ft_data.trial{i} = ft_data.trial{i} - mean(ft_data.trial{i},2);
    ft_data.trial{i} = detrend(ft_data.trial{i});
end
save([save_dir, subj, '/ft_data.mat'], 'ft_data')