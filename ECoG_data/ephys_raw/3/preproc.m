%% Preprocessing for subj 1
clc
clear

%%
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current/'))
addpath(genpath('/Users/stiso/Documents/Code/graph_learning/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')

% define variables
RID = '502';
subj = '3';
sessions = [{'1'}];
save_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_raw/';

% get regions regions
regions = readtable(['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_raw/', subj, '/RID', RID, '/electrodenames_coordinates_native_and_T1.csv']);
regions = regions(:,1:2); % get only names and labels

load([save_dir, subj, '/header.mat'], 'elec_labels', 'srate', 'HUP_ID', 'subj') % these should be the same
save([save_dir, subj, '/header.mat'], 'elec_labels', 'srate', 'HUP_ID', 'subj', 'regions') % these should be the same
if numel(sessions) > 1
    data_all = struct('sess', []);
    for i = 1:numel(sessions)
        data_all(i).sess = load([save_dir, subj, '/raw_data_sess', sessions{i}, '.mat'], 'data');
        data_all(i).sess = data_all(i).sess.data;
    end
else
   data_all = struct('sess', []);
    data_all(1).sess = load([save_dir, subj, '/raw_data.mat'], 'data');
    data_all(1).sess = data_all(1).sess.data; 
end


%% Remove noisy elecs

% sidebar: align the region names
[regions] = region_match(elec_labels, regions);

% remove elecs that are out of the brain
out_of_brain = cellfun(@(x) isempty(x), regions);
for i = 1:numel(sessions)
    data_all(i).sess = data_all(i).sess(~out_of_brain,:);
end
regions = regions(~out_of_brain);
elec_labels = elec_labels(~out_of_brain);

%data-driven: find channels with large kurtosis
rmv = false(numel(elec_labels,1));
for i = 1:numel(sessions)
    rmv = rmv | reject_elecs(data_all(i).sess, 2, srate); % thr of 2 stds
end
rmv(113) = 1;
eegplot(data_all(1).sess(rmv,:), 'srate', srate)

for i = 1:numel(sessions)
    data_all(i).sess = data_all(i).sess(~rmv,:);
end
elec_labels = elec_labels(~rmv,:);
regions = regions(~rmv, :);

%% Look at line noise
% check for things outside of 60Hz...especially in theta/alpha range

fft_plot(data_all(1).sess', 500, srate)

% filter out 60 Hz harmonics
for i = 1:numel(sessions)
    [b,a] = butter(4, [59/(srate/2), 61/(srate/2)], 'stop');
    data_all(i).sess = filtfilt(b,a,data_all(i).sess')';
    
    [b,a] = butter(4, [119/(srate/2), 121/(srate/2)], 'stop');
    data_all(i).sess = filtfilt(b,a,data_all(i).sess')';
    
    [b,a] = butter(4, [179/(srate/2), 181/(srate/2)], 'stop');
    data_all(i).sess = filtfilt(b,a,data_all(i).sess')';
end

%% CAR
% if data looks to have different levels of noise based on elec, might want
% to do this in groups

for i = 1:numel(sessions)
    data_all(i).sess = data_all(i).sess - mean(data_all(i).sess,2); % demean
    data_all(i).sess = detrend(data_all(i).sess')'; % as opposed to low pass filtering
    data_all(i).sess = get_CAR(data_all(i).sess, elec_labels); % CAR by group
end

%% Get ictal spikes
% doing this using algorithm from Kate's paper
min_chan = 4; % minimum number of channels that need to be recruited
win = 0.05; % size of the window to look for the minimum number of channels, in seconds
seq = 0.015; % 15ms for spikes within a sequence, taken from Erin Conrads Brain paper
thr = 0.002; % reject spike if they spread to a lot of channels in less than 2ms (larger than paper), also from Erins paper
out_all = struct('sess', []);
marker_all = struct('sess', []);
discharge_tol=0.005;

for s = 1:numel(sessions)
    [out,MARKER] = ...
spike_detector_hilbert_v16_byISARG(data_all(s).sess', srate);
% sort spikes by onset
[sort_pos,I] = sort(out.pos, 'ascend');
out.pos = sort_pos;
out.chan = out.chan(I);
out.con = out.con(I);
out.dur = out.dur(I);
out.weight = out.weight(I);
out.pdf = out.pdf(I);
out.seq = nan(size(out.pos));
seq_cnt = 0; % counter for sequences, you can reset it here because the same number will never be in the same 1s window later on

% eliminate some spikes
include_length = win;
nSamp = size(MARKER.d,1);
nSpike = numel(out.pos);
kept_spike = false(size(out.pos));

for i = 1:nSpike
curr_pos = out.pos(i);

if kept_spike(i) == 0
    win_spike = (out.pos > curr_pos & out.pos < (curr_pos + win));
    % add spikes within 15ms of the
    % last one
    last_spike = max(out.pos(win_spike));
    curr_sum = 0;
    while curr_sum < sum(win_spike) % stop when you stop adding new spikes
        curr_sum = sum(win_spike);
        win_spike = win_spike | (out.pos > last_spike & out.pos < (last_spike + seq));
        last_spike = max(out.pos(win_spike));
    end
    win_chan = out.chan(win_spike);

    % set sequence ID for all spikes
    % in this window
    seqs = struct('idx',[],'chan',[]);
    if sum(win_spike) > 0
        [~,leader_idx] = min(out.pos(win_spike));
        leader = win_chan(leader_idx);
        if sum(win_chan == leader) == 1
            out.seq(win_spike) = seq_cnt;
            seqs(1).idx = find(win_spike);
            seqs(1).chan = win_chan;
            seq_cnt = seq_cnt + 1;
        else
            % parse sequence at the
            % leader.
            end_pts = find(win_chan == leader);
            end_pts = [end_pts; numel(win_chan) + 1];
            win_idx = find(win_spike);
            for m = 1:(numel(end_pts)-1)
                seqs(m).chan = win_chan(end_pts(m):(end_pts(m+1)-1));
                seqs(m).idx = win_idx(end_pts(m):(end_pts(m+1)-1));
                out.seq(seqs(m).idx) = seq_cnt;
                seq_cnt = seq_cnt + 1;
            end
        end
    end
    % for each sequence, remove events that generalize
    % to 80% of elecs within 2ms,
    % keep if it spread to at least
    % 3 channels
    for m = 1:numel(seqs)
        time_between = diff(out.pos(seqs(m).idx));        
        if (numel(unique(seqs(m).chan)) >= min_chan) && ~(numel(unique(seqs(m).chan(time_between <= thr))) >= numel(elec_labels)*.5)
            kept_spike(seqs(m).idx) = true;
        end
    end
end
end
% select only good spikes
out_clean.pos = out.pos(kept_spike);
out_clean.dur = out.dur(kept_spike);
out_clean.chan = out.chan(kept_spike);
out_clean.weight = out.weight(kept_spike);
out_clean.con = out.con(kept_spike);
out_clean.seq = out.seq(kept_spike);
fprintf('This dataset had %d IEDs\n', numel(kept_spike))

% get new M
m_clean = zeros(size(MARKER.M));
for i=1:size(out_clean.pos,1)
m_clean(round(out_clean.pos(i)*MARKER.fs:out_clean.pos(i)*MARKER.fs+discharge_tol*MARKER.fs),...
    out_clean.chan(i))=out_clean.con(i);
end
    MARKER.m_clean = m_clean;
    
    out_all(s).sess = out_clean;
    marker_all(s).sess = MARKER;
end

% visualize - check that these match
eegplot(marker_all(1).sess.d', 'srate', MARKER.fs)
eegplot(marker_all(1).sess.m_clean', 'srate', MARKER.fs)

save([save_dir, subj, '/data_clean.mat'], 'data_all')
save([save_dir, subj, '/header_clean.mat'], 'elec_labels', 'srate', 'HUP_ID', 'subj', 'regions', 'sessions')

%% Get extra event fields

load([save_dir, subj, '/events.mat']) % in samples
load([save_dir, subj, '/task_data.mat'])

nTrial = size(events,1);
trans_nodes = [4,5,9,0];
module0 = [0 1 2 3 4];
module1 = [5 6 7 8 9];

% get trials with IEDs
ictal_trials = zeros(1,nTrial);
for i = 1:nTrial
    if any((out_clean.pos > (events(i,1)/srate)) & (out_clean.pos < (events(i,2)/srate)))
        ictal_trials(i) = 1;
    else
        ictal_trials(i) = 0;
    end
end

%get only good trials
good_trials = logical(correct) & cutoff & ~ictal_trials;
good_events = events(logical(correct) & cutoff & ~ictal_trials,:);

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
trans_idx = trans_idx(logical(correct) & cutoff & ~ictal_trials);
module_idx = module_idx(logical(correct) & cutoff & ~ictal_trials);

save([save_dir, subj, '/good_events.mat'], 'good_events', 'trans_idx', 'module_idx', 'good_trials', 'ictal_trials');

%% Get into fieldtrip format

for i = 1:numel(sessions)
    % get currrent sessions events
    sess_idx = good_events(:,3) == i;
    trl = [good_events(sess_idx,1:2), zeros(sum(sess_idx),1), good_events(sess_idx,3)];
    % epoch
    curr = fieldtrip_format(data_all(i).sess, srate, elec_labels, trl);
    % concetenate
    if i == 1
        ft_data = curr;
    else
        ft_data.trial = [ft_data.trial, curr.trial];
        ft_data.time = [ft_data.time, curr.time];
        ft_data.trialinfo = [ft_data.trialinfo; curr.trialinfo];
        ft_data.sampleinfo = [ft_data.sampleinfo; curr.sampleinfo];
    end
end

save([save_dir, subj, '/ft_data.mat'], 'ft_data')