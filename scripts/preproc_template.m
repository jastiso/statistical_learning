%% Preprocessing for subj 6, HUP187
clear

addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current/'))
addpath(genpath('/Users/stiso/Documents/Code/graph_learning/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')

% define variables
RID = '529';
subj = '6';
sessions = [{'1'}, {'2'}];
save_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_raw/';

% get regions regions
regions = readtable(['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_raw/', subj, '/RID', RID, '/electrodenames_coordinates_native_and_T1.csv']);
regions = regions(:,1:2); % get only names and labels

load([save_dir, subj, '/header_sess1.mat'], 'elec_labels', 'srate', 'HUP_ID', 'subj') % these should be the same
save([save_dir, subj, '/header_sess1.mat'], 'elec_labels', 'srate', 'HUP_ID', 'subj', 'regions') % these should be the same
data_all = struct('sess', []);
for i = 1:numel(sessions)
    data_all(i).sess = load([save_dir, subj, '/raw_data_sess', sessions{i}, '.mat'], 'data');
    data_all(i).sess = data_all(i).sess.data;
end



%% Remove noisy elecs

% sidebar: align the region names
[regions] = region_match(elec_labels, regions);

% remove elecs that are out of the brain
out_of_brain = cellfun(@(x) isempty(x), regions);
data = data(~out_of_brain, :);
regions = regions(~out_of_brain);
elec_labels = elec_labels(~out_of_brain);

%data-driven: find channels with large kurtosis
rmv = false(numel(elec_labels,1));
for i = 1:numel(sessions)
    rmv = rmv | reject_elecs(data_all(i).sess, 2, srate); % thr of 2 stds
end

eegplot(data_all(1).sess(rmv,:), 'srate', srate)

eegplot(data_all(2).sess(~rmv,:), 'srate', srate)

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

out_all = struct('sess', []);
marker_all = struct('sess', []);

for s = 1:numel(sessions)
    [out,MARKER] = spike_detector_hilbert_v16_byISARG(data_all(s).sess', srate);
    
    % select for only spikes in many channels
    win = 0.05;
    min_chan = 4;
    nSamp = size(MARKER.d,1);
    nSpike = numel(out.pos);
    kept_spike = false(size(out.pos));
    all_soz = []; % when I get SOZ channels, change this
    discharge_tol=0.005; % taken from spike function
    
    
    for i = 1:nSpike
        curr_pos = out.pos(i);
        curr_chan = out.chan(i);
        
        if kept_spike(i) == 0
            win_spike = (out.pos > curr_pos & out.pos < curr_pos + win);
            win_chan = out.chan(win_spike);
            
            if ~isempty(all_soz) % eventually we want only spikes that generalize to other channels
                if sum(intersect(win_chan,all_soz)) >= min_chan
                    kept_spike(win_spike) = true;
                end
            else
                if sum(unique(win_chan)) >= min_chan+1 % if not SOZ marked have slightly stricter cutoff
                    kept_spike(win_spike) = true;
                end
            end
        end
    end
    % slect only good spikes
    out_clean.pos = out.pos(kept_spike);
    out_clean.dur = out.dur(kept_spike);
    out_clean.chan = out.chan(kept_spike);
    out_clean.weight = out.weight(kept_spike);
    out_clean.con = out.con(kept_spike);
    fprintf('Found %d spikes\n"', sum(kept_spike));
    
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

% get concatenated events
load([save_dir, subj, '/events.mat']); % in samples

load([save_dir, subj, '/task_data.mat']) % should already have both sessions

nTrial = size(walk,2);
trans_nodes = [4,5,9,0];
module0 = [0 1 2 3 4];
module1 = [5 6 7 8 9];

%get only good trials
good_trials = logical(correct) & cutoff;
good_events = events(good_trials,:);

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
% save walk with ECoG data, for most people this will just be the full walk
rec_walk = walk(1:514);

save([save_dir, subj, '/good_events.mat'], 'good_events', 'trans_idx', 'module_idx', 'good_trials', 'rec_walk');

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