%% Get trial indices, and ictal periods
clear 

addpath(genpath('/Users/stiso/Documents/MATLAB/ieeg-matlab-1.13.2/'))
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current/'))
% define variables
subj = '18';
sess = ''; % can be empty
sess_flag = ~isempty(sess);
save_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_raw/';

% load stuff
if ~sess_flag
    load([save_dir, subj, '/raw_pd_sess1.mat'], 'pd')
    load([save_dir, subj, '/header_sess1.mat'], 'srate', 'HUP_ID', 'subj')
else
    load([save_dir, subj, '/raw_pd_sess', sess, '.mat'], 'pd')
    load([save_dir, subj, '/header_sess', sess, '.mat'], 'srate', 'HUP_ID', 'subj')
end

%% Interactive plot for getting photodiode events

% events = find_events(data, srate, minITI, minDuration, event_thresh, show)

% for this task, ITI is always 50, event duration can range as long as RT
% can. All times are in ms, not samples
%
% you might have to play around with the threshold to get something that
% works

try 
    if sess_flag
        load([save_dir, subj, '/event_params_sess', sess, '.mat'])
    else
        load([save_dir, subj, '/event_params.mat'])
    end
catch
    minITI = 40;
    minDuration = 50;
    event_thresh = 5000000;
end

events = find_events(pd(1,:), srate, minITI, minDuration, event_thresh, "show") ;
events = [events.onsets events.offsets];

%% Remove events that weren't marked well, or are demo
% here, only one bad trial at the end when the PD is unplugged

% if all trials are present, there should be 1000 at the end
events = events(11:end-1,:);

if ~sess_flag
    save([save_dir, subj, '/event_params.mat'], 'minITI', 'minDuration','event_thresh')
    save([save_dir, subj, '/events.mat'], 'events')
else
    save([save_dir, subj, '/event_param_sess', sess, '.mat'], 'minITI', 'minDuration','event_thresh')
    save([save_dir, subj, '/events_sess', sess, '.mat'], 'events')
end

%% Optional: combine sessions
sessions = ['1', '2'];

events = [];
if sess_flag
    for i = 1:numel(sessions)
        curr_events = load([save_dir, subj, '/events_sess', num2str(i), '.mat']);
        curr_events.events(:,3) = i;
        events = [events; curr_events.events];
    end
end

save([save_dir, subj, '/events.mat'], 'events')
