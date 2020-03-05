%% Read raw data from ieeg.org

clear
clc

%% 
% make sure you have saved a password file before running. Save password
% file with:

%yourPath = IEEGSession.createPwdFile(?userName?,?yourPassword?);

% for ieeg toolbox documentation, see https://www.ieeg.org/main.html

addpath(genpath('/Users/stiso/Documents/MATLAB/ieeg-matlab-1.13.2/'))

% define variables
HUP_ID = 'HUP198';
subj = '10';
sess = '';
save_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_raw/';

% save stuff
% make directories
if ~exist([save_dir, subj], 'dir')
    mkdir([save_dir, subj])
end

if ~isempty(sess)
    session = IEEGSession([HUP_ID, '_typingTask_session', sess], 'jastiso', '/Users/stiso/Documents/MATLAB/jas_ieeglogin.bin');
else
    session = IEEGSession([HUP_ID, '_typingTask'], 'jastiso', '/Users/stiso/Documents/MATLAB/jas_ieeglogin.bin');
end
nElec = size(session.data.channelLabels,1);
% duration in ms
srate = session.data.sampleRate;
dur = (session.data.rawChannels(1).get_tsdetails.getDuration/1e6)*srate;
save_flag = 0; % save start time or not

try
    if ~isempty(sess)
        load([save_dir, subj, '/sess_', sess, 'start_time.mat'])
    else
        load([save_dir, subj, '/start_time.mat'])
    end
catch
    save_flag = 1;
    error('You dont have a start time saved for this subject! Go into iEEG and find it')
    % if you haven't saved the number, you have to go into iEEG portal and
    % find it
end
% this should work but hasnt been...
%st = (session.data.rawChannels(1).get_tsdetails.getStartTime/1e6)*srate;

%% Save start time
% only run this if you get an error above

if save_flag
    st = 19939.42*srate;
    if ~isempty(sess)
        save([save_dir, subj, '/sess_', sess, 'start_time.mat'], 'st')
    else
        save([save_dir, subj, '/start_time.mat'], 'st')
    end
end

%% Back to preprocessing

% there are limits on how much data you can request at once.
% The current (hz * channels * seconds) max is (500 * 130 * 2000)
data1 = session.data.getvalues(st:dur,1:nElec/4);
data1 = data1';

data2 = session.data.getvalues(st:dur,(nElec/4 + 1):(nElec/4 + nElec/4));
data2 = data2';

data3 = session.data.getvalues(st:dur,(2*nElec/4 + 1):(2*nElec/4 + nElec/4));
data3 = data3';

data4 = session.data.getvalues(st:dur,(3*nElec/4 + 1):nElec);
data4 = data4';

data = [data1; data2; data3; data4];

clear data1 data2 data3 data4
%% Save raw data before individual preprocessing

% get other useful stuff from header
labels = session.data.channelLabels(:,1); % only get the first row, these are not bipolar referenced yet
pd_idx = cellfun(@(x) contains(x,'DC'),labels);
% get only iEEG elecs
scalp_eeg_idx = ~cellfun(@(x) strcmpi(x(1),'r') | strcmpi(x(1),'l') ,labels);
ref = cellfun(@(x) strcmpi(x,'roc') | strcmpi(x,'loc'), labels);
rm_idx = pd_idx | scalp_eeg_idx | ref;
elec_labels = labels(~rm_idx);
pd = data(pd_idx,:);
data = data(~rm_idx,:);

pd = pd(:,~(isnan(pd(1,:))));
data = data(:,~(isnan(data(1,:))));

%check that pd and data are the same length
if size(pd,2) ~= size(data,2)
    warning('The length of your PD and data are not the same...something went wrong')
else
    if ~isempty(sess)
        save([save_dir, subj, '/raw_data_sess', sess, '.mat'], 'data', '-v7.3')
        save([save_dir, subj, '/raw_pd_sess', sess, '.mat'], 'pd')
        save([save_dir, subj, '/header_sess', sess, '.mat'], 'elec_labels', 'srate', 'HUP_ID', 'subj')
    else
        save([save_dir, subj, '/raw_data.mat'], 'data', '-v7.3')
        save([save_dir, subj, '/raw_pd.mat'], 'pd')
        save([save_dir, subj, '/header.mat'], 'elec_labels', 'srate', 'HUP_ID', 'subj')
    end
end

