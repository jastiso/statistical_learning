%% Read raw data from ieeg.org

% make sure you have saved a password file before running. Save password
% file with:

%yourPath = IEEGSession.createPwdFile(?userName?,?yourPassword?);

% for ieeg toolbox documentation, see https://www.ieeg.org/main.html

addpath(genpath('/Users/stiso/Documents/MATLAB/ieeg-matlab-1.13.2/'))

% define variables
HUP_ID = 'HUP187';
subj = '2';
save_dir = '/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_raw/';


session = IEEGSession([HUP_ID, '_typingTask'], 'jastiso', '/Users/stiso/Documents/MATLAB/jas_ieeglogin.bin');
nElec = size(session.data.channelLabels,1);
% duration in ms
srate = session.data.sampleRate;
dur = (session.data.rawChannels(1).get_tsdetails.getDuration/1e6)*srate;
%st = 27578.11*srate;
st = (session.data.rawChannels(1).get_tsdetails.getStartTime/1e6)*srate;


% there are limits on how much data you can request at once.
% The current (hz * channels * seconds) max is (500 * 130 * 2000)
data = session.data.getvalues(st:dur,1:nElec);
data = data';

% get other useful stuff from header
labels = session.data.channelLabels(:,1); % only get the first row, these are not bipolar referenced yet
pd_idx = cellfun(@(x) contains(x,'DC'),labels);
elec_labels = labels(~pd_idx);
pd = data(pd_idx,:);
data = data(~pd_idx,:);
 
% save stuff
% make directories
if ~exist([save_dir, subj], 'dir')
    mkdir([save_dir, subj])
end

save([save_dir, subj, '/raw_data.mat'], 'data')
save([save_dir, subj, '/raw_pd.mat'], 'pd')
save([save_dir, subj, '/header.mat'], 'labels', 'srate', 'HUP_ID', 'subj')