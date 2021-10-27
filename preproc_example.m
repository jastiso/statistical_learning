%% Preprocessing for subj X, HUPXXX

% I use this toolbox fo rplotting, but you can just comment out the plots
% if you want
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current/'))

%% Remove noisy elecs

%data-driven: find channels with large kurtosis or linelength
rmv = reject_elecs(data, 2); % thr of 2 stds

data = data(~rmv,:);
elec_labels = elec_labels(~rmv,:);

% data should look clean, with the exception of 60 Hz noise
eegplot(data, 'srate', srate) % srate is the samplig rate

%% Look at line noise

% filter out 60 Hz harmonics
% first input is order, then bounds of filter, and type of filter
[b,a] = butter(4, [59/(srate/2), 61/(srate/2)], 'stop');
data = filtfilt(b,a,data')'; % filtfilt prevents phase distortion

[b,a] = butter(4, [119/(srate/2), 121/(srate/2)], 'stop');
data = filtfilt(b,a,data')';

[b,a] = butter(4, [179/(srate/2), 181/(srate/2)], 'stop');
data = filtfilt(b,a,data')';

% check for things outside of 60Hz...especially in theta/alpha range
% should look like a decreasing exponential with notches
spectopo(data, 0, srate)

%% CAR - common average reference
% if data looks to have different levels of noise based on elec, might want
% to do this in groups

data = data - mean(data,2); % demean
data = detrend(data')'; % as opposed to low pass filtering
data = get_CAR(data, elec_labels); % CAR by group - it expects elec labels to be strings, where the letters define the electrode and the number defines the contact

%% Get into fieldtrip format
% last argument is for events structure
ft_data = fieldtrip_format(data, srate, elec_labels, []);

