function [filter_data] = filter_resp(filter_data, l, srate)
% make data for the given filter
step = [ones(1,l/2), ones(1,l/2)*2];
impulse = ones(1,l);
impulse(l/2) = 2;
% get new data strucutre
filter_data.trial = {[impulse; step]};
filter_data.label = [{'impulse'}, {'step'}];
filter_data.sampleinfo = [1*srate, l*srate];
filter_data.time = {0.001:0.001:l/1000};

end

