%% Get atlas labels for all subjects
clear

top_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/';
data_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_raw/';
subjs = [{'2'}, {'4'}, {'6'}];
RIDs = [{'RID522'}, {'RID454'}, {'RID529'}];

%Current labels
load([top_dir, 'AAL_regions_by_subject.mat'])

%% get new labels

for i = 1:numel(subjs)
   if ~isfield(data, RIDs{i}) 
       coord_dir = [data_dir, subjs{i}, '/', RIDs{i}, '/electrode_coordinates_mni.csv'];
       [regions,D,atlas] = mniCoord2Label(coord_dir,'AAL116');
       eval(['data.', RIDs{i}, ' = regions;']);
   end
end

save([top_dir, 'AAL_regions_by_subject.mat'], 'data')
