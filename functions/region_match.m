function [regions_aligned] = region_match(elec_labels, regions)
% Match up region and clinical labels, sometimes there are more included in
% the csv file of regions than the data
% 

% parse regions and labels from table
region_labels = table2cell(regions(:,1));
region_labels = cellfun(@(x) [x(1:2), sprintf('%02d', str2double(x(3:end)))], region_labels, 'UniformOutput', 0);

% remove ROC and LOC if present
ref = cellfun(@(x) strcmpi(x,'roc') | strcmpi(x,'loc'), elec_labels);
elec_labels = elec_labels(~ref);

% get index of regions not in elec labels
aligned = cellfun(@(x) contains(x, elec_labels), region_labels, 'UniformOutput', 1);
regions_aligned = table2cell(regions(aligned,2));



end

