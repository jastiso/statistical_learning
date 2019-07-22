function [unaligned_elecs] = AAL_match(clinical)
% Match up AAL and clinical labels
%   basically you have to manually figure out when you are missing elecs in
%   the sequence. i.e. you have RB1-9, and 11-12, but not 10.

% initialize
unaligned_elecs = {};

% remove ROC and LOC if present
ref = cellfun(@(x) strcmpi(x,'roc') | strcmpi(x,'loc'), clinical);
clinical = clinical(~ref);

% get electrodes
prefix = cellfun(@(x) x(isstrprop(x,'alpha')),clinical, 'UniformOutput', 0);
groups = unique(prefix);
nGroup = numel(groups);

% contacts are typically groups of 8 or 12
cnt = 1;
for i = 1:nGroup
    curr = groups(i);
    group_idx = cellfun(@(x) strcmpi(curr,x), prefix);
    curr_elec = clinical(group_idx);
    last_elec = curr_elec{end};
    % assuming 2 positions for numbers
    if str2double(last_elec(end-1:end)) ~= sum(group_idx)
        unaligned_elecs{cnt} = curr{:};
        cnt = cnt + 1;
    end
end


end

