function [new_names] = fix_names(names)
% changes strings to add 0 padding before

% parse extension from string
prefix = cellfun(@(x) x(isstrprop(x,'alpha')),names, 'UniformOutput', 0);
suffix = cellfun(@(x) x(~isstrprop(x,'alpha')),names, 'UniformOutput', 0);

% add 2 digits
suffix = cellfun(@(x) num2str(str2double(x),'%02.f'), suffix, 'UniformOutput', false);

% stick back together
new_names = cell(size(names));
for i = 1:numel(new_names)
   new_names{i} = [prefix{i}, suffix{i}]; 
end
end

