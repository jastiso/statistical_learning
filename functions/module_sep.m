function [mod_sep] = module_sep(Y, mod_idx,n)
% get the module separation based on low dimensional data points
%   Inputs
% Y:        N x 2 matrix where N is the number of nodes. These are the
%           coordinates in soe low dimensional space
% mod_idx:  a vector of numerical or logical coordinates selecting the
%           current module
% n:        A scalar or logical vector indicating the position of the node
%           of interest

num = mean(sqrt((Y(mod_idx,1) - Y(n,1)).^2 + (Y(mod_idx,2) - Y(n,2)).^2));
denom = mean(pdist(Y));
mod_sep = num/denom;
end

