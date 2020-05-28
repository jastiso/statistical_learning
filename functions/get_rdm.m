function [D,N] = get_rdm(feats, train, test, labels, nClass)
%get a representational dissiimilarity matrix using normalized euclidean
%distance
% Inputs:
% feats:    the features size N x M, where N is number of features and M is
%           number of observations
%
% train:    A logical vector size N x 1 that indexes the training data
%
% test:     A logical vector size N x 1 that indexed the testing data
%
% labels:   A vector of size N x 1 that indexes the classes of each
%           observation
%
% nClass:   A scalar that gives the number of classes

% Ouputs:
% D:        The RDM for this fold
% N:        The number of times each transition was seen

D = zeros(nClass);
N = zeros(nClass);
curr_node = labels(test);
for j = 1:nClass
    if j ~= curr_node
        X = feats(train,:);
        X = X(labels(train) == j,:);
        Y = feats(test,:);
        X = zscore(mean(X));
        Y = zscore(Y);
        d = pdist([Y;X], 'Euclidean');
        D(curr_node,j) = D(curr_node,j) + d;
        D(j,curr_node) = D(j,curr_node) + d;
        N(curr_node,j) = N(curr_node,j) + 1;
        N(j,curr_node) = N(j,curr_node) + 1;
    end
end

end

