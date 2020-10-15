%% Get the similarity of stimuli

% define constants
stim_dir = '/Users/stiso/Documents/Code/tasks/ECoG_graph_learn_blocks/img/';
nNode = 10;
prefix = 'target_';

%% Load all images
stim = zeros(214*935, nNode); % I got the dimensions from loading a sample image
for i = 1:nNode
    tmp = imread([stim_dir,prefix, num2str(i), '.png']);
    imshow(tmp); pause(0.1)
    tmp = sum(tmp,3);
    stim(:,i) = reshape(tmp,1,[])';
end

%% Get normalized euclidean distance

D = zeros(nNode);
for i = 1:nNode
    for j = i:nNode
        X = zscore(stim(:,i))';
        Y = zscore(stim(:,j))';
        d = pdist([Y;X], 'Euclidean');
        D(i,j) = d;
    end
end
D = D+D';
D(logical(eye(nNode))) = NaN;

% plot
figure(1); clf
imagesc(D); colorbar

save('/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_analysis/RSA_null.mat')
