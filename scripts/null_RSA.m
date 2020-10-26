%% Get the similarity of stimuli
clear
% define constants
stim_dir = '/Users/stiso/Documents/Code/tasks/ECoG_graph_learn_blocks/img/';
nNode = 10;
prefix = 'target_';

%% Load all images
stim = zeros(214*935, nNode); % I got the dimensions from loading a sample image
stim_sq = [];
for i = 1:nNode
    tmp = imread([stim_dir,prefix, num2str(i), '.png']);
    imshow(tmp); pause(0.1)
    tmp = sum(tmp,3);
    stim(:,i) = reshape(tmp,1,[])';
    stim_sq(:,:,i) = tmp;
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

%% get distance from center

x = size(stim_sq,1);
y = size(stim_sq,2);
ctr = [(x/2),(y/2)];

null_dist = zeros(nNode);
d = zeros(nNode,2);
for i = 1:nNode
   % get center point of current square
   [row,col] = find(stim_sq(:,:,i)==255);
   sq_x = mean(row);
   sq_y = mean(col);
   tmp = stim_sq(:,:,i);
   figure(i); clf
   tmp(round(sq_x),round(sq_y)) = 1000;
   imagesc(tmp)
   pause(0.1)
   
   d(i,:) = [sq_x,sq_y];
end

close all

for i = 1:nNode
    for j = (i+1):nNode
        D(i,j) = pdist([d(i,:);d(j,:)], 'euclidean');
    end
end
D=D+D';
imagesc(D)
save('/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_analysis/RSA_dist_null.mat')

