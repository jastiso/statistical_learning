%% Plot electrodes for each contrast
clear

addpath(genpath('/Users/stiso/Documents/MATLAB/Colormaps/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
subjs = [{'2'},{'4'}];

% define variables
save_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_raw/';
r_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_analysis/';

ext = '';

% initialize
x = [];
y = [];
z = [];
c = [];
s = [];
labels = {};

for i = 1:numel(subjs)
    
subj = subjs{i};
img_dir = ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/subj', subj];

% load stuff
load([r_dir, 'subj' subj, '/theta_peaks', ext, '.mat']); %ps tells you if it was >3std above 1/f
load([save_dir, subj, '/header_clean.mat']);
fname = dir([save_dir, subj, '/*/electrodenames_coordinates_mni.csv']);
coords = readtable([fname.folder, '/', fname.name]);
coords.Var1 = fix_names(coords.Var1);

if sum(ps) ~= 0
    coord_idx = ismember(cell2mat(coords.Var1),cell2mat(elec_labels(logical(ps))),'rows');
else
    coord_idx = false(size(cell2mat(coords.Var1)));
end

% concatenate
x = [x; coords.Var2(coord_idx)];
y = [y; coords.Var3(coord_idx)];
z = [z; coords.Var4(coord_idx)];
c = [c; peaks(logical(ps))];
s = [s; repmat(i, sum(ps),1)];
labels = vertcat(labels, AAL(logical(ps),1));

% save node file, with x, y, z, beta, p, AAL label
subj_node_file = [r_dir, 'subj', subj, '/theta_node', ext, '.node'];
write_bv_node( subj_node_file, coords.Var2(coord_idx), coords.Var3(coord_idx), coords.Var4(coord_idx)...
    , peaks(logical(ps)), repmat(i, sum(ps),1), AAL(logical(ps),1));

end

% save node file, with x, y, z, beta, p, AAL label
node_file = [r_dir, 'theta_node', ext,'.node'];
write_bv_node( node_file, x, y, z, c, s, labels);


%% Make Plot

BrainNet_MapCfg('/Users/stiso/Documents/MATLAB/BrainNetViewer_20171031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv',...
    [r_dir, 'theta_node', ext,'.node'],[r_dir, 'theta_inferno.mat'], ...
    '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/theta_all.jpg');

