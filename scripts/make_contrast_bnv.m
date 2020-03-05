%% Plot electrodes for each contrast
clear

addpath(genpath('/Users/stiso/Documents/MATLAB/BrainNetViewer_20171031/'))
addpath(genpath('/Users/stiso/Documents/MATLAB/Colormaps/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
subjs = [{'2'},{'4'},{'6'}];

% define variables
save_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_raw/';
r_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_analysis/';

ext = '';

% initialize
x_ramp = [];
y_ramp = [];
z_ramp = [];
c_ramp = [];
s_ramp = [];
labels_ramp = {};
x_mod = [];
y_mod = [];
z_mod = [];
c_mod = [];
s_mod = [];
labels_mod = {};
x_max_ent = [];
y_max_ent = [];
z_max_ent = [];
c_max_ent = [];
s_max_ent = [];
labels_max_ent = {};

for s = 1:numel(subjs)
    
subj = subjs{s};
img_dir = ['/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_img/subj', subj];

% make diractories
if ~exist(img_dir, 'dir')
    mkdir(img_dir);
end

% load stuff
ramp = readtable([r_dir, 'subj' subj, '/ramp_stats', ext, '.csv']);
ramp = ramp(:,2:end);
mod = readtable([r_dir, 'subj' subj, '/mod_stats', ext, '.csv']);
mod = mod(:,2:end);
max_ent = readtable([r_dir, 'subj' subj, '/max_ent_stats', ext, '.csv']);
max_ent = max_ent(:,2:end);
fname = dir([save_dir, subj, '/*/electrodenames_coordinates_mni.csv']);
coords = readtable([fname.folder, '/', fname.name]);
coords.Var1 = fix_names(coords.Var1);

% find significant elecs
sig_idx_ramp = logical(ramp.p < 0.05);
sig_idx_mod = logical(mod.p < 0.05);
sig_idx_max_ent = logical(max_ent.p < 0.05);

mod = mod(sig_idx_mod,:);
ramp = ramp(sig_idx_ramp,:);
max_ent = max_ent(sig_idx_max_ent,:);

% get the x y and z coordinates of the significant elecs
if size(mod,1) ~= 0
    mod_coord_idx = ismember(cell2mat(coords.Var1),cell2mat(mod.elecs),'rows');
else
    mod_coord_idx = false(size(cell2mat(coords.Var1),1),1);
end

if size(ramp,1) ~= 0
    ramp_coord_idx = ismember(cell2mat(coords.Var1),cell2mat(ramp.elecs),'rows');
else
    ramp_coord_idx = false(size(cell2mat(coords.Var1)));
end

if size(max_ent,1) ~= 0
    max_ent_coord_idx = ismember(cell2mat(coords.Var1),cell2mat(max_ent.elecs),'rows');
else
    max_ent_coord_idx = false(size(cell2mat(coords.Var1)));
end

% concatenate
x_ramp = [x_ramp; coords.Var2(ramp_coord_idx)];
y_ramp = [y_ramp; coords.Var3(ramp_coord_idx)];
z_ramp = [z_ramp; coords.Var4(ramp_coord_idx)];
c_ramp = [c_ramp; ramp.betas];
s_ramp = [s_ramp; -log10(ramp.p)];
curr_regions = cellfun(@(x) strsplit(x), ramp.region, 'UniformOutput', false);
labels_ramp = vertcat(labels_ramp, cellfun(@(x) x{2}, curr_regions, 'UniformOutput', false));

x_mod = [x_mod; coords.Var2(mod_coord_idx)];
y_mod = [y_mod; coords.Var3(mod_coord_idx)];
z_mod = [z_mod; coords.Var4(mod_coord_idx)];
c_mod = [c_mod; mod.betas];
s_mod = [s_mod; -log10(mod.p)];
curr_regions = cellfun(@(x) strsplit(x), mod.region, 'UniformOutput', false);
labels_mod = vertcat(labels_mod, cellfun(@(x) x{2}, curr_regions, 'UniformOutput', false));

x_max_ent = [x_max_ent; coords.Var2(max_ent_coord_idx)];
y_max_ent = [y_max_ent; coords.Var3(max_ent_coord_idx)];
z_max_ent = [z_max_ent; coords.Var4(max_ent_coord_idx)];
c_max_ent = [c_max_ent; max_ent.betas];
s_max_ent = [s_max_ent; -log10(max_ent.p)];
curr_regions = cellfun(@(x) strsplit(x), max_ent.region, 'UniformOutput', false);
labels_max_ent = vertcat(labels_max_ent, cellfun(@(x) x{2}, curr_regions, 'UniformOutput', false));

end

% save node file, with x, y, z, beta, p, AAL label
ramp_file = [r_dir, 'ramp_node', ext,'.node'];
write_bv_node( ramp_file, x_ramp, y_ramp, z_ramp, c_ramp, s_ramp, labels_ramp);

mod_file = [r_dir, 'mod_node', ext, '.node'];
write_bv_node( mod_file, x_mod, y_mod, z_mod, c_mod, s_mod, labels_mod);

max_ent_file = [r_dir, 'max_ent_node', ext, '.node'];
write_bv_node( max_ent_file, x_max_ent, y_max_ent, z_max_ent, c_max_ent, s_max_ent, labels_max_ent);

%% Make Plot

cont = 'ramp';

BrainNet_MapCfg('/Users/stiso/Documents/MATLAB/BrainNetViewer_20171031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv',...
    [r_dir, cont, '_node', ext,'.node'],[r_dir, 'label_contrast.mat'], ...
    ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/', cont, ext, '.jpg']);

