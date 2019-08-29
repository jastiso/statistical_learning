%% Plot electrodes for each contrast
clear

addpath(genpath('/Users/stiso/Documents/MATLAB/IRASA/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
subjs = [{'2'},{'4'}];

% define variables
save_dir = '/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_raw/';
r_dir = '/Users/stiso/Documents/Python/graphLearning/ECoG data/ephys_analysis/';

ext = '_hg';

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
fname = dir([save_dir, subj, '/*/electrodenames_coordinates_mni.csv']);
coords = readtable([fname.folder, '/', fname.name]);
coords.Var1 = fix_names(coords.Var1);

% find significant elecs
sig_idx_ramp = logical(ramp.p < 0.05);
sig_idx_mod = logical(mod.p < 0.05);

mod = mod(sig_idx_mod,:);
ramp = ramp(sig_idx_ramp,:);

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

% concatenate
x_ramp = [x_ramp; coords.Var2(ramp_coord_idx)];
y_ramp = [y_ramp; coords.Var3(ramp_coord_idx)];
z_ramp = [z_ramp; coords.Var4(ramp_coord_idx)];
c_ramp = [c_ramp; ramp.betas];
s_ramp = [s_ramp; -log10(ramp.p)];
labels_ramp = vertcat(labels_ramp, ramp.region);

x_mod = [x_mod; coords.Var2(mod_coord_idx)];
y_mod = [y_mod; coords.Var3(mod_coord_idx)];
z_mod = [z_mod; coords.Var4(mod_coord_idx)];
c_mod = [c_mod; mod.betas];
s_mod = [s_mod; -log10(mod.p)];
labels_mod = vertcat(labels_mod, mod.region);

end

% save node file, with x, y, z, beta, p, AAL label
ramp_file = [r_dir, 'ramp_node', ext,'.node'];
write_bv_node( ramp_file, x_ramp, y_ramp, z_ramp, c_ramp, s_ramp, labels_ramp);

mod_file = [r_dir, 'mod_node', ext, '.node'];
write_bv_node( mod_file, x_mod, y_mod, z_mod, c_mod, s_mod, labels_mod);


