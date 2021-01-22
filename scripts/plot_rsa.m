%% Plot combined searchlight results

clear
clc

addpath(genpath('/Users/stiso/Documents/MATLAB/IRASA/'))
addpath(genpath('/Users/stiso/Documents/MATLAB/Colormaps/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
addpath(genpath('/Users/stiso/Documents/MATLAB/BrainNetViewer_20171031/'))
addpath(genpath('/Users/stiso/Documents/Code/graph_learning/functions/'))
save_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_raw/';
r_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_analysis/';

subj = [{'1'}, {'2'}, {'3'}, {'4'}, {'5'}, {'6'}, {'8'}, {'10'}, {'12'}];
nSim = 100;
feat_type = 'lfp';
flag = 'all';
if strcmp(flag,'mod')
    idx = cellfun(@(x) mod(str2double(x),2) == 0, subj);
    subj = subj(idx);
elseif strcmp(flag,'lat')
    idx = cellfun(@(x) mod(str2double(x),2) == 1, subj);
    subj = subj(idx);
end

%%
all_a = [];
all_d = [];
for i = 1:numel(subj)
   subj_node_file = [r_dir, 'subj', subj{i}, '/A_hat_', feat_type, '.node'];
   f = fopen(subj_node_file);
   a = textscan(f, '%10f %10f %d %d %d %s\n');
   all_a = [all_a; a];
   fclose(f);
   
   subj_node_file = [r_dir, 'subj', subj{i}, '/D_', feat_type, '.node'];
   f = fopen(subj_node_file);
   d = textscan(f, '%10f %10f %d %d %d %s\n');
   all_d = [all_d; d];
   fclose(f);
end

% reformat
x = [];xd = [];
y = []; yd = [];
z = []; zd = [];
labels = {}; labelsd = {};
cnt = 1; cntd = 1;
for i = 1:numel(subj)
   x = [x;all_a{i,1}];
   y = [y;all_a{i,2}];
   z = [z;all_a{i,3}];
   labels(cnt:(cnt+numel(all_a{i,6})-1)) = all_a{i,6};
   cnt = cnt + numel(all_a{i,6});
   
   xd = [xd;all_d{i,1}];
   yd = [yd;all_d{i,2}];
   zd = [zd;all_d{i,3}];
   labelsd(cntd:(cntd+numel(all_d{i,6})-1)) = all_d{i,6};
   cntd = cntd + numel(all_d{i,6});
end
labels = labels';
labelsd = labelsd';
labels_all = [labels;labelsd];
write_bv_node( [r_dir, 'sig_a_hat', flag, '.node'], x, y, z, [], [], labels);
write_bv_node( [r_dir, 'sig_d', flag, '.node'], xd, yd, zd, [], [], labelsd);
write_bv_node( [r_dir, 'sig_all', flag, '.node'], [x;xd], [y;yd], [z;zd], [ones(numel(x),1);ones(numel(xd),1)+1], [], labels_all);

%%

BrainNet_MapCfg('/Users/stiso/Documents/MATLAB/BrainNetViewer_20171031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv',...
        [r_dir, '/sig_all', flag, '.node'],[r_dir, 'sig_all.mat'], ...
        ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/sig_all_',flag,'.jpg']);
BrainNet_MapCfg('/Users/stiso/Documents/MATLAB/BrainNetViewer_20171031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv',...
        [r_dir, '/sig_a_hat', flag, '.node'],[r_dir, 'sig_all.mat'], ...
        ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/sig_a_hat_', flag, '.jpg']);
BrainNet_MapCfg('/Users/stiso/Documents/MATLAB/BrainNetViewer_20171031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv',...
        [r_dir, '/sig_d', flag, '.node'],[r_dir, 'sig_all.mat'], ...
        ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/sig_d_', flag, '.jpg']);
       
%% Regions barplot

% for each subject, load regions and sig_idx
all_reg = {};
all_reg_d = {};
subj_list = []; subjd = [];
cnt = 1; cntd = 1;
for s = 1:numel(subj)
    load([save_dir, subj{s}, '/header_clean.mat'], 'regions')
    load([r_dir, 'subj' subj{s}, '/searchlight_corrs.mat'], 'sig_idx', 'sig_null_idx')

    all_reg(cnt:(cnt+sum(sig_idx)-1),1) = regions(sig_idx);
    all_reg_d(cntd:(cntd+sum(sig_null_idx)-1),1) = regions(sig_null_idx);
    subj_list = [subj_list; repmat(str2double(subj{s}), sum(sig_idx),1)];
    subjd = [subjd; repmat(str2double(subj{s}), sum(sig_null_idx),1)];
    cnt = cnt + sum(sig_idx);
    cntd = cntd + sum(sig_null_idx);
end

searchlight_regions = table([subj_list; subjd], [repmat('latent', numel(all_reg),1); ...
    repmat('euclid', numel(all_reg_d),1)], [all_reg; all_reg_d], 'VariableNames', {'subj', 'space', 'region'});
writetable(searchlight_regions, [r_dir, 'searchlight_regions.csv'])
   