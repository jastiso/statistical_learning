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

subj = [{'1'}, {'2'}, {'3'}, {'4'}, {'8'}, {'10'}];
feat_type = 'lfp';
flag = 'mod';
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
        [r_dir, '/sig_all.node'],[r_dir, 'sig_all.mat'], ...
        ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/sig_all_.jpg']);
BrainNet_MapCfg('/Users/stiso/Documents/MATLAB/BrainNetViewer_20171031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv',...
        [r_dir, '/sig_a_hat', flag, '.node'],[r_dir, 'sig_all.mat'], ...
        ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/sig_a_hat_', flag, '.jpg']);
BrainNet_MapCfg('/Users/stiso/Documents/MATLAB/BrainNetViewer_20171031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv',...
        [r_dir, '/sig_d', flag, '.node'],[r_dir, 'sig_all.mat'], ...
        ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/sig_d_', flag, '.jpg']);
       
%% Repeat fo anticorrelations

%save node files
for i = 1:numel(subj)
    load([r_dir, 'subj' subj, '/searchlight_corrs.mat'], 'sig_idx', 'sig_null_idx')
    subj_node_file = [r_dir, 'subj', subj, '/A_hat_', feat_type, '_anti.node'];
    write_bv_node( subj_node_file, coords.Var2(~sig_idx), coords.Var3(~sig_idx), coords.Var4(~sig_idx),...
        A_hat_corr(~sig_idx), [], elec_labels(~sig_idx));
    
    subj_node_file = [r_dir, 'subj', subj, '/D_', feat_type, '_anti.node'];
    write_bv_node( subj_node_file, coords.Var2(~sig_null_idx), coords.Var3(~sig_null_idx), coords.Var4(~sig_null_idx)...
        , D_corr(~sig_null_idx), [], elec_labels(~sig_null_idx));
end

% combine node_files
all_a = [];
all_d = [];
for i = 1:numel(subj)
   subj_node_file = [r_dir, 'subj', subj{i}, '/A_hat_', feat_type, '_anti.node'];
   f = fopen(subj_node_file);
   a = textscan(f, '%10f %10f %d %d %d %s\n');
   all_a = [all_a; a];
   fclose(f);
   
   subj_node_file = [r_dir, 'subj', subj{i}, '/D_', feat_type, '_anti.node'];
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
write_bv_node( [r_dir, 'sig_a_hat', flag, '_anti.node'], x, y, z, [], [], labels);
write_bv_node( [r_dir, 'sig_d', flag, '_anti.node'], xd, yd, zd, [], [], labelsd);
write_bv_node( [r_dir, 'sig_all', flag, '_anti.node'], [x;xd], [y;yd], [z;zd], [ones(numel(x),1);ones(numel(xd),1)+1], [], labels_all);

%%

BrainNet_MapCfg('/Users/stiso/Documents/MATLAB/BrainNetViewer_20171031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv',...
        [r_dir, '/sig_all_', flag, 'anti.node'],[r_dir, 'sig_all.mat'], ...
        ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/sig_all_', flag, '_anti.jpg']);
BrainNet_MapCfg('/Users/stiso/Documents/MATLAB/BrainNetViewer_20171031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv',...
        [r_dir, '/sig_a_hat', flag, '_anti.node'],[r_dir, 'sig_all.mat'], ...
        ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/sig_a_hat_', flag, '_anti.jpg']);
BrainNet_MapCfg('/Users/stiso/Documents/MATLAB/BrainNetViewer_20171031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv',...
        [r_dir, '/sig_d', flag, '_anti.node'],[r_dir, 'sig_all.mat'], ...
        ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/sig_d_', flag, '_anti.jpg']);
       
