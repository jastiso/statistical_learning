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

subj = [{'1'}, {'2'}, {'3'}, {'4'}, {'5'}, {'6'}, {'7'}, {'8'}, {'10'}, {'12'}];
nSim = 100;
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
    [r_dir, '/sig_all', flag, '.node'],[r_dir, 'sig_all.mat'], ...
    ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/sig_all_',flag, '_', feat_type, '.jpg']);
BrainNet_MapCfg('/Users/stiso/Documents/MATLAB/BrainNetViewer_20171031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv',...
    [r_dir, '/sig_a_hat', flag, '.node'],[r_dir, 'sig_all.mat'], ...
    ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/sig_a_hat_', flag,'_', feat_type, '.jpg']);
BrainNet_MapCfg('/Users/stiso/Documents/MATLAB/BrainNetViewer_20171031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv',...
    [r_dir, '/sig_d', flag, '.node'],[r_dir, 'sig_all.mat'], ...
    ['/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_img/sig_d_', flag, '_', feat_type,'.jpg']);

%% Regions barplot

% for each subject, load regions and sig_idx
if strcmp(flag, 'all')
    all_reg = {};
    all_reg_d = {};
    subj_list = []; subjd = [];
    cnt = 1; cntd = 1;
    clear localization
    
    for s = 1:numel(subj)
        load([save_dir, subj{s}, '/header_clean.mat'], 'elec_labels')
        load([r_dir, 'subj' subj{s}, '/searchlight_corrs', '_', feat_type, '.mat'], 'sig_idx', 'sig_null_idx')
        
        % load localization
        loc_folder = dir([save_dir, subj{s}]);
        loc_folders = loc_folder([loc_folder(:).isdir]);
        loc_folders = loc_folders(~ismember({loc_folders(:).name},{'.','..'}));
        if exist([save_dir, subj{s}, '/', loc_folders.name, '/electrode_localization.csv'],'file')
            if s == 1
                localization = readtable([save_dir, subj{s}, '/', loc_folders.name, '/electrode_localization.csv'], 'TreatAsEmpty',{'NaN'});
                idx = cellfun(@(x) any(strcmp(x,elec_labels)), localization.electrode_name);
                % essentially we want to skip all numeric types
                col_idx = cellfun(@(x) any([contains(x, 'label'), contains(x,'coord'), contains(x,'electrode_name')]) &...
                    ~contains(x,'Random') & ~contains(x,'MMP_in_MNI') & ~contains(x, 'cc') & ~contains(x,'AAL600')...
                    & ~contains(x,'brodman'), localization.Properties.VariableNames);
                localization = localization(idx,col_idx);
                localization.subj = repmat(str2double(subj{s}), size(localization,1),1);
                curr_loc = localization;
                if size(localization,1) ~= numel(elec_labels)
                    disp('Some contacts arent in the localization')
                end
            else
                curr_loc = readtable([save_dir, subj{s}, '/', loc_folders.name, '/electrode_localization.csv'], 'TreatAsEmpty',{'NaN'});
                idx = cellfun(@(x) any(strcmp(x,elec_labels)), curr_loc.electrode_name);
                col_idx = cellfun(@(x) any([contains(x, 'label'), contains(x,'coord'), contains(x,'electrode_name')])&...
                    ~contains(x,'Random') & ~contains(x,'MMP_in_MNI') & ~contains(x,'cc') & ~contains(x,'AAL600')...
                    & ~contains(x, 'brodman'), curr_loc.Properties.VariableNames);
                curr_loc = curr_loc(idx,col_idx);
                curr_loc.subj = repmat(str2double(subj{s}), size(curr_loc,1),1);
                if size(curr_loc,1) ~= numel(elec_labels)
                    disp('Some contacts arent in the localization')
                end

                if any(strcmp('tissue_segmentation_distance_from_label_3', curr_loc.Properties.VariableNames))
                   curr_loc = removevars(curr_loc, 'tissue_segmentation_distance_from_label_3');
                end
                idx = cellfun(@(x) any(strcmp(x,localization.Properties.VariableNames)), curr_loc.Properties.VariableNames);
                curr_loc = removevars(curr_loc, curr_loc.Properties.VariableNames(~idx));
                localization = [localization; curr_loc];
            end
            
            %sig idx
            idx = cellfun(@(x) any(strcmp(x,elec_labels(sig_idx))), curr_loc.electrode_name);
            all_reg(cnt:(cnt+sum(sig_idx)-1),1) = table2cell(curr_loc(idx,'electrode_name'));
            null_idx = cellfun(@(x) any(strcmp(x,elec_labels(sig_null_idx))), curr_loc.electrode_name);
            all_reg_d(cntd:(cntd+sum(sig_null_idx)-1),1) = table2cell(curr_loc(null_idx,'electrode_name'));
            subj_list = [subj_list; repmat(str2double(subj{s}), sum(sig_idx),1)];
            subjd = [subjd; repmat(str2double(subj{s}), sum(sig_null_idx),1)];
            cnt = cnt + sum(sig_idx);
            cntd = cntd + sum(sig_null_idx);
            
        end
    end
    searchlight_regions = table([subj_list; subjd], [repmat('latent', numel(all_reg),1); ...
        repmat('euclid', numel(all_reg_d),1)], [all_reg; all_reg_d], 'VariableNames', {'subj', 'space', 'electrode_name'});
    writetable(searchlight_regions, [r_dir, 'searchlight_regions', '_', feat_type, '.csv'])
    writetable(localization, [r_dir, 'localization', '_', feat_type, '.csv']);
end
