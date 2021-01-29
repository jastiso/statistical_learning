% Script to compute rate-distortion curves and other measurements for
% A_hat from different betas

clear

addpath(genpath('/Users/stiso/Documents/Code/graph_learning/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
% define variables
save_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/behavior_preprocessed/';
ephys_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_raw/';
img_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/behavior_preprocessed/images/';
r_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_analysis/';

subjs = [{'1'}, {'2'}, {'3'}, {'4'}, {'5'}, {'6'}, {'8'}, {'10'}, {'12'},{'18'}];
nNode = 10;
M = [0 1 1 1 0 0 0 0 0 1;
    1 0 1 1 1 0 0 0 0 0;
    1 1 0 1 1 0 0 0 0 0;
    1 1 1 0 1 0 0 0 0 0;
    0 1 1 1 0 1 0 0 0 0;
    0 0 0 0 1 0 1 1 1 0;
    0 0 0 0 0 1 0 1 1 1;
    0 0 0 0 0 1 1 0 1 1;
    0 0 0 0 0 1 1 1 0 1
    1 0 0 0 0 0 1 1 1 0].*0.25;
L = [0 1 1 0 0 0 0 0 1 1;
    1 0 1 1 0 0 0 0 0 1;
    1 1 0 1 1 0 0 0 0 0;
    0 1 1 0 1 1 0 0 0 0;
    0 0 1 1 0 1 1 0 0 0;
    0 0 0 1 1 0 1 1 0 0;
    0 0 0 0 1 1 0 1 1 0;
    0 0 0 0 0 1 1 0 1 1;
    1 0 0 0 0 0 1 1 0 1
    1 1 0 0 0 0 0 1 1 0].*0.25;
graphs = [{'mod'},{'lat'}];
nGraph = numel(graphs);

trans1 = [253,224,239]./255;
trans2 = [230,245,208]./255;
within1 = [233,163,201]./255;
within2 = [161,215,106]./255;
center1 = [197,27,125]./255;
center2 = [77,146,33]./255;
cluster_colors = [trans1; within1; center1; within1; trans1; trans2; within2; center2; within2; trans2];

nSim = 1000;
test_betas = (logspace(-3,3,1000));

norm = false;

%% get compressibility

% Structure to save:
measures_struct = [];

for i = 1:nGraph
    
    if strcmp(graphs{i},'mod')
        A = M;
        colors = cluster_colors;
    else
        A = L;
        colors = viridis(nNode);
    end
    
    % Different things to compute:
    rate_distortion_upper = cell(nSim, 1);
    rate_distortion_lower = cell(nSim, 1);
    compressibility = cell(nSim, 1);
    module_sep = cell(nSim,1);

    
    % Loop over beta values:
    for j = 1:nSim
        
        A_hat = (1 - exp(-test_betas(j)))*A*(eye(nNode) - exp(-test_betas(j))*A)^(-1);
        %A_hat = round(A_hat,4);
        ks = sum(A_hat ~= 0,1);
        E = numel(nonzeros(A_hat));
        A_hat_dist = A_hat - diag(diag(A_hat));
        
        % plot A_hat and get separability of modules
        [Y, E] = cmdscale(A_hat_dist,2);
        figure(2); clf
        scatter(Y(:,1), Y(:,2), 10000, colors, '.', 'MarkerFaceAlpha', 0.4)
        if strcmp(graphs{i},'mod')
            module_sep{j} = pdist([mean(Y(2:4,1)),mean(Y(2:4,2));mean(Y(7:9,1)),mean(Y(7:9,2))]);
        end
        
        % Use all pairs heuristic:
        [rd_upper, rd_lower, Cs, ~] = rate_distortion(A_hat, 1, 100);
        compressibility_temp = mean(rd_upper(end) - rd_upper);       

        % Record measures for this sample:
        rate_distortion_upper{j} =  rd_upper;
        rate_distortion_lower{j} = rd_lower;
        N = size(A_hat,1);

        [p_ss, D] = eigs(A_hat');

        [~, ind] = max(diag(D));

        p_ss = p_ss(:,ind)/sum(p_ss(:,ind));
    
        S = -sum(nonzeros(repmat(p_ss,1,N).*A_hat).*log2(nonzeros(A_hat)));
        
        if norm
            compressibility{j} = compressibility_temp/S;
        else
            compressibility{j} = compressibility_temp;
        end
    end


    % Add measures to structure:
    measures_struct.([graphs{i}, '_rate_distortion_upper']) = rate_distortion_upper;
    measures_struct.([graphs{i}, '_rate_distortion_lower']) = rate_distortion_lower;
    measures_struct.([graphs{i}, '_compressibility']) = compressibility;
    measures_struct.([graphs{i}, '_module_sep']) = module_sep;
end

% Save measures:
save([r_dir, graphs{i}, '_rate_distortion_samples'], '-struct', 'measures_struct');

%% plot

load(['/Users/stiso/Documents/Python/graphLearning/old_tasks/mTurk-10-node-breaks/data/preprocessed/', 'max_ent.mat'], 'beta')
beta_mturk = beta;
beta_mturk(beta_mturk > 1000) = 1000;
beta_mturk(beta_mturk < 0) = 0;
load([save_dir, 'max_ent.mat'], 'beta')
figure(1); clf
for b = 1:numel(beta)
    plot([log10(beta_mturk(b)), log10(beta_mturk(b))], ...
        [min([measures_struct.lat_compressibility{:}]),max([measures_struct.mod_compressibility{:}])], 'color', rgb('gray')); hold on
end
for b = 1:numel(beta)
    plot([log10(beta(b)), log10(beta(b))],...
        [min([measures_struct.lat_compressibility{:}]),max([measures_struct.mod_compressibility{:}])], 'r'); hold on
end
plot(log10(test_betas), [measures_struct.mod_compressibility{:}], 'color', rgb('darkorchid'), 'linewidth', 4);
plot(log10(test_betas), [measures_struct.lat_compressibility{:}], 'color', rgb('steelblue'), 'linewidth', 4);
xlabel('log10( beta )'); ylabel('Compressibility')
saveas(gca, [img_dir, 'compressibility.png']);


figure(2); clf
plot(log10(test_betas), [measures_struct.mod_module_sep{:}], 'linewidth', 2); hold on
for b = 1:numel(beta)
    if (mod(str2double(subjs{b}),2) == 0) && (str2double(subjs{b}) ~= 6)
        plot([log10(beta(b)), log10(beta(b))],...
            [min([measures_struct.mod_module_sep{:}]),max([measures_struct.mod_module_sep{:}])], 'r');
    end
end
