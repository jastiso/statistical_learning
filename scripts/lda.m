% Script to compute rate-distortion curves and other measurements for
% A_hat from different betas

clear

addpath(genpath('/Users/stiso/Documents/Code/graph_learning/'))
addpath(genpath('/Users/stiso/Documents/MATLAB/Colormaps/'))
addpath(('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/'))
% define variables
save_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/behavior_preprocessed/';
ephys_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_raw/';
img_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/behavior_preprocessed/images/';
r_dir = '/Users/stiso/Documents/Code/graph_learning/ECoG_data/ephys_analysis/';

subjs = [{'1'}, {'2'}, {'3'}, {'4'}, {'5'}, {'6'}, {'7'}, {'8'}, {'10'}, {'12'},{'18'}];
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
    1 0 0 0 0 0 1 1 0 1;
    1 1 0 0 0 0 0 1 1 0].*0.25;
graphs = [{'mod'},{'lat'}];
nGraph = numel(graphs);

trans1 = [239,169,186]./255;
trans2 = [177,191,146]./255;
within1 = [174,116,133]./255;
within2 = [126,138,96]./255;
cluster_colors = [trans1; within1; within1; within1; trans1; trans2; within2; within2; within2; trans2];
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
    losses = cell(nSim,1);
    losses_pca = cell(nSim,1);
    
    % Loop over beta values:
    for j = 1:nSim
        A_hat = (1 - exp(-test_betas(j)))*A*(eye(nNode) - exp(-test_betas(j))*A)^(-1);
        A_hat_dist = A_hat - diag(diag(A_hat));
        
        % plot A_hat and get separability of modules with two methods
        [Y, ~] = cmdscale(A_hat_dist,2);
        
        % now with pca
        coeff = pca(A_hat_dist);
        x = A_hat_dist - mean(A_hat_dist);
        y_pca = x*coeff;
        y_pca = y_pca(:,1:2);
        
        
        if strcmp(graphs{i},'mod')
            % discriminability
            ytab = table(Y(:,1), Y(:,2), 'VariableNames', [{'x1'}, {'x2'}]);
            fit = fitcdiscr(ytab, [{'1'},{'1'},{'1'},{'1'},{'1'},{'2'},{'2'},{'2'}, {'2'},{'2'}]);
            losses{j} = loss(fit, ytab, [{'1'},{'1'},{'1'},{'1'},{'1'},{'2'},{'2'},{'2'}, {'2'},{'2'}]);
            %pca
            ytab_pca = table(y_pca(:,1), y_pca(:,2), 'VariableNames', [{'x1'}, {'x2'}]);
            fit_pca = fitcdiscr(ytab_pca, [{'1'},{'1'},{'1'},{'1'},{'1'},{'2'},{'2'},{'2'}, {'2'},{'2'}]);
            losses_pca{j} = loss(fit_pca, ytab_pca, [{'1'},{'1'},{'1'},{'1'},{'1'},{'2'},{'2'},{'2'}, {'2'},{'2'}]);
        else
            
        end
    end
    
    
    % Add measures to structure:
    measures_struct.([graphs{i}, 'losses']) = losses;
    measures_struct.([graphs{i}, 'losses_pca']) = losses_pca;
end

% Save measures:
save([r_dir, graphs{i}, '_rate_distortion_samples'], '-struct', 'measures_struct');

