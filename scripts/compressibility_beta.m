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
nSubj = numel(subjs);
nNode = 10;
beta = zeros(nSubj,1);
r0 = zeros(nSubj,1);
r1 = zeros(nSubj,1);
E = zeros(nSubj,1);
diff = zeros(nSubj,1);
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

nSim = 1000;
test_betas = (logspace(-3,3,1000));


%% get compressibility

% Structure to save:
measures_struct = [];

for i = 1:nGraph
    
    if strcmp(graphs{i},'mod')
        A = M;
    else
        A = L;
    end
    
    % Different things to compute:
    rate_distortion_upper = cell(nSim, 1);
    rate_distortion_lower = cell(nSim, 1);
    compressibility = cell(nSim, 1);
    

    
    % Loop over beta values:
    for j = 1:nSim
        
        A_hat = (1 - exp(-test_betas(j)))*A*(eye(nNode) - exp(-test_betas(j))*A)^(-1);
        
        ks = sum(A_hat,1);

        % Use all pairs heuristic:
        [rd_upper, rd_lower, Cs, ~] = rate_distortion(A_hat, 1, 100);
        compressibility_temp = mean(rd_upper(end) - rd_upper);       

        % Record measures for this sample:
        rate_distortion_upper{j} =  rd_upper;
        rate_distortion_lower{j} = rd_lower;
        compressibility{j} = compressibility_temp;
    end


    % Add measures to structure:
    measures_struct.([graphs{i}, '_rate_distortion_upper']) = rate_distortion_upper;
    measures_struct.([graphs{i}, '_rate_distortion_lower']) = rate_distortion_lower;
    measures_struct.([graphs{i}, '_compressibility']) = compressibility;
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
    plot([log10(beta_mturk(b)), log10(beta_mturk(b))], [0.6,1.5], 'color', rgb('gray')); hold on
end
for b = 1:numel(beta)
    plot([log10(beta(b)), log10(beta(b))], [0.6,1.5], 'r'); hold on
end
plot(log10(test_betas), [measures_struct.mod_compressibility{:}], 'color', rgb('darkorchid'), 'linewidth', 4);
plot(log10(test_betas), [measures_struct.lat_compressibility{:}], 'color', rgb('steelblue'), 'linewidth', 4);
xlabel('log10( beta )'); ylabel('Compressibility')
saveas(gca, [img_dir, 'compressibility.png']);

