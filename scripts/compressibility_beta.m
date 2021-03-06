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
K=500; %number of random slelections of modules
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
    module_sep_mds = cell(nSim,1);
    module_sep_pca = cell(nSim,1);
    losses = cell(nSim,1);
    losses_pca = cell(nSim,1);
    
    % Loop over beta values:
    for j = 1:nSim
        null_sep = zeros(1,K);
        null_sep_pca = zeros(1,K);
        A_hat = (1 - exp(-test_betas(j)))*A*(eye(nNode) - exp(-test_betas(j))*A)^(-1);
        %A_hat = round(A_hat,4);
        ks = sum(A_hat ~= 0,1);
        E = numel(nonzeros(A_hat));
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
            
            mod1 = 2:4;
            mod2 = 7:9;
            curr_dist = 0;
            curr_dist_pca = 0;
            for n = 1:nNode
                mod_flag = any(n == mod1);
                
                if mod_flag
                    tmp_dist = module_sep(Y,mod1,n);
                    curr_dist = curr_dist + tmp_dist;
                    % pca
                    tmp_dist = module_sep(y_pca,mod1,n);
                    curr_dist_pca = curr_dist_pca + tmp_dist;
                elseif any(n == mod2)
                    tmp_dist = module_sep(Y,mod2,n);
                    curr_dist = curr_dist + tmp_dist;
                    % pca
                    tmp_dist = module_sep(y_pca,mod2,n);
                    curr_dist_pca = curr_dist_pca + tmp_dist;
                end
            end
            module_sep_mds{j} = curr_dist/numel([mod1,mod2]);
            module_sep_pca{j} = curr_dist_pca/numel([mod1,mod2]);
        else
            for u = 1:K
                nodes = datasample(1:nNode, 6, 'Replace', false);
                mod1 = nodes(1:3);
                mod2 = nodes(4:6);
                
                % check if these overlap
                if any(mod1 == mod2)
                    warning('Your set selection process is wrong')
                else
                    curr_dist = 0;
                    curr_dist_pca = 0;
                    for n = 1:nNode
                        mod_flag = any(n == mod1);
                        
                        if mod_flag
                            tmp_dist = module_sep(Y,mod1,n);
                            curr_dist = curr_dist + tmp_dist;
                            % pca
                            tmp_dist = module_sep(y_pca,mod1,n);
                            curr_dist_pca = curr_dist_pca + tmp_dist;
                        elseif any(n == mod2)
                            tmp_dist = module_sep(Y,mod2,n);
                            curr_dist = curr_dist + tmp_dist;
                            % pca
                            tmp_dist = module_sep(y_pca,mod2,n);
                            curr_dist_pca = curr_dist_pca + tmp_dist;
                        end
                    end
                    null_sep(u) = curr_dist/numel(nodes);
                    null_sep_pca(u) = curr_dist_pca/numel(nodes);
                end
            end
            module_sep_mds{j} = null_sep;
            module_sep_pca{j} = null_sep_pca;
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
    measures_struct.([graphs{i}, '_module_sep']) = module_sep_mds;
    measures_struct.([graphs{i}, '_module_sep_pca']) = module_sep_pca;
    measures_struct.([graphs{i}, 'losses']) = losses;
    measures_struct.([graphs{i}, 'losses_pca']) = losses_pca;
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
plot(log10(test_betas), [measures_struct.mod_compressibility{:}], 'color', [174/255,116/255,133/255], 'linewidth', 4);
plot(log10(test_betas), [measures_struct.lat_compressibility{:}], 'color', [101/255,111/255,147/255], 'linewidth', 4);
xlabel('log10( beta )'); ylabel('Compressibility')
saveas(gca, [img_dir, 'sim_compressibility.pdf']);

null_dists = reshape([measures_struct.lat_module_sep{:}],K,nSim);
null_dists_pca = reshape([measures_struct.lat_module_sep_pca{:}],K,nSim);
figure(2); clf
plot(log10(test_betas), null_dists, 'color', [101/255,111/255,147/255], 'linewidth', 0.1); hold on
plot(log10(test_betas), [measures_struct.mod_module_sep{:}], 'color', [174/255,116/255,133/255], 'linewidth', 4);

for b = 1:numel(beta)
    if (mod(str2double(subjs{b}),2) == 0) && (str2double(subjs{b}) ~= 6)
        plot([log10(beta(b)), log10(beta(b))],...
            [min([measures_struct.lat_module_sep{:}]),max([measures_struct.mod_module_sep{:}])], 'k');
    end
end
saveas(gca, [img_dir, 'sim_mod_dist.png']);

figure(2); clf
plot(log10(test_betas), null_dists_pca, 'color', [101/255,111/255,147/255], 'linewidth', 0.5); hold on
plot(log10(test_betas), [measures_struct.mod_module_sep_pca{:}], 'color', [174/255,116/255,133/255], 'linewidth', 4);
for b = 1:numel(beta)
    if (mod(str2double(subjs{b}),2) == 0) && (str2double(subjs{b}) ~= 6)
        plot([log10(beta(b)), log10(beta(b))],...
            [min([measures_struct.lat_module_sep_pca{:}])-0.2,max([measures_struct.mod_module_sep_pca{:}])+0.2], 'k');
    end
end
saveas(gca, [img_dir, 'sim_mod_dist_pca.png']);

figure(2); clf
plot(log10(test_betas), [measures_struct.modlosses{:}], 'color', [.4,.4,.4], 'linewidth', 3);
saveas(gca, [img_dir, 'linear_disc.png']);

figure(2); clf
plot(log10(test_betas), [measures_struct.modlosses_pca{:}], 'color', [.4,.4,.4], 'linewidth', 4);
saveas(gca, [img_dir, 'linear_disc_pca.png']);
