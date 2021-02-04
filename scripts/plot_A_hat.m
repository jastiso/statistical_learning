% Loop over beta values:
nSim = 6;
betas = [0.00000000000000001, 0.00001, .001, .1, .5, 1];
nNode = 10;
A_hats = zeros(nNode,nNode,nSim);
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
flag = 'mod';
if strcmp(flag, 'mod')
    A = M;
else
    A=L;
end
for j = 1:nSim
    A_hat = (1 - exp(-betas(j)))*A*(eye(nNode) - exp(-betas(j))*A)^(-1);
    A_hat = A_hat - diag(diag(A_hat));
    A_hats(:,:,j) = A_hat;
    imagesc(A_hat)
    pause(0.001)
    
end

save(['/Users/stiso/Documents/Code/graph_learning/ECoG_data/behavior_preprocessed/a_hat_ex_', flag, '.mat'],'A_hats')
