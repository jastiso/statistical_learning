function beta = learn_linear_search_func(S, rt, trials)
% Inputs: Sequence S of nodes of length Ttot and sequence of reaction times 
% rt of length T. We assume T < Ttot and that rt contains reaction times to
% the trial indices 'trials' in S. S, rt, and trials are assumed column 
% vectors.
%
% Output: Calculate RMSE for different values of beta and return the beta
% that gives the lowest error in predicting people's reaction times.

% Betas to sweep over in search:
betas = [0, logspace(-4, 1, 100), 1000];

% Trials to include in parameter estimation:
T_full = length(S) - min(trials) + 1;
indsOK = trials - min(trials) + 1;
T = length(indsOK);

% Loop over different betas:

RMSEs = zeros(1, length(betas));
    
for j = 1:length(betas)
        
    beta = betas(j);
    
    [B, ~] = belief_db_sequence2(S, T_full, beta, 1);
    
    R = [T sum(B(indsOK)); sum(B(indsOK)) sum(B(indsOK).^2)]\...
        [sum(rt); sum(rt.*B(indsOK))];
    r0 = R(1);
    r1 = R(2);
    
    E = sqrt(sum((rt - r0 - r1*B(indsOK)).^2)/T);
    
    RMSEs(j) = E;

end
    
[~, mInd] = min(RMSEs);
beta = betas(mInd);
    
end