function [beta, r0, r1, E, diff] = learn_linear_real_full(S, rt, trials)
% Inputs: Sequence S of nodes of length Ttot and sequence of reaction times 
% rt of length T. We assume T < Ttot and that rt contains reaction times to
% the trial indices 'trials' in S. S, rt, and trials are assumed column 
% vectors.
%
% Outputs: Using our free energy model, we learn the inverse temperature
% beta and the linear reaction time coefficients r0 and r1 by minimizing
% the RMSE between model predictions and the true reaction times in rt. In
% this function, beta, r0, r1, and RMSE are numbers that represent the
% values at the end of gradient descent estimation.
%
% NOTE: S cannot be longer than about 65,000.
%
% NOTE: We learn parameters based on the model r(t) = r0 + r1*b(t), where
% b(t) is the person's "belief" or anticipation at time t.
%
% NOTE: The difference between this function and 'learn_linear_real', is
% that here we search for the optimal beta parameter from which to start
% gradient descent.

% Other model parameters:
alpha = 1; % Ignore this parameter
diffThreshold = 10^(-6); % Precision with which we learn beta parameter

% Initialize beta parameter:
beta = learn_linear_search_func(S, rt, trials);

% Inertia for gradient descent (can usually ignore):
inertia = 0;

% Step size for gradient descent:
if abs(beta) == 0
    stepSize = .0000000001;
elseif abs(beta) < .001
    stepSize = .000001;
elseif abs(beta) < .005
    stepSize = .00001;
elseif abs(beta) < .01
    stepSize = .0001;
elseif abs(beta) < .05
    stepSize = .0005;
elseif abs(beta) < .1
    stepSize = .001;
elseif abs(beta) < .5
    stepSize = .01;
elseif abs(beta) < 1
    stepSize = .1;
elseif abs(beta) < 3
    stepSize = 5;
elseif abs(beta) <= 10
    stepSize = 100;
else
    stepSize = 100000;
end

% stepSize = stepSize*10; % Do this if algorithm converges too slowly
% stepSize = stepSize/10; % Do this if algorithm fluctuates wildly

% Trials to include in parameter estimation:
T_full = length(S) - min(trials) + 1;
indsOK = trials - min(trials) + 1;
T = length(indsOK);

% Learn parameters:
db_new = 0;

[B, dBdb] = belief_db_sequence2(S, T_full, beta, alpha);

R = [T sum(B(indsOK)); sum(B(indsOK)) sum(B(indsOK).^2)]\...
	[sum(rt); sum(rt.*B(indsOK))];
r0 = R(1);
r1 = R(2);

E = sqrt(sum((rt - r0 - r1*B(indsOK)).^2)/T);
dEdb = -r1/(T*E)*sum((rt - r0 - r1*B(indsOK)).*dBdb(indsOK));
    
diff = abs(dEdb);

if (beta == 1000) || (beta == 0)
    return;
end

t = 1;

% Only allow up to 400 gradient steps, but make sure there are at least 10
% gradient descent steps
while t < 400 && (t < 10 || diff > diffThreshold)
    
    % Make parameter step:
    db_old = db_new;
    
    %db_new = (stepSize_b/t)*((1-inertia)*dEdb + inertia*db_old);
    db_new = stepSize*((1-inertia)*dEdb + inertia*db_old);
    
    %beta = max(beta - db_new, 0);
    beta = beta - db_new;
    
    if beta < 0
        break;
    end
    
    % New beliefs:
    [B, dBdb] = belief_db_sequence2(S, T_full, beta, alpha);
    
    R = [T sum(B(indsOK)); sum(B(indsOK)) sum(B(indsOK).^2)]\...
        [sum(rt); sum(rt.*B(indsOK))];
    r0 = R(1);
    r1 = R(2);
    
    % Change in prediction error:
    E = sqrt(sum((rt - r0 - r1*B(indsOK)).^2)/T);
    
    % Calculate derivatives of RMSE wrt r0 and beta:
    dEdb = -r1/(T*E)*sum((rt - r0 - r1*B(indsOK)).*dBdb(indsOK));
    
    diff = abs(dEdb);
    
    % Print some things:
    t = t + 1
    beta
%     r0
%     r1
    %E
	dEdb
    
end