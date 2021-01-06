function Ahats = Ahat_sequence(S, beta)
% Inputs: Sequence (col vector) of nodes S of length T, and inverse 
% temperature parameter beta.
%
% Output: NxNxT array of internal transition probability estimates
% Ahats, where Ahats(:,:,t) represents the transition probability estimates
% at time t.
%
% NOTE: The transition probability estimates at time t (i.e. Ahats(:,:,t))
% only depend on the sequence of nodes that occured up to time t-1.
%
% NOTE: This only works for sequences S of length less than about 65,000.

T = length(S);
N = max(S);

% Initialize:
Ahats = zeros(N,N,T);

% Compute matrix of shuffling probabilities:
dt = repmat(1:T, T, 1);
dt = abs(dt' - dt);
LowT = tril(ones(T));

W = exp(-beta*dt) + realmin;
Z = sum(W.*LowT,2);
P = W.*LowT./repmat(Z, 1, T);

% Compute transition probability estimates:
n = ones(N);
Ahats(:,:,1) = n./repmat(sum(n,2), 1, N);

% Loop over all steps in the sequence:
for t = 1:(T-1)

    j = S(t+1);
    
    % Belief that node i appeared at time t:
    Q = zeros(N,1);
    for i = 1:N
        Qi = sum(P(1:t, S(1:t) == i),2);
        Q(i) = Qi(end);
    end
    
    % Update transition counts:
    n(:,j) = n(:,j) + Q;
    
    % Compute transition probability estimates:
    Ahats(:,:,t+1) = n./repmat(sum(n,2), 1, N);
    
end

