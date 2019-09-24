function [B, dBdb] = belief_db_sequence2(S, T, beta, alpha)
% Inputs: Sequence (col vector) of nodes S of length Ttot, length of final
% transitions T for which we want beliefs, inverse temperature beta, and
% energy exponent alpha.
%
% Output: Tx1 vector of beliefs B, where B(t) represents the belief of
% transitioning from node L(t + Ttot - T - 1) to L(t+1 + Ttot - T - 1). We 
% also output Tx1 vector of derivatives of B with respect to inverse 
% temperature beta (dBdb).
%
% NOTE: The belief at time t (i.e. B(t)) only depends on the sequence of
% nodes that occured up to time t-1.
%
% NOTE: Only include sequences L of length less than about 65,000.
%
% NOTE: We also assume T < Ttot
%
% NOTE: The only difference between this function and 'belief_db_sequence'
% is that here we add 1 to every entry of the count matrix as a uniform
% prior

Ttot = length(S);
Tinit = Ttot - T - 1;

B = zeros(T,1);
dBdb = zeros(T,1);

N = max(S);

% Matrix of probabilities:

dt = repmat(1:Ttot, Ttot, 1);
dt = abs(dt' - dt);
LowT = tril(ones(Ttot, Ttot));

if beta >= 0
    E = dt.^(alpha);
else
    E = dt.^(alpha).*LowT;
    E = E - repmat(max(E,[],2), 1, Ttot);
    %E = E.*LowT;
end

%W = exp(-beta*E);
W = exp(-beta*E) + realmin;
%W = exp(-beta*E) + 10^(-10);

Z = sum(W.*LowT,2);

P = W.*LowT./repmat(Z, 1, Ttot);
dPdb = P.*(-E + repmat(sum(P.*E,2), 1, Ttot));

% Calculate beliefs:

% A = zeros(1,N);
% for i = 1:N
%     A(i) = length(unique(S(find(S(1:500) == i) + 1)));
% end

for t = 1:T
    
    ti = Tinit + t;

    i = S(ti);
    j = S(ti+1);
    
    %indsi = find(S(1:(ti-1)) == i);
    %Qi = sum(P(1:(ti-1),indsi),2);
    Qi = sum(P(1:(ti-1), S(1:(ti-1)) == i),2);
    dQidb = sum(dPdb(1:(ti-1), S(1:(ti-1)) == i),2);
    
    %indsj = find(S(2:ti) == j);
    %B(t) = sum(Qi(indsj))/sum(Qi);
    
%     B(t) = sum(Qi(S(2:ti) == j))/sum(Qi);
%     dBdb(t) = (sum(dQidb(S(2:ti) == j)) - B(t)*sum(dQidb))/sum(Qi);
%     B(t) = (sum(Qi(S(2:ti) == j)) + 1)/(sum(Qi) + A(i));
%     dBdb(t) = (sum(dQidb(S(2:ti) == j)) - B(t)*sum(dQidb))/(sum(Qi) + A(i));
    B(t) = (sum(Qi(S(2:ti) == j)) + 1)/(sum(Qi) + N);
    dBdb(t) = (sum(dQidb(S(2:ti) == j)) - B(t)*sum(dQidb))/(sum(Qi) + N);
    
end