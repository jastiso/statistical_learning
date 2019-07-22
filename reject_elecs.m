function [rmv] = reject_elecs(data,thr)
%identify bad elecs by kurtosis
%   this is good for elecs with large artifacts, but not for persistently
%   noisy elecs
% Inputs:
%       thr         the threshold for rejection in z-scores
%       data        data, nElec x nTime
% Outputs
%       rmv         logical index of elecs to remove
% @author JStiso 07/2019


k = kurtosis(data');                % kurtosis of ts -- are there weird peaks? if so, maybe artifact?
zk = zscore(k);                     % zscore kurtosis

rmv = zk >= thr;                    % remove channels that show more kurtosis
end

