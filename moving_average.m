
function mvg = moving_average(sig, N)
% Moving average of signal
%
% Parameters
% ----------
%
% sig : array_like
%     Signal for moving average
% N: number of averaging samples
%
% Returns
% -------
%
% mvg : array_like
%     Average signal of length length(sig)-n+1

sig = sig(:); % convert to column vector
m = length(sig);
ret = cumsum([0; sig]); % compute the cumulative sum of the signal
mvg = (ret(N+1:end) - ret(1:end-N))/N; % compute the moving average by subtracting the previous N values
end