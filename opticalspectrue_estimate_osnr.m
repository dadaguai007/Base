function OSNRdB = opticalspectrue_estimate_osnr(E, lambda, f,resolution)
% 光谱函数中的估计OSNR
tolVec = [1e-3 1e-2 1e-1 1]; % tolerances for measuring noise floor

[YdBm, ly]=Opticalspectrum(E, lambda, f,resolution);

[power_sig_noise, idx] = max(YdBm); % signal + noise power

% Measure noise floor assuming different tolerances
for k = 1:length(tolVec)
    tol = tolVec(k);
    idx1 = find(abs(diff(YdBm(idx:-1:1))) < tol, 1); % Select points where difference was below tol
    idx2 = find(abs(diff(YdBm(idx:1:end))) < tol, 1); % Select points where difference was below tol
    if not(isempty(idx1)) && not(isempty(idx2))
        break
    end
end
if isempty(idx1) || isempty(idx2)
    error('OSA estimate_osnr fail: it was not possible to find noise floor')

end
%调制坐标位置
idx1 = idx - idx1(1);
idx2 = idx + idx2(1);
%noise power
NdB = interp1([ly(idx1), ly(idx2)], [YdBm(idx1), YdBm(idx2)], ly(idx));

% W
SN = 10^(power_sig_noise/10);
N = 10^(NdB/10);
% dB
OSNRdB = 10*log10(SN/N - 1); % unbias OSNR estimate