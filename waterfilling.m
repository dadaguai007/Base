function p = waterfilling(P0,N0,lambda)
% P0：总可用功率（线性单位，非分贝单位）。
% N0：噪声功率谱密度（每个子载波的噪声功率）。
% lambda：信道增益（每个子载波的信道增益）。
% Bisection search for mu
mu_low = 0; % Initial low
mu_high = (P0 + sum(N0./lambda)); % Initial high

stop_threshold = 1e-7; % Stop threshold

% Iterate while low/high bounds are further than stop_threshold
% 使用二分法进行找到合适的mu值
while(abs(mu_low - mu_high) > stop_threshold)
    mu = (mu_low + mu_high) / 2; % Test value in the middle of low/high

    % Solve the power allocation
    p = 1/mu - N0./lambda; 
    p = max(p,0); % Consider only positive power allocation
    %检查功率约束
    % Test sum-power constraints
    if (sum(p) > P0) % Exceeds power limit => lower the upper bound
        mu_low = mu;
    else % Less than power limit => increase the lower bound
        mu_high = mu;
    end
end

    