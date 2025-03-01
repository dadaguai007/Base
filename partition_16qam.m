function [class1_mask, class2_mask] = partition_16qam(E)
    % Partition a 16-QAM signal into the inner and outer circles.

    % Separates a 16-QAM signal into the inner and outer rings, which have
    % different phase orientations. Detailed in [1].

    % Parameters
    % ----------
    %     E : electric field of the signal

    % Returns
    % -------
    %     class1_mask : array_like
    %         A mask designating the class 1 symbols which are the smallest and
    %         largest rings.
    %     class2_mask : array_like
    %         A mask designating the class 2 symbols which lie on the middle ring

    % References
    % ----------
    % [1] R. Muller and D. D. A. Mello, “Phase-offset estimation for
    %    joint-polarization phase-recovery in DP-16-QAM systems,” Photonics
    %    Technol. Lett. …, vol. 22, no. 20, pp. 1515–1517, 2010.

    S0 = cal_s0(E, 16);
    inner = (sqrt(S0 / 5) + sqrt(S0)) / 2;
    outer = (sqrt(9 * S0 / 5) + sqrt(S0)) / 2;
    Ea = abs(E);
    % 赋值，小于小圈，大于外圈
    class1_mask = (Ea < inner) | (Ea > outer);
    class2_mask = ~class1_mask;

end


function S0 = cal_s0(E, M)
    % Calculate the signal power S0 using the formula from Gao and Tepedelenlioglu (2005).

    % Parameters:
    %   E: Input field (array-like)
    %   M: Modulation order (int)
    % Returns:
    %   S0: Signal power estimate (float)

%     M=16;
    z = (0:M-1)';
    y = qammod(z,M);

    x=pnorm(y);
    scater=y./x;
    [C,~,ic] = unique(scater);
    a_counts = accumarray(ic,1);
    gamma=sum(C.^4.*a_counts/M);
    
    % 内圈
    r2 = mean(abs(E).^2);
    % 外圈
    r4 = mean(abs(E).^4);
    S1 = 1 - 2 * r2.^2 / r4 - sqrt((2 - gamma) * (2 * r2.^4 / r4.^2 - r2.^2 / r4));
    S2 = gamma * r2.^2 / r4 - 1;
    % S0 = r2/(1+S2/S1) because r2=S0+N and S1/S2=S0/N
    S0 = r2 / (1 + S2 / S1);
end


