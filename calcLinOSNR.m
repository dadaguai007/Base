function OSNR = calcLinOSNR(Ns, Pin, alpha, Ls, OSNRin, NF, Fc, Bref)
    % Calculate the OSNR evolution in a multi-span fiber transmission system.
% Ns 为光纤的段落
% Pin为输入光功率 dBm
% alpha损耗系数
% Ls每段的长度
% OSNRin表示输入OSNR值  dB
% NF EDFA的噪声因子 4.5
% Fc光频率 193.1e12
% Bref为OSNR的参考带宽 12.5e9

    h=6.62607004e-34;
    G = alpha * Ls;
    NF_lin = 10^(NF / 10);  % Amplifier noise figure (linear)
    G_lin = 10^(G / 10);  % Amplifier gain (linear)
    nsp = (G_lin * NF_lin - 1) / (2 * (G_lin - 1));

    % ASE noise power calculation:
    % Ref. Eq.(54) of R. -J. Essiambre,et al, "Capacity Limits of Optical Fiber
    % Networks," in Journal of Lightwave Technology, vol. 28, no. 4,
    % pp. 662-701, Feb.15, 2010, doi: 10.1109/JLT.2009.2039464.
    N_ase = (G_lin - 1) * nsp * h * Fc;
    P_ase = (2 * N_ase * Bref) / 1e-3;  % In mW

    P_ase_dBm = 10 * log10(P_ase);  % ASE power in dBm generated per EDFA
    
    %  pass 1st fiber， and then sent to the first fiber channel
    %  dBm - dB
    Pn_in_edfa = (Pin - OSNRin) - alpha * Ls;  % ASE power sent to the 1st EDFA
    
    OSNR = zeros(1, Ns + 1);
    OSNR(1) = OSNRin;

    % Calculate OSNR at the output of each span
    for spanN = 2:Ns + 1
        Pn_out_edfa = 10 * log10(10^((Pn_in_edfa + G) / 10) + 10^(P_ase_dBm / 10));
        % Total ASE power at the output of the spanN-th EDFA
       
        % power dBm - ASE dBm = OSNR dB
        OSNR(spanN) = Pin - Pn_out_edfa;  % Current OSNR
        Pn_in_edfa = Pn_out_edfa - alpha * Ls;  % ASE power sent to the next EDFA
    end
end


%
% function OSNR = calcLinOSNR(Ns, Pin, alpha, Ls, OSNRin, NF, Fc, Bref)
%     % Calculate the OSNR evolution in a multi-span fiber transmission system.
% 
%     if nargin < 6
%         NF = 4.5;
%     end
%     if nargin < 7
%         Fc = 193.1e12;
%     end
%     if nargin < 8
%         Bref = 12.5e9;
%     end
% 
%     G = alpha * Ls;
%     NF_lin = 10^(NF / 10);  % amplifier noise figure (linear)
%     G_lin = 10^(G / 10);    % amplifier gain (linear)
%     nsp = (G_lin * NF_lin - 1) / (2 * (G_lin - 1));
% 
%     % ASE noise power calculation:
%     % "Capacity Limits of Optical Fiber
%     % Networks," in Journal of Lightwave Technology, vol. 28, no. 4,
%     % pp. 662-701, Feb.15, 2010, doi: 10.1109/JLT.2009.2039464.
%     N_ase = (G_lin - 1) * nsp * const.h * Fc;
%     P_ase = (2 * N_ase * Bref) / 1e-3;  % in mW
% 
%     P_ase_dBm = 10 * log10(P_ase);  % ASE power in dBm generated per EDFA
% 
%     Pn_in_edfa = (Pin - OSNRin) - alpha * Ls;  % ASE power sent to the 1st EDFA
%     OSNR = zeros(1, Ns + 1);
%     OSNR(1) = OSNRin;
% 
%     % Calculate OSNR at the output of each span
%     for spanN = 2:Ns + 1
%         Pn_out_edfa = 10 * log10(10^((Pn_in_edfa + G) / 10) + 10^(P_ase_dBm / 10));  % Total ASE power at the output of the spanN-th EDFA
%         OSNR(spanN) = Pin - Pn_out_edfa;  % current OSNR
%         Pn_in_edfa = Pn_out_edfa - alpha * Ls;  % ASE power sent to the next EDFA
%     end
% end
