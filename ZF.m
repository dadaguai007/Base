 %% ZF 
    % 零强迫（ZF）频域均衡
function S_zf=ZF(H,ycp_mat,Nfft)
    %FDE
    % 频域均衡器系数计算
    W_zf = 1./H; % 零强迫频域均衡器的系数，公式为W_zf = 1 / H

    %DFT
    % 离散傅里叶变换
    ycp_fft = fft(ycp_mat,Nfft); % 对去除循环前缀后的信号进行Nfft点的快速傅里叶变换，转换到频域

    %Egalisation
    % 均衡操作
    S_zf = diag(W_zf)*ycp_fft; % 使用零强迫均衡器系数对频域信号进行均衡处理
end