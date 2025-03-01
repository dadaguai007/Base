% %% ============================ Preemphasis ===============================
%预加重
function out=preemphasis(xd,N,SpS,Fs,preemphRange)
% xd 输入，N信号长度，Fs采样，preemphRange预加重范围，SpS上采样
%preemphRange = 25e9;
    femph = abs(freq_time(N*SpS, Fs));
    femph(femph >= preemphRange) = 0;
    %
    preemphasis_filter = 10.^(polyval([-0.0013 0.5846 0], femph/1e9)/20);  % Coefficients were measured in the lab  
% 
    out = real(ifft(fft(xd).*ifftshift(preemphasis_filter)));
    function [f, t] = freq_time(N, fs)
        %% Generate frequency and time measures
        dt = 1/fs;
        t = 0:dt:(N-1)*dt;
        df = 1/(dt*N);
        f = -fs/2:df:fs/2-df;
    end
end