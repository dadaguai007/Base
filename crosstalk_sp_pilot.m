function sig_XT = crosstalk_sp_pilot(sig_h,stepnum,stepsize,hxt,XT_ave,Snm,omega,omega_p)
% call format: sig_XT = crosstalk_new(sig,stepnum,stepsize,fs,K,Snm,omega)
% from test_crosstalk_3.m
% function: adding crosstalk and intercore group delay
% stepnum & stepsize(m): stepnum & stepsize of crosstalk calculation
% hxt: average distance between PMPs
% XT_ave: average crosstalk value between two channels
% K: average coupling strength at each PMP, K = sqrt(XT_ave * hxt)
% e.g. XT_ave = -10dB/100km = -30dB/km(power), hxt = 1e-3km, then K = 1e-3
% Snm(m) : group velocity difference (typical:1.35ps/m)
% omega: angular frequency of electrical signal

[leng,ch2] = size(sig_h);
% sig_f = fft(sig);
% sig_f = sig_f.';

sig = [sig_h(:, 1), sig_h(:, 3)];
pilot = [sig_h(:, 2), sig_h(:, 4)];
sig = sig.';
pilot = pilot.';

ch = ch2 / 2;

K = sqrt(10^(XT_ave/10) * hxt * 1e-3); 

I = eye(2);
Z = zeros(2,2);

switch ch
    case 2
        dim = 1;
    case 4
        dim = 4;
end

phase_total = zeros(stepnum, 1);
phase_total(1) = exp(-1j*2*pi*rand);

for z1 = 1:stepnum
%     % case A: original(incorrect串扰太快且每次都完全随机）
%     phase_total(1) = 0;
%     for k = 1 : stepsize / hxt
%         zk = stepsize*rand;
%         del_Phi = 2*pi*rand(1, dim);
%         phase_total(z1) = phase_total(z1) + exp(-1j*del_Phi).*exp(-1j*Snm*zk*omega.');
%     end
    
    % case B:Gaussian random walk (changing rate controlled by sigma)
    sigma = sqrt(stepsize / hxt * K^2 / 2);
    phase_total(z1+1) = sigma*(randn(1, dim) + 1j*randn(1, dim)) + phase_total(z1);
    zk = stepsize * rand; % !!!

    switch ch
        case 2
            for n = 1:leng
                C12 = -1j * abs(K) * phase_total(z1) * exp(-1j*Snm*zk*omega(n)); 
                C21 = -conj(C12);
                C{z1} = [sqrt(1-abs(C12)^2) C12; C21 sqrt(1-abs(C12)^2)];
    %             sig_f = C{z1} * sig_f; 
                sig(:, n) = C{z1} * sig(:, n);
            end
            C12_p = -1j * abs(K) * phase_total(z1) * exp(-1j*Snm*zk*omega_p); 
            C21_p = -conj(C12_p);
            C_p{z1} = [sqrt(1-abs(C12_p)^2) C12_p; C21_p sqrt(1-abs(C12_p)^2)];
            pilot = C_p{z1} * pilot;
    end
%     disp(C{z1});
    h12(z1, 1) = C12;
    h21(z1, 1) = C21;
end

assignin('base', 'C', C);
assignin('base', 'h12', h12);
assignin('base', 'h21', h21);

% sig_f = sig_f.';
% sig_XT = ifft(sig_f);
sig = sig.';
pilot = pilot.';

sig_XT = [sig(:, 1), pilot(:, 1), sig(:, 2), pilot(:, 2)];