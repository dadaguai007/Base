clear;
close all;
clc;

% constant
hp = 6.6260657e-34;                        % [J*s]
c0 = 299792458;                            % [m/s] 
lambda = 1550*1e-9;                        % [m]
CenterFrequency = c0/lambda;               % [Hz]
rng(1234);
M = 4;
lengthSequence = 2^15; 
sps = 2 ;
roll_off = 0.1;
% RSOP set
% r = 4e-3;
% psi = rand*pi/2;
% chi = rand*pi/2;
% rate = 1e10;                             % 符号速率
% fs = rate*sps;  

% CMA set
mu_CMA = 6e-3;
ntaps = 5;
% DD-LMS
mu_DDLMS = 3e-3;
ntaps_DDLMS = 31;

% ==========================Tx=========================%
data = randi([0,M-1], lengthSequence,2);
sig_ori = pskmod(data,M,pi/M);
sig_x = sig_ori(:,1);
sig_y = sig_ori(:,2);
% -----------------------------------------------------%
figure;
plot(real(sig_x), imag(sig_x),'.','MarkerSize',15);
hold on;
plot(real(sig_y), imag(sig_y),'.','MarkerSize',15);
title('QPSK');
legend('X-Pol', 'Y-Pol');
% ---------------------upsamping-----------------------%
sig_resample = upsample(sig_ori, sps);
% ---------------------pulse shaping-------------------%
filterSymbolLength = 32; % 2*10*16+1
filterCoeffs2 = rcosdesign(roll_off,filterSymbolLength, sps,'sqrt').';

sig_ps_x = conv(sig_resample(:,1), filterCoeffs2, 'same');
sig_ps_x = (sig_ps_x-mean(sig_ps_x))./sqrt(mean(abs(sig_ps_x).^2));
sig_ps_y = conv(sig_resample(:,2), filterCoeffs2, 'same');
sig_ps_y = (sig_ps_y-mean(sig_ps_y))./sqrt(mean(abs(sig_ps_y).^2));
sig_ps = [sig_ps_x, sig_ps_y];

% -----------------------------------------------------%
figure;
plot(real(sig_ps(:,1)), imag(sig_ps(:,1)),'.');
hold on;
plot(real(sig_ps(:,2)), imag(sig_ps(:,2)),'.');
title('rrc');
legend('X-Pol', 'Y-Pol');

%% RSOP 
t_spacing =5e-11;
T = length(sig_ps);
t = 0:t_spacing:T*t_spacing;
N =10;
w_alpha = 53e3; %krad/s---rad/s
w_phi = 51e3;
w_kappa =163e3;
alpha= zeros(N,1);
phi = zeros(N,1);
kappa = zeros(N,1);
U = cell(1,T);

kappa0 = pi+randn(N,1)*2*pi;
alpha0 = pi+randn(N,1)*2*pi;
phi0 = pi+randn(N,1)*2*pi;

for a = 1:T
    alpha = w_alpha.*t(a)+alpha0;
    phi = w_phi.*t(a)+phi0;
    kappa = w_kappa.*t(a)+kappa0;
    oo11 = cos(kappa).*exp(1j*(alpha));
    oo12 = -sin(kappa).*exp(-1j*(phi));
    oo21 = -conj(oo12);
    oo22 = conj(oo11);
    y=1;
    for b = 1:N
        y = y*[oo11(b) oo12(b);oo21(b) oo22(b)];
    end
    U{a} = y;
Eoutx(a) = U{a}(1,1).*sig_ps_x(a)+U{a}(1,2).*sig_ps_y(a);
Eouty(a) = U{a}(2,1).*sig_ps_x(a)+U{a}(2,2).*sig_ps_y(a);
end
Eoutx = Eoutx';
Eouty = Eouty';
Eout = [Eoutx Eouty];
figure;
plot(real(Eout(:,1)), imag(Eout(:,1)),'.');                                         
hold on;
plot(real(Eout(:,2)), imag(Eout(:,2)),'.');
title('经过rsop后的星座图');
legend('X-Pol', 'Y-Pol');
axis equal;
axis([-2.1,2.1,-2.1,2.1]); 



% %================ 匹配滤波 =======================%
sig_hps_x = conv(Eoutx, filterCoeffs2, 'same');
sig_hps_x = (sig_hps_x-mean(sig_hps_x))./sqrt(mean(abs(sig_hps_x).^2));
sig_hps_y = conv(Eouty, filterCoeffs2, 'same');
sig_hps_y = (sig_hps_y-mean(sig_hps_y))./sqrt(mean(abs(sig_hps_y).^2));
sig_hps = [sig_hps_x, sig_hps_y];
%--------------------------------------------------------%
figure;
plot(real(sig_hps_x), imag(sig_hps_x),'.');
hold on;
plot(real(sig_hps_y), imag(sig_hps_y),'.');
title('RRC');
legend('X-Pol', 'Y-Pol');
% 




%% CMA
sig_hps_x1 = sig_hps_x;
sig_hps_y1 = sig_hps_y;
sig_hps_1 = [sig_hps_x1, sig_hps_y1];
sig_CMA_in = sig_hps_1;
[sig_out,~,~] = BFE_CMA_T2(sig_CMA_in,mu_CMA,ntaps);
sig_cma_xo = sig_out(:,1);
sig_cma_yo = sig_out(:,2); 
% --------------------------------------------------------%
figure;
plot(real(sig_cma_xo),imag(sig_cma_xo),'.');
hold on;
plot(real(sig_cma_yo),imag(sig_cma_yo),'.');
title('CMA');
legend('X-Pol', 'Y-Pol');
% 
% 
% 
%% ================ DD-LMS ===========================%

% P = pskmod(0:M-1,M,pi/M).';
% sig_LMS_xi = (sig_cma_xo-mean(sig_cma_xo))./sqrt(mean(abs(sig_cma_xo).^2));
% sig_LMS_yi = (sig_cma_yo-mean(sig_cma_yo))./sqrt(mean(abs(sig_cma_yo).^2));
% [sig_LMS_xo, sig_LMS_yo] = DDLMS_v0(sig_LMS_xi,sig_LMS_yi,mu_DDLMS,ntaps_DDLMS,1,P);
% sig_LMS_x1 = [zeros(33,1);sig_LMS_xo];
% sig_LMS_y1 = [zeros(33,1);sig_LMS_yo];
% figure;
% plot(real(sig_LMS_xo(end-10000:end)), imag(sig_LMS_xo(end-10000:end)),'.');
% hold on;
% plot(real(sig_LMS_yo(end-10000:end)), imag(sig_LMS_yo(end-10000:end)),'.');
% title('DD-LMS');
% legend('X-Pol', 'Y-Pol'); 
% rxSig = [sig_LMS_x1 sig_LMS_y1];
% dataOut = pskdemod(rxSig,4,0);
% [ber,~,~] = CalcBER(dataOut(:,1),data(:,1));

% MM = 4;
% NN = 2;
% % row = size(sig_out,1);
% % coloumn = size(sig_out,2);
% Eo=zeros(size(sig_out,1),size(sig_out,2));
% nModes = size(sig_out, 2);%column of E
% zeroPad = zeros(NN, nModes);
% x = [zeroPad; sig_out; zeroPad]; 
% D=zeros(size(x));
% theta=zeros(size(sig_out));
% for model = 1:nModes
% for index = 2:length(x)-1
%     d=x(index,model)*conj(x(index-1,model));
%     D(index,model)=d.^MM;
% 
%     if index >= 2*NN
%         % sum row
%         sumD = mean(D(:,model));
%         % get the best theata
%         theta(index-2 * NN+1,model) = angle(sumD)/4;
%     end
% end
% Eo(:,model)=sig_out(:,model).* exp(-1i * theta(:,model));                                                                                                                                           
% end
% figure;
% plot(real(Eo(:,1)),imag(Eo(:,1)),'.');
% hold on;
% plot(real(Eo(:,2)),imag(Eo(:,2)),'.');
% title('VVF');
% legend('X-Pol', 'Y-Pol');

%% phase compensation

V_p2xr = real(sig_cma_xo).^2-imag(sig_cma_xo).^2;
V_p2xi = 2*real(sig_cma_xo).*imag(sig_cma_xo);
V_p2x = V_p2xr+1j*V_p2xi;

V_p4xr = V_p2xr.^2-V_p2xi.^2;
V_p4xi = 2*V_p2xr.*V_p2xi;
V_p4x = V_p4xr+1j*V_p4xi;

V_p2yr = real(sig_cma_yo).^2-imag(sig_cma_yo).^2;
V_p2yi = 2*real(sig_cma_yo).*imag(sig_cma_yo);
V_p2y = V_p2yr+1j*V_p2yi;

V_p4yr = V_p2yr.^2-V_p2yi.^2;
V_p4yi = 2*V_p2yr.*V_p2yi;
V_p4y = V_p4yr+1j*V_p4yi;

Np = 158;
for i = Np:2:length(sig_cma_xo)

in_xr = V_p4xr(i:-1:(i-Np+1));
sum_in_xr = sum(in_xr);
in_xi = V_p4xi(i:-1:(i-Np+1));
sum_in_xi = sum(in_xi);
                                                                                                     
in_yr = V_p4yr(i:-1:(i-Np+1));
sum_in_yr = sum(in_yr);
in_yi = V_p4yi(i:-1:(i-Np+1));
sum_in_yi = sum(in_yi);
PEx14 = atan(sum_in_xi/sum_in_xr);
PEy14 = atan(sum_in_yi/sum_in_yr);

if (sum_in_xr <0 && sum_in_xi <0)
PEx4 = PEx14-pi;
end
if (sum_in_xr <0 && sum_in_xi >0)
PEx4 = PEx14+pi;
end
if (sum_in_xr >0 && sum_in_xi >0)
PEx4 = pi/2;
end
if (sum_in_xr >0 && sum_in_xi <0)
PEx4 = -pi/2;
end

if (sum_in_yr <0 && sum_in_yi <0)
PEy4 = PEy14-pi;
end
if (sum_in_yr <0 && sum_in_yi >0)
PEy4 = PEy14+pi;
end
if (sum_in_yr >0 && sum_in_yi >0)
PEy4 = pi/2;
end
if (sum_in_yr >0 && sum_in_yi <0)
PEy4 = -pi/2;
end

PEx = PEx4/4;
PEy = PEy4/4;   
tt = ceil(i/2);
out_xr(i-Np+1) = in_xr(tt).*cos(PEx)+in_xi(tt).*sin(PEx);
out_xi(i-Np+1) = in_xi(tt).*cos(PEx)-in_xr(tt).*sin(PEx);
out_yr(i-Np+1) = in_yr(tt).*cos(PEy)+in_yi(tt).*sin(PEy);
out_yi(i-Np+1) = in_yi(tt).*cos(PEy)-in_yr(tt).*sin(PEy);
end

out_x = out_xr+1j*out_xi;
out_y = out_yr+1j*out_yi;


figure;
plot(out_xr,out_xi,'.');
hold on;
plot(out_yr,out_yi,'.');
title('phase');
legend('X-Pol', 'Y-Pol');



