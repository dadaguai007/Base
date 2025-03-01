% intal parameter
SpS = 10;
Rs  = 10e6;          
Ts  = 1/Rs ;         
Fs  = SpS*Rs;    
Ta  = 1/Fs;
%Laser 
Plo_dBm  = 10;   
f_lo = 0 ;  
phi_lo  = 0 ;    
lw = 1e6;
Plo =10.^(Plo_dBm/10)*1e-3 ;
% Test
N=80000;
t_lo = (0:N-1) * Ta;
pn_lo  = phaseNoise(lw,N,Ta);
sigLO = sqrt(Plo)*exp(1i*(2*pi*f_lo*t_lo + phi_lo+pn_lo));
plot_spectrum(sigLO,Fs);

%
f= Fs * (-0.5:1/N:0.5-1/N);
S = lw./(2*pi*(f.^2+lw.^2));
figure;
plot(f,S)