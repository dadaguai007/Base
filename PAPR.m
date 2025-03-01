function [PAPR_dB, AvgP_dB, PeakP_dB] = PAPR(x)
% PAPR_dB  : PAPR[dB]
% AvgP_dB  : Average power[dB]
% PeakP_dB : Maximum power[dB]

Nx=length(x); 
xI=real(x); 
xQ=imag(x);

Power = xI.*xI + xQ.*xQ;

PeakP = max(Power); 
PeakP_dB = 10*log10(PeakP);

AvgP = sum(Power)/Nx; 
AvgP_dB = 10*log10(AvgP);

PAPR_dB = 10*log10(PeakP/AvgP);