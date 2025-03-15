function [Eout] = MZM_DualDrive(Ein,VRF1,VRF2,Vbias1,Vbias2,Vpi,gamma)
if nargin < 7
    gamma = 1;
end
Eout = Ein/2.*(exp(1j*(VRF1+Vbias1)/Vpi*pi)+gamma*exp(1j*(VRF2+Vbias2)/Vpi*pi));
end

