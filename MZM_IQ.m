function [Eout] = MZM_IQ(Ein,VRF_i,Vbias_i,VRF_q,Vbias_q,Vbias_p,Vpi_i,Vpi_q,Vpi_p,gamma_i,gamma_q,gamma_p)
if nargin < 10
    gamma_i = 1;
    gamma_q = 1;
    gamma_p = 1;
end
% MZM_I
Ein_i = Ein / 2;
MZM_i = MZM_SingleDrive(Ein_i,VRF_i,Vbias_i,Vpi_i,gamma_i);
% MZM_Q
Ein_q = Ein / 2;
MZM_q = MZM_SingleDrive(Ein_q,VRF_q,Vbias_q,Vpi_q,gamma_q);
% phase
phase = Vbias_p/Vpi_p*pi;

Eout = MZM_i + gamma_p*MZM_q.*exp(1i.*phase);

end

