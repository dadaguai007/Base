clc;clear;
% Rz pattern
%First ， Generate an NRZ pattern
%Use MZM handle NRZ with some specific parameters

% time vector RZ use NRZ
% Define RZ pulse type (33%, 50%, 67%)
RZ = 33;

% Parameters of the MZM
Vpi = 2;
Ai = 1;
% Parameters of signal
Rs=20;
sps=10;
% NRZ


% Parameters for the formatter for each RZ pulse type
if RZ == 33
    Vb = 0;
    % Parameters for the sinusoidal signal
    fs = Rs / 2;
    Vs = Vpi;
    phi_s = pi/2;
elseif RZ == 50
    Vb = Vpi/2;
    % Parameters for the sinusoidal signal
    fs = Rs;
    Vs = Vpi/2;
    phi_s = pi;
elseif RZ == 67
    Vb = Vpi;
    % Parameters for the sinusoidal signal
    fs = Rs / 2;
    Vs = Vpi;
    phi_s = 0;
end

% Generate a sinusoidal signal
t = (0:(length(sigTxo_NRZ)-1)) * (1/(Rs*sps));
% RF_modularRZ
senoideRF = Vs * cos(2 * pi * fs * t + phi_s);

% MZM used as a pulse formatter (pulse carver)
sigTxo_Rz = mzm(sigTxo_NRZ, senoideRF, Vb, Vpi);

Nsamples = 10000;
