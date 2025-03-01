clc;clear;
% received optical power
Pin_dBm = -1;   % Average optical power in dBm
Rd  = 0.5;      % Responsivity in A/W
Id  = 100e-9;   % Dark current in nA
B   = 10e9;     % Receiver bandwidth in Hz
q  = 1.60217663e-19; % Elementary charge in coulombs
Pin = 10^(Pin_dBm/10)*1e-3; % Average optical power in W
Ip  = Rd*Pin;

% Shot noise (white Gaussian noise)
Nsamples = 4300000;

sigma_s = 2*q*(Ip + Id)*B;  % Variance
mu    = 0;                % Mean
sigma = sqrt(sigma_s);

%gaussian
Is = normrnd(mu, sigma, [1,Nsamples]);

figure;hold on;
plot(Is)
% gaussian
x = -6 * sigma:sigma/10:6 * sigma;
fdp = gaussian(x, mu, sigma);
figure;hold on;grid on;
histogram(Is, 51, 'DisplayName', 'histograma');
plot(x,fdp,'Color', 'red','LineWidth', 2)
box on;
set(gca, 'LineWidth', 2);
set(gca, 'FontName', 'Arial', 'FontSize', 16);

% Tempurture noise
clc;clear;
Tc = 25;           % Temperature in Celsius
B  = 10e9;         % Receiver bandwidth
RL = 50;           % RL in Ohms
T = Tc + 273.15;   % Temperature in Kelvin
kB  = 1.380649e-23; % Boltzmann constant (in J/K)

% Thermal noise (white Gaussian noise)
Nsamples = 100000;
sigma2_T = 4 * kB * T * B / RL; % Variance
mu = 0;                % Mean
sigma = sqrt(sigma2_T);
It = normrnd(mu, sigma, 1, Nsamples);

% Plot the first 1000 samples
figure;hold on;
plot(It, 'LineWidth', 0.8);
xlabel('Sample');
grid on;
