clc;
clear;
close all;

% Estimate model order using decay of reflection coefficients.
rng default;
y=filter(1,[1 -0.75 0.5],0.2*randn(1024,1));
% Create AR(2) process
[ar_coeffs,NoiseVariance,reflect_coeffs] = aryule(y,10);

% Fit AR(10) model
stem(reflect_coeffs);
axis([-0.05 10.5 -1 1]);
title('Reflection Coefficients by Lag');
xlabel('Lag');ylabel('Reflection Coefficent');


%%
clc;
clear;
close all;

randn('state',0); 
noise = randn(50000,1);  % Normalized white Gaussian noise
x = filter(1,[1 1/2 1/3 1/4],noise);
x = x(45904:50000);
a = lpc(x,3);
est_x = filter([0 -a(2:end)],1,x);  % Estimated signal
e = x - est_x;                      % Prediction error
[acs,lags] = xcorr(e,'coeff');      % ACS of prediction error

%   Compare the predicted signal to the original signal
figure;
plot(1:97,x(4001:4097),1:97,est_x(4001:4097),'--');
title('Original Signal vs. LPC Estimate');
xlabel('Sample Number'); ylabel('Amplitude'); grid on;
legend('Original Signal','LPC Estimate')

%   Look at the autocorrelation of the prediction error.
figure;
plot(lags,acs); 
title('Autocorrelation of the Prediction Error');
xlabel('Lags'); ylabel('Normalized Value'); grid on;

%%
clc;
clear;
close all;

%  Estimate input noise variance for AR(4) model.
A = [1 -2.7607 3.8106 -2.6535 0.9238]; 
% Generate noise standard deviations
rng default;
noise_stdz = rand(1,50)+0.5;

% Generate column vectors that have corresponding standard deviation
x = bsxfun(@times,noise_stdz,randn(1024,50));

% filter each column using the AR model.
y = filter(1,A,x);

% Compute the estimated coefficients and deviations for each column
[ar_coeffs,NoiseVariance]=arburg(y,4);

% Display the mean value of each estimated polynomial coefficient
estimatedA = mean(ar_coeffs);

% Compare actual vs. estimated variances
plot(noise_stdz.^2,NoiseVariance,'k*');
xlabel('Input Noise Variance');
ylabel('Estimated Noise Variance');

