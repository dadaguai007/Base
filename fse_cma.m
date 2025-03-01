function [Y, W, MSE] = fse_cma(X, eq)
%% Fractionally spaced equalization (FSE) using constant modulus algorithm (CMA)
%% Equalization is done in time domain
% Inputs:
% - X : input samples in the two polarizations [2 x N] at rate ros x Rs
% - eq : equalizer parameters {eq.ros = oversampling ratio, eq.Ntrain = 
% number of traning symbols, eq.mu = adaptation rate}
% - verbose (optional, default=false): whether to plot equalizer frequency
% response and convergence
% Output:
% - Y : Equalized symbols at rate 1 x Rs
% - W : cell containing Wx, Wy, Wmix filters coefficients 
% - MSE : mean square error
ros = eq.ros;
mu = eq.mu;
Nsymb = floor(length(X)/ros);

if mod(eq.Ntaps, 2) == 0 % always uses odd number of taps
    eq.Ntaps = eq.Ntaps + 1;
end
Ntaps = eq.Ntaps;

% Number of filters
[Nsamp, Nfilters]  = rat(ros);

if strcmpi(eq.structure, '4 filters')
    %% Implementation using 4 filters
    % Initialize filters coefficients
    Wx = zeros(2*Ntaps, Nfilters); % used for x polarization
    Wy = zeros(2*Ntaps, Nfilters); % used for y polarization

    % Initialize filters as if there were no ISI and there was no mixing between the two pols
    Wx(ceil(Ntaps/2), :) = 1; 
    Wy(ceil(Ntaps/2)+Ntaps, :) = 1;
    
    ex = zeros(1, Nsymb);
    ey = zeros(1, Nsymb);
    yx_hat = zeros(1, Nsymb); % output in x pol stream
    yy_hat = zeros(1, Nsymb); % output in y pol stream
    filt = 1;
    window = (1:Ntaps) - ceil(Ntaps/2); % {-(Ntaps-1)/2, ..., 0, ..., (Ntaps-1)/2}.
    kstart = ceil(Ntaps/(2*Nsamp))*Nsamp+1; 
    k = kstart; % runs over samples
    for n = ((kstart-1)/ros+1):Nsymb-ceil(Ntaps/2) % n runs over symbols
        % k indexes samples and n indexes symbols
        z = [X(1, k + window), X(2, k + window)]; % input of filter [x pol, y pol]

        % Filter output
        yx_hat(n) = z*Wx(:, filt);
        yy_hat(n) = z*Wy(:, filt);

        % Calculate error
        ex(n) = 2 - abs(yx_hat(n))^2;
        ey(n) = 2 - abs(yy_hat(n))^2;

        % Update filters coefficients
        Wx(:, filt) = Wx(:, filt) + mu*ex(n)*yx_hat(n)*z';
        Wy(:, filt) = Wy(:, filt) + mu*ey(n)*yy_hat(n)*z';

        % Increment filter index
        filt = mod(filt, Nfilters) + 1;

        % ensures that the center of the window of samples remains 
        % close to the nth symbol
        if abs((k+1)/ros - n-1) > 0.5
            k = k + 2;
        else
            k = k + 1;
        end
    end
    W = {Wx, Wy};
    
elseif strcmpi(eq.structure, '2 filters')
    %% Implmentation using 2 filters
    % Number of operations is halved
    % Initialize filters coefficients
    Wx = zeros(Ntaps, Nfilters); % used for x polarization
    Wy = zeros(Ntaps, Nfilters); % used for y polarization    
    Wmix = zeros(2, Nfilters);  % crossover filters
    
    % Initialize filters as if there were no ISI and there was no mixing between the two pols
    Wx(ceil(Ntaps/2), :) = 1; 
    Wy(ceil(Ntaps/2), :) = 1;
        
    ex = zeros(1, Nsymb);
    ey = zeros(1, Nsymb);
    yx_hat = zeros(1, Nsymb); % output in x pol stream
    yy_hat = zeros(1, Nsymb); % output in y pol stream
    filt = 1;
    window = (1:Ntaps) - ceil(Ntaps/2); % {-(Ntaps-1)/2, ..., 0, ..., (Ntaps-1)/2}.
    kstart = ceil(Ntaps/(2*Nsamp))*Nsamp+1; 
    k = kstart; % runs over samples
    for n = ((kstart-1)/ros+1):Nsymb-ceil(Ntaps/2) % n runs over symbols
        % k indexes samples and n indexes symbols
        zx = X(1, k + window);
        zy = X(2, k + window); 

        % Filter output
        yx = zx*Wx(:, filt);
        yy = zy*Wy(:, filt);

        % Filter output
        yx_hat(n) = yx + yy*Wmix(1, filt);
        yy_hat(n) = yx*Wmix(2, filt) + yy;

        % Calculate error
        ex(n) = 2 - abs(yx_hat(n))^2;
        ey(n) = 2 - abs(yy_hat(n))^2;

        % Update filters coefficients
        Wx(:, filt) = Wx(:, filt) + mu*ex(n)*yx_hat(n)*zx';
        Wy(:, filt) = Wy(:, filt) + mu*ey(n)*yy_hat(n)*zy';

        Wmix(1, filt) = Wmix(1, filt) + mu*ex(n)*yx_hat(n)*conj(yy);    
        Wmix(2, filt) = Wmix(2, filt) + mu*ey(n)*yy_hat(n)*conj(yx);    

        % Increment filter index
        filt = mod(filt, Nfilters) + 1;

        % ensures that the center of the window of samples remains 
        % close to the nth symbol
        if abs((k+1)/ros - n-1) > 0.5
            k = k + 2;
        else
            k = k + 1;
        end
    end
    
    W = {Wx, Wy, Wmix};
end
    
% Build outputs
Y = [yx_hat; yy_hat];
MSE = [abs(ex).^2; abs(ey).^2];