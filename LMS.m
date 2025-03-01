function [out, h, squaredError] = LMS(y, x, Ntaps, mu)
    % Apply the LMS (Least Mean Squares) algorithm for equalization.

    % y - input signal
    % x - reference signal

    % Initialize the equalizer filter coefficients
    h = zeros(1, Ntaps);
    L = floor(length(h) / 2); % Decision delay

    % Apply the LMS algorithm
    squaredError = zeros(size(y));
    out = zeros(size(y));
    ind = -L:L;
    
    % 在两侧各自填充L个0
    y = padarray(y, [0,L], 0, 'both');
    

    % Iterate through each sample of the signal
    for i = (L+1):(length(y)-L)
        yVec = y(i+ind);
        yVec=yVec(end:-1:1);

        % Generate the estimated signal using the equalizer filter
        xhat = yVec * h';

        % Compute the error between the estimated signal and the reference signal
        error = x(i-L) - xhat;

        % Update the filter coefficients using the LMS update rule
        h = h + mu * yVec * error;

        squaredError(i-L) = error.^2;
        out(i-L) = xhat;
    end
end
