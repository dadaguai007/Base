function d = finddelay(x, y)
    % Find delay between x and y.

    % Ensure that the input signals are column vectors
    x = x(:);
    y = y(:);

    % Compute the cross-correlation between x and y
    correlation = xcorr(x, y);
    % Find the index of the maximum value in the cross-correlation
    [~, maxIndex] = max(correlation);

    % Calculate the delay in samples
    d = maxIndex - length(x) + 1;
end
