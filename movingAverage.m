function y = movingAverage(x, N)
    % Calculate the sliding window moving average of a 2D NumPy array along each column.

    % Input:
    % x: Input 2D array with shape (M, N), where M is the number of samples and N is the number of columns.
    % N: Size of the sliding window.

    % Output:
    % y: 2D array containing the sliding window moving averages along each column.

    % Notes:
    % The function pads the signal with zeros at both ends to compensate for the lag between the output
    % of the moving average and the original signal.

if size(x,1)<size(x,2)
x=x.';
end

    [~, nCol] = size(x);
    y = zeros(size(x));

    startInd = floor(N/2);

    if mod(N, 2) == 0
        endInd = -floor(N/2) + 1;
    else
        endInd = -floor(N/2);
    end
    % endInd 设置为负数，直接进行相减即可
    for indCol = 1:nCol
        % 在数据前部和后部都补零，来执行一维卷积，实现
        % Pad the signal with zeros at both ends
        padded_x = [zeros(startInd, 1); x(:, indCol); zeros(-endInd, 1)];

        % Calculate moving average using convolution
        h = ones(N, 1) / N;
        out = conv(padded_x, h, 'same');
        y(:, indCol) = out(startInd+1:end+endInd);
    end
end
