% 符号同步
function tx = symbolSync(rx, tx, SpS, mode)
    % Symbol synchronizer.
    if nargin < 4
        mode = 'amp';
    end
    %tx rx should be N ×（1,2）
    nModes = size(rx, 2);
    %downsample
    rx = rx(1:SpS:end, :);

    % Calculate time delay
    delay = zeros(1, nModes);

    corrMatrix = zeros(nModes, nModes);

    if strcmp(mode, 'amp')
        for n = 1:nModes
            for m = 1:nModes
                b=xcorr(abs(tx(:, m)-mean(tx(:, m))), abs(rx(:, n)-mean(rx(:, n))));
                corrMatrix(m, n) = max(abs(b));
            end
        end
        %1*1 为tx model1 对应 rx model1
        %2*1 为tx model1 对应 rx model2
        %1*2 为tx model2 对应 rx model1
        %2*2 为tx model2 对应 rx model2
        %找到对于每一个 rx 模式，与之相关性最高的 tx 模式
        % 对列 求最大值， 找到相对应的 rx model
        [~, swap] = max(corrMatrix,[], 1);
        % 进行重排
        tx = tx(:, swap);

        for k = 1:nModes
            %用tx相关性最大的模式找两个模式之间的延迟，实现同步
            delay(k) = finddelay(abs(tx(:, k)-mean(tx(:, k))), abs(rx(:, k)-mean(rx(:, k))));
        end
    elseif strcmp(mode, 'real')
        for n = 1:nModes
            for m = 1:nModes
                corrMatrix(m, n) = max(abs(xcorr(real(tx(:, m)), real(rx(:, n)))));
            end
        end

        %1*1 为tx model1 对应 rx model1
        %2*1 为tx model1 对应 rx model2
        %1*2 为tx model2 对应 rx model1
        %2*2 为tx model2 对应 rx model2
        %找到对于每一个 rx 模式，与之相关性最高的 tx 模式
        [~, swap] = max(corrMatrix,[],1);
        % 进行重排
        tx = tx(:, swap);

        for k = 1:nModes
            delay(k) = finddelay(real(tx(:, k)-mean(tx(:, k))), real(rx(:, k))-mean(rx(:, k)));
        end
    end

    % Compensate time delay
    % 为负数，向右移动
    for k = 1:nModes
        % delay is the index , so we should cut the one to remove
        delay(k)=delay(k)-1;
        tx(:, k) = circshift(tx(:, k), -round(delay(k)));
    end


    function d = finddelay(x, y)
        % Find delay between x and y.
    %计算峰值位置与输入信号序列长度的差，并加1，以获得信号同步的偏移量。
        [~, d] = max(xcorr(x-mean(x), y-mean(y)));
        d = d - length(x) + 1;
    end
end

