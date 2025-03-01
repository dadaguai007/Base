function hammingDistance=hanmming(originalBits,decodedBits)

% 原始比特序列
% originalBits = [1 0 1 1 0 1];

% 解码后的比特序列（模拟）
% decodedBits = [1 1 0 1 0 1];

% 计算汉明距离，计算原始比特序列和解码后的比特序列中不同位置的比特数,两种实现方式
% hammingDistance = sum(originalBits ~= decodedBits);
hammingDistance = sum(xor(originalBits,decodedBits));
disp(['汉明距离：', num2str(hammingDistance)]);
