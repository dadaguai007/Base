function flatness = flatness_cal(data)
flatness = var(data/sqrt(bandpower(data)));
end

%这个 MATLAB 函数flatness_cal用于计算输入数据data的平坦度（flatness）。它通过一种特定的方式来量化数据在功率分布上的平坦程度。