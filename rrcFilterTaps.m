function coeffs = rrcFilterTaps(t, alpha, Ts)
% Generate Root-Raised Cosine (RRC) filter coefficients.
coeffs = zeros(length(t),1);

for i = 1:length(t)
    t_i = t(i);
    t_abs = abs(t_i);

    if t_i == 0
        coeffs(i) = (1 / Ts) * (1 + alpha * (4 / pi - 1));
    elseif t_abs == Ts / (4 * alpha)
        term1 = (1 + 2 / pi) * sin(pi / (4 * alpha));
        term2 = (1 - 2 / pi) * cos(pi / (4 * alpha));
        coeffs(i) = (alpha / (Ts * sqrt(2))) * (term1 + term2);
    else
        t1 = pi * t_i / Ts;
        t2 = 4 * alpha * t_i / Ts;
        coeffs(i) = (1 / Ts) * (sin(t1 * (1 - alpha)) + (4 * alpha * t_i / Ts) * cos(t1 * (1 + alpha))) / (pi * t_i * (1 - t2^2));
    end
end
end

%note:
% RRC滤波器系数的计算分为三种方式，是为了处理滤波器的不同时间域区间，以满足滤波器的设计要求和特性。每个计算方式对应不同的时间点和滤波器特性。
% 1. **t_i == 0：** 这是计算RRC滤波器系数的中心点。在这一点上，系数的计算基于RRC滤波器的中心值。这个中心值通常用于补偿信号的延迟，以确保信号保持时间对齐。
% 2. **t_abs == Ts / (4 * alpha)：** 这是计算系数的过渡点之一，也被称为RRC滤波器的过渡点。在这一点上，系数的计算基于RRC滤波器的过渡特性，以确保信号的平滑过渡。
% 3. **其他情况：** 这是计算系数的一般情况，涵盖了滤波器的其余部分。在这些情况下，系数的计算基于RRC滤波器的一般性质，以满足滤波器的频率响应要求。
% 这三种计算方式涵盖了RRC滤波器在不同时间点的不同特性，以确保滤波器能够在时间域和频率域上具有所需的特性。这种方法能够实现滤波器的平滑过渡、时域对齐和频域响应，从而适应不同的通信系统需求。不同的RRC滤波器参数（如滚降因子 `alpha`）可以导致不同的过渡点和滤波器特性，因此需要不同的系数计算方式来满足这些不同的参数设置。