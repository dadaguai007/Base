function coeffs = rcFilterTaps(t, alpha, Ts)
    % Generate Raised Cosine (RC) filter coefficients.
    coeffs = zeros(length(t), 1);

    for i = 1:length(t)
        t_i = t(i);
        t_abs = abs(t_i);
        if t_abs == Ts/(2 * alpha)
            coeffs(i) = pi / (4 * Ts) * sinc(1 / (2 * alpha));
        else
            coeffs(i) = (1/Ts) * sinc(t_i/Ts) * cos(pi*alpha*t_i/Ts) / (1-4*alpha^2*t_i.^2/Ts.^2);
        end
    end
    
end
