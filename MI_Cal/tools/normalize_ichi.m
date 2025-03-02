function out = normalize_ichi(in, type)
    if isreal(in)
        switch type
            case 'max'
                out = in ./ max(abs(in));
            case 'gauss'
                avg = mean(in);
                sigma2 = var(in);
                out = (in-avg) ./ sqrt(sigma2);
            case 'meanpwr'
                out = in*length(in) / (sum(sqrt(in.^2)));
        end
    else
        in_real = real(in);
        in_imag = imag(in);
        switch type
            case 'max'
                in_real = in_real ./ max(abs(in_real));
                in_imag = in_imag ./ max(abs(in_imag));
                out = in_real + 1j*in_imag;
            case 'gauss'
                avg_real = mean(in_real);
                avg_imag = mean(in_imag);
                sigma2_real = var(in_real);
                sigma2_imag = var(in_imag);
                in_real = (in_real-avg_real) ./ sqrt(sigma2_real);
                in_imag = (in_imag-avg_imag) ./ sqrt(sigma2_imag);
                out = in_real + 1j*in_imag;
            case 'meanpwr'
                out = in*length(in) / (sum(sqrt(in.*conj(in))));
        end
    end
end

