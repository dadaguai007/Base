function pxy=Calc_px_given_y(x, y, sigma2, chi)
py = 0;
for i =1:length(chi)
    z=calculate_py_given_x(chi(i), y, sigma2);
    py= py+0.25 * z;
end
Z=calculate_py_given_x(x, y, sigma2);
pxy = 0.25 * Z ./ py;
    function y=calculate_py_given_x(x, y, sigma2)
        y=(1/(sqrt(2*pi*sigma2))) * exp(-(y-x).^2/sigma2/2);
    end
end


%sigma2=mean(abs(X) ^ 2)/snr, chi=ref_symbols
