function u1 = ssprop(u0, fs, dz, nz, alpha, beta, gamma)
    dt = 1/fs;
    nt = length(u0);
    w = 2*pi*[(0:nt/2-1),(-nt/2:-1)]'/(dt*nt);

    linearoperator = -alpha/2;
    for ii = 0:length(beta)-1
        linearoperator = linearoperator - 1j*beta(ii+1)*(w).^ii/factorial(ii);
    end
    % linearoperator = linearoperator - 1j*beta(3)*(w).^2/factorial(2) - 1j*beta(4)*(w).^3/factorial(3);
    % linearoperator = linearoperator - 1j*beta(3)*(w).^2/factorial(2);
    halfstep = exp(linearoperator*dz/2);
    ufft = fft(u0);
    for i = 1:nz
        uhalf = ifft(halfstep.*ufft);
        nonlinearoperator = -1j * gamma * abs(uhalf).^2;
        uhalfnl = uhalf .* exp(nonlinearoperator*dz);
        ufft = halfstep .* fft(uhalfnl);
    end
    u1 = ifft(ufft);
end

