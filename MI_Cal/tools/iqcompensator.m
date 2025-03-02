function out = iqcompensator(in)
    xi = real(in);
    xq = imag(in);
    rho = mean(xi.*xq);
    
    Pi = pwr.meanpwr(xi);
    xq = xq - rho*xi/Pi;
    Pq = pwr.meanpwr(xq);
    out = pwr.normpwr(xi/sqrt(Pi) + 1j*xq/sqrt(Pq));
end