function p = gaussian(x, mu, sigma)
    p = exp(-((x - mu).^2) / (2 * sigma.^2)) / sqrt(2 * pi * sigma.^2);
end
