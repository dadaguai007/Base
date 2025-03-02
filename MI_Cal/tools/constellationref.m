function constellation = constellationref(modulationFormat, M, varargin)
    
    switch lower(modulationFormat)
        case 'ask'
            if mod(M, 1) == 0
                constellation = 0:M-1;
            else
                slog('Constellation order wrong', 'ERR');
            end
        case 'ook'
            if M == 2
                constellation = [-1 1];
            else
                slog('Constellation order wrong', 'ERR');
            end
        case 'pam'
            if mod(M, 1) == 0
                constellation = 1-M:2:M-1;
            else
                slog('Constellation order wrong', 'ERR');
            end
        case 'qam'
            if mod(M, 2)
                slog('The modulation order is incorrect, and therefore, the correct constellation diagram cannot be generated.', 'ERR');
            end
            if M == 2
                constellation = constellationref('psk', 2);
            else
                tempM = nextsquarenumber2(M);
                if M == tempM
                    X = 1-sqrt(M):2:sqrt(M)-1;
                    [I, Q] = meshgrid(X);
                    constellation = I + 1j*Q;
                else
                    d = (tempM - M) / 4;
                    d = sqrt(d);
                    mask = ones(sqrt(tempM));
                    mask(1:d, 1:d) = 0;
                    mask(1:d, end-d+1:end) = 0;
                    mask(end-d+1:end, 1:d) = 0;
                    mask(end-d+1:end, end-d+1:end) = 0;
                    constellation = constellationref('qam', tempM);
                    constellation = constellation(logical(mask(:)));
                end
            end
        case 'psk'
            if ~mod(M, 1)
                constellation = exp(2*pi*1j*(0:M-1)/M);
            else
                slog('Constellation order wrong', 'ERR');
            end
        otherwise
            slog('Incorrect input modualtion format, this func only provides "QAM", "PAM", "OOK", "ASK", "PSK" modualtion format', 'ERR');
    end
    constellation = constellation(:).';
%     constellation = constellation / sqrt(pwr.meanpwr(constellation));
    if ~isempty(varargin)
        if length(varargin) == 1
            constellation = constellation.*exp(1j*varargin{1});
        else
            slog('Too many parameters are being input. The reserved parameter can only accept angles.', 'ERR');
        end
    end
end