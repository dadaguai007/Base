% 脉冲成型滤波器

filterSymbolLength = 202;
switch lower(obj.pulseShape)
    case 'rz33'
        obj.filterCoeffs = cos((pi/2)*sin(pi*(0:(1/obj.samplesPerSymbol):1)-pi/2));
    case 'rz50'
        obj.filterCoeffs = cos((pi/4)*sin(2*pi*(0:(1/obj.samplesPerSymbol):1)-pi/2)-pi/4);
    case 'rz66'
        obj.filterCoeffs = cos((pi/2)*sin(pi*(0:(1/obj.samplesPerSymbol):1))-pi/2);
    case 'nrz'
        obj.filterCoeffs = ones(1,obj.samplesPerSymbol);
    case 'rc'
        if ~isfield(param, 'rollOff')
            robolog('Please set a roll-off factor', 'ERR')
        end
        filterFreqs = linspace(-(obj.filterSymbolLength/2), (obj.filterSymbolLength/2), obj.samplesPerSymbol*obj.filterSymbolLength+1);
        obj.filterCoeffs = sinc(filterFreqs).*cos(pi*obj.rollOff*filterFreqs)./(1-4*obj.rollOff^2*filterFreqs.^2);
        obj.filterCoeffs(abs(filterFreqs) == 1/(2*obj.rollOff)) = (pi/4)*sinc(1/(2*obj.rollOff));
    case 'rrc'
        if ~isfield(param, 'rollOff')
            robolog('Please set a roll-off factor', 'ERR')
        end
        filterFreqs = linspace(-(obj.filterSymbolLength/2), (obj.filterSymbolLength/2), obj.samplesPerSymbol*obj.filterSymbolLength+1);
        obj.filterCoeffs = (sin(pi*filterFreqs*(1-obj.rollOff)) + 4*obj.rollOff*filterFreqs.*cos(pi*filterFreqs*(1+obj.rollOff)))./(pi*filterFreqs.*(1-(4*obj.rollOff*filterFreqs).^2));
        obj.filterCoeffs(abs(filterFreqs) == 1/(4*obj.rollOff)) = (obj.rollOff/sqrt(2))*( (1+2/pi)*sin(pi/(4*obj.rollOff)) + (1-2/pi)*cos(pi/(4*obj.rollOff)));
        obj.filterCoeffs(filterFreqs == 0) = 1 - obj.rollOff + 4*obj.rollOff/pi;
    otherwise
        robolog('Please, define one of the allowed pulse shapes ("NRZ", "RZ33", "RZ50", "RZ67", "RC" or "RRC".', 'ERR')
end


%执行上采样、滤波和截断操作
function out = ps(in, samplesPerSymbol, filterCoeffs)
% Upsampling
in = upsample(in, samplesPerSymbol);
% Shaping 成型滤波
out = conv([in(:) ; in(1:10,1)], filterCoeffs, 'same');
% truncating
%截断
out = out(1:length(in));
end