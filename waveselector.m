% IMDD WDM wavechannel selecter
function  dataout=waveselector(Es,param)
%type is 'match' or 'lpf'
% output is the wavelength of the index of the wave channel
if iscolumn(Es)
    Es=Es.';
end
Fs=param.Fs;
freqGrid=param.freqGrid;
chIndex=param.chIndex;
hsqrt=param.hsqrt;
type=param.type;
Bandwidth=param.Bandwidth;
f_lo=(freqGrid(chIndex)/Fs);
t= (0:length(Es)-1);


Elo = exp(1i*2*pi*f_lo*t);

% the row vector
s= Es.*Elo;

% filter
if strcmp(type,'match')
    % match filter
    dataout=conv(s,hsqrt,'same');
elseif strcmp(type,'lpf')
    % ideal filter LPF
    dataout=LPF(s,Fs,Bandwidth);
end

end