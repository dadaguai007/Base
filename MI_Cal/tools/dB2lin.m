% file dB2lin.m
% This is a function that convert dB to linear units
% Only for power convertion.

function lin = dB2lin(dB, varargin)

type = defaultunit('dB', varargin);

switch type
	case "dB"
		fact = 0;
	case 'dBm'
		fact = -30;
	case 'dBu'
		fact = -60;
	otherwise
		error('Conversion type must be ''dB'', ''dBm'', ''dBu''.')
end

lin = 10.^((dB+fact)/10);