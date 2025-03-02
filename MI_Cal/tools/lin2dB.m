% file lin2dB
% This is a function that convert linear units to dB.
% Only for power convertion.

function dB = lin2dB(lin, varargin)

type = defaultunit('W', varargin);

switch type
	case "W"
		fact = 1;
	case "mW"
		fact = 1e-3;
	case "uW"
		fact = 1e-6;
	otherwise

end

dB = 10*log10(lin*fact);