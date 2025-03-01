function [Et, Pt]=eam(Ein,V,param)

% filt = design_filter('bessel', 5, 0.7*ModFormat.Rs/(sim.fs/2));
% Mod.H = filt.H(f/fs);
if isfield(param, 'f')
    f = param.f;
end
if isfield(param, 'alpha')
    alpha = param.alpha;
end
if isfield(param, 'H')
    H = param.H;
end

% ifftshift is sure the zero is under the center
Pt = abs(Ein).^2./mean(abs(Ein).^2).*real(ifft(fft(V).*ifftshift(H(f))));

% make sure  all pulse is not negative number
Pt(Pt <= 0) = 0;

Et = sqrt(Pt);

% Add transient chirp
% Calculate electric field including chirp
% the eam do not consider adiabatica chirp
dphi = alpha/2*log(Pt + realmin);         % Phase variation due to chirp (only transient chirp is considered)
%note: realmin,smallest positive normalized number.
% to avoid log(0)
Et = Et.*exp(1j*dphi);



end