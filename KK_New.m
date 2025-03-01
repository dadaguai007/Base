function [compSig] = KK_New(rxSig,fOsc,fUp)
   % up-sampling
   upRxSig = resample(rxSig,fUp,fOsc); % 5 sps
   % square root
   hn = sqrt(upRxSig);
   % ln of abs
   lnSig = log(abs(hn));
   % fft
   fftSig = fft(lnSig);
   % i*sign(w)
   N = length(fftSig);
   freq = [0:N/2-1,-N/2:-1].';
   signFFTSig = -1i.*sign(freq).*fftSig;
   % ifft
   ifftSig = ifft(signFFTSig);
   % exp
   phi = exp(1i.*(ifftSig));
   % phase recovered signal
   compSig = hn.*phi;
end