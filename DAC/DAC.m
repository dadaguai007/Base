function [wout, snr_dB] = DAC(win, param)
    [wout, snr_dB] = quantizer(win, param.quantizer.bits, param.quantizer.ENoB, param.quantizer.type, param.snr_dB);
    wout = resample(wout, param.dac.targetFs, param.dac.Fs);
    wout = efilter(wout, param.efilter.type, param.efilter.order, param.efilter.bwl, param.dac.targetFs, param.efilter.sps, param.efilter.outputVoltage);
end