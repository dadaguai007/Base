% 
function outparams = paramsref(objname, inparams)
    basenamestruct   = {'signal_interface'};
    txnamestruct     = {'ChannelCombiner', 'Dfilter', 'IntensityModulator', 'IQModulator', 'Mapper', 'PatternGenerator', 'PulseShaper_v1', 'PulseShaper_v2'};
    chnamestruct     = {'EDFA', 'Gain', 'NonlinearCh', 'setOSNR'};
    rxnamestruct     = {'AdaptiveEqualizer', 'BaseBandFilter', 'BERT', 'CDCompensator', 'DBPCompensator', 'DCBlock', 'DDPLL', 'IQCompensator'};
    devicenamestruct = {'BPD', 'Delayer', 'Efilter', 'Laser', 'OpticalHybrid', 'PBC', 'PBS', 'PD', 'Polarizer', 'Quantizer', 'Replicator', 'Resampler'};
    othernamestruct  = {'PhaseNoise', 'Kmeans'};
    namestruct = [basenamestruct, txnamestruct, chnamestruct, rxnamestruct, devicenamestruct, othernamestruct];
    if ismember(objname, namestruct)
        props = properties(objname);
        outparams = rmfield(inparams, setdiff(fieldnames(inparams), props));
    else
        slog('The information about the class is not available in the reference library. Please update the library if new classes are added.', 'ERR');
    end
end