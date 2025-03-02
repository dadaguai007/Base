function out = splineResampler(in, outputSamplingRate, inputSamplingRate)
    k = outputSamplingRate / inputSamplingRate;
    inputSamplingTime = 0:1/inputSamplingRate:((length(in)-1)/inputSamplingRate);
    outputSamplingTime = 0:1/outputSamplingRate:inputSamplingTime(end)+(k-1)/outputSamplingRate;

    out = interp1(inputSamplingTime, in.', outputSamplingTime, 'spline', 0);
    out = out.';
end