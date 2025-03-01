
function [rmsEVM_symbol,rmsEVM_subcarrier,rmsEVM_martix]=EVM_Measure(data_remove_ssbi,qam_signal_mat)


    %EVM measure
    evm = comm.EVM(AveragingDimensions=1);
    rmsEVM_symbol = evm(data_remove_ssbi,qam_signal_mat);


    evm = comm.EVM(AveragingDimensions=2);
    rmsEVM_subcarrier = evm(data_remove_ssbi,qam_signal_mat);


    evm = comm.EVM(AveragingDimensions=[1 2]);
    rmsEVM_martix = evm(data_remove_ssbi,qam_signal_mat);
end