function sigRxo=Attenuator(sigTxo,Att_dB)
% 衰减，以dB为单位值
sigRxo=sigTxo*10.^(-Att_dB/20);

end