function sigTotal=Combiner_Pow(sigTxo_1,sigTxo_2)

% CombinerPow
signal_Array=[sigTxo_1;sigTxo_2];
sigTotal=(1/sqrt(min(size(signal_Array))))*sum(signal_Array);

end