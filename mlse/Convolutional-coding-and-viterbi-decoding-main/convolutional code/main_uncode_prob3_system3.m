% ======================================================================= %
% SSY125 Project
% ======================================================================= %
clc
clear

% ======================================================================= %
% Simulation Options
% ======================================================================= %
N = 1e5;  % simulate N bits each transmission (one block)
maxNumErrs = 100; % get at least 100 bit errors (more is better)
maxNum = 1e6; % OR stop if maxNum bits have been simulated
EbN0 = -1:8; % power efficiency range

% ======================================================================= %
% Other Options
% ======================================================================= %
% ...

% ======================================================================= %
% Simulation Chain
% ======================================================================= %
 
 BER = zeros(1, length(EbN0)); % pre-allocate a vector for BER results 
 Pb = zeros(1, length(EbN0));
for i = 1:length(EbN0) % use parfor ('help parfor') to parallelize  
  totErr = 0;  % Number of errors observed
  num = 0; % Number of bits processed

  while((totErr < maxNumErrs) && (num < maxNum))
  % ===================================================================== %
  % Begin processing one block of information
  % ===================================================================== %
  % [SRC] generate N information bits 
  
  tx_bit = randi([0 1],1,N);

  % [ENC] convolutional encoder

  % [MOD] symbol mapper
  
  tx_sym =modulation(tx_bit,'AMPM');


  % scatterplot(tx_sym);

  % [CHA] add Gaussian noise

  EbN0_value = 10^(EbN0(i)/10);
  R = length(tx_bit)/length(tx_sym);
  SNR = R * EbN0_value; 

  sigma2 = 1/(SNR*2);
  noise = sqrt(sigma2)*(randn(1,length(tx_sym)) + 1i*randn(1,length(tx_sym)));
  rx_data = tx_sym +noise;
  
  % scatterplot(rx_data);


  % demodulation

  bit_out = demodulation(rx_data,'AMPM');

  % ===================================================================== %
  % End processing one block of information
  % ===================================================================== %
  
  BitErrs = length(find(bit_out ~= tx_bit)); % count the bit errors and evaluate the bit error rate
  totErr = totErr + BitErrs;
  num = num + N; 

  disp(['+++ ' num2str(totErr) '/' num2str(maxNumErrs) ' errors. '...
      num2str(num) '/' num2str(maxNum) ' bits. Projected error rate = '...
      num2str(totErr/num, '%10.1e') '. +++']);
  end 
  BER(i) = totErr/num; 

  % theoritical BER for QPSK gray map
  Pb(i) = qfunc(sqrt(2*EbN0_value));

end

% plots
figure(1)  
a= semilogy(EbN0, BER,'--^','DisplayName','uncoded transmission(system 3)');
a.Color = 	'#7E2F8E';
a.LineWidth =1;
ylim([1e-4 1]);
hold on
grid on
xlabel('EbN0(dB)');
ylabel('BER');
legend

% ======================================================================= %
% End
% ======================================================================= %