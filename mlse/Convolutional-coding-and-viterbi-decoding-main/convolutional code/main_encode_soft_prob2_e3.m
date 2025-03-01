% ======================================================================= %
% SSY125 Project
% ======================================================================= %
clc
clear

% ======================================================================= %
% Simulation Options
% ======================================================================= %
N = 1e7;  % simulate N bits each transmission (one block)
maxNumErrs = 500; % get at least 100 bit errors (more is better)
maxNum = 1e8; % OR stop if maxNum bits have been simulated
EbN0 = -1:8; % power efficiency range

% ======================================================================= %
% Other Options
% ======================================================================= %
% 

% ======================================================================= %
% Simulation Chain
% ======================================================================= %
BER = zeros(1, length(EbN0)); % pre-allocate a vector for BER results
berub = zeros(1,length(EbN0));

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


  trellis = poly2trellis(5,{'x4+x3+1','x4+x3+x+1'});
  
  G3 = [1,1,0,0,1;1,1,0,1,1];
 
  bits_enc1 = mod(conv(tx_bit,G3(1,:)),2);
  bits_enc2 = mod(conv(tx_bit,G3(2,:)),2);
  bits_enc = reshape([bits_enc1;bits_enc2],1,[]);
  bits_enc = bits_enc(1:2*N);
  

  % [MOD] symbol mapper
  tx_sym = modulation(bits_enc,'QPSK');

  % scatterplot(tx_sym);

  % [CHA] add Gaussian noise

  EbN0_value = 10^(EbN0(i)/10);
  R = length(tx_bit)/length(tx_sym);
  SNR = R * EbN0_value; 


  sigma2 = 1/(SNR*2);
  noise = sqrt(sigma2)*(randn(1,length(tx_sym)) + 1i*randn(1,length(tx_sym)));
  rx_data = tx_sym + noise;
  

  % soft receiver
  


  rx_dec = conv_soft_e3_qpsk(rx_data,G3);
  bits_out = rx_dec(1:N);

  % ===================================================================== %
  % End processing one block of information
  % ===================================================================== %
  
  BitErrs = length(find(bits_out ~= tx_bit)); % count the bit errors and evaluate the bit error rate
  totErr = totErr + BitErrs;
  num = num + N; 

  disp(['+++ ' num2str(totErr) '/' num2str(maxNumErrs) ' errors. '...
      num2str(num) '/' num2str(maxNum) ' bits. Projected error rate = '...
      num2str(totErr/num, '%10.1e') '. +++']);
  end 
  BER(i) = totErr/num; 


  N_d = 15;
  spect = distspec(trellis,N_d);
  df = spect.dfree;
  A_d = spect.event;
  w = spect.weight;
  
  for d_i = 1 : N_d
      berub(i) =  berub(i) + w(d_i)*qfunc(sqrt((df+d_i-1)*EbN0_value));

  end

end


% plots
figure(1)  
e = semilogy(EbN0, BER,'-^','DisplayName','coded transmission(e3)'); 
e.Color = '#EDB120';
e.LineWidth =2;
ylim([1e-4 1]);
hold on
grid on
f = semilogy(EbN0,berub,'--^','DisplayName','union bound(e3)');
f.Color = '#EDB120';
f.LineWidth =2;
legend
xlabel('EbN0(dB)');
ylabel('BER');
title('Eps 1&2&3 with QPSK Comparison');