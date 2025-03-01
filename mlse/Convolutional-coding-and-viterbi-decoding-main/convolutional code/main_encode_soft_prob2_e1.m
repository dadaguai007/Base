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


  trellis = poly2trellis(3,{'x2+1','x2+x+1'});
  
  G1 = [1,0,1;1,1,1];
 
  bits_enc1 = mod(conv(tx_bit,G1(1,:)),2);
  bits_enc2 = mod(conv(tx_bit,G1(2,:)),2);
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
 
  rx_dec = conv_soft_e1_qpsk(rx_data,G1);
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


  N_d = 20;
  spect = distspec(trellis,N_d);
  df = spect.dfree;
  e = spect.event;
  w = spect.weight;
  
 for d_i = 1 : N_d

     berub(i) =  berub(i) + w(d_i)*qfunc(sqrt((df+d_i-1)*EbN0_value));

 end


end



% plots
figure(1)  
a = semilogy(EbN0, BER,'-+','DisplayName','coded transmission(e1)'); 
a.Color = 	'#D95319';
a.LineWidth =2;
ylim([1e-4 1]);
hold on
grid on
b = semilogy(EbN0,berub,'--+','DisplayName','union bound(e1)');
b.Color = '#D95319';
b.LineWidth =2;
legend
xlabel('EbN0(dB)');
ylabel('BER');
% title('System 1: eps4 and AMPM');