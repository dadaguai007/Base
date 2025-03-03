% COHERENT DETECTION IN COLOURED NOISE USING A PREDICTIVE VITERBI ALGORITHM
% See Section 2.7 in K Vasudevan, "Digital Communications and Signal Processing,
% Second edition (CD-ROM included), Universities Press (Hyderabad), India.

clear all
close all
clc
SNR_dB = -10; % SNR PER BIT in dB
FRAME_SIZE = 10^3; % SYMBOLS IN EACH FRAME
NUM_BIT = 2*FRAME_SIZE; % NUMBER OF BITS IN EACH FRAME
NUM_FRAMES = 10^2; % SIMULATION RUNS
DECODING_DELAY = 20; % DECODING DELAY OF THE VITERBI ALGORITHM

% SNR PARAMETERS (OVERALL RATE IS 2)
SNR = 10^(0.1*SNR_dB); % SNR IN LINEAR SCALE
NOISE_VAR_1D = 0.5*2/(2*SNR); % 1D NOISE VARIANCE
NOISE_STD_DEV = sqrt(NOISE_VAR_1D); % NOISE STANDARD DEVIATION

% IIR FILTER PARAMETERS (USED TO GENERATE COLOURED NOISE)
a = 0.99;
B = sqrt(1-a^2); % INPUT COEFFICIENTS IN THE DIFFERENCE EQUATION
A = [1 -a]; % OUTPUT COEFFICIENTS IN THE DIFFERENCE EQUATION
AUTOCORR_SEQ = [NOISE_VAR_1D; a*NOISE_VAR_1D]; % AUTOCORRELATION OF COLOURED NOISE

% GENERATE PREDICTION FILTER COEFFICIENTS USING LEVINSON DURBIN ALGORITHM
PRED_COEF = Gen_Coef(AUTOCORR_SEQ,1);

%--------------------------------------------------------------------------
tic()
C_BER = 0; % CHANNEL ERRORS
for FRAME_CNT = 1:NUM_FRAMES
% SOURCE
DATA = randi([0 1],1,NUM_BIT);

% QPSK MAPPER
QPSK_SEQ = 1-2*DATA(1:2:end) + 1i*(1-2*DATA(2:2:end));

% ADDITIVE COLOURED GAUSSIAN NOISE
% AWGN
AWGN = normrnd(0,NOISE_STD_DEV,1,1000+FRAME_SIZE)+1i*normrnd(0,NOISE_STD_DEV,1,1000+FRAME_SIZE);
ACGN = filter(B,A,AWGN);
ACGN(1:1000) = []; % DISCARDING TRANSIENT SAMPLES

% CHANNEL OUTPUT
CHAN_OP = QPSK_SEQ + ACGN;

% PREDICTIVE VITERBI ALGORITHM BASED RECEIVER
QPSK_SYM = [1+1i 1-1i -1+1i -1-1i];
% BRANCH METRIC FOR THE PREDICTIVE VA
BRANCH_METRIC = zeros(16,FRAME_SIZE);
BRANCH_METRIC(1,2:end) = (abs((CHAN_OP(2:end)-QPSK_SYM(1))+(PRED_COEF*(CHAN_OP(1:end-1)-QPSK_SYM(1))))).^2;
BRANCH_METRIC(2,2:end) = (abs((CHAN_OP(2:end)-QPSK_SYM(2))+(PRED_COEF*(CHAN_OP(1:end-1)-QPSK_SYM(1))))).^2;
BRANCH_METRIC(3,2:end) = (abs((CHAN_OP(2:end)-QPSK_SYM(3))+(PRED_COEF*(CHAN_OP(1:end-1)-QPSK_SYM(1))))).^2;
BRANCH_METRIC(4,2:end) = (abs((CHAN_OP(2:end)-QPSK_SYM(4))+(PRED_COEF*(CHAN_OP(1:end-1)-QPSK_SYM(1))))).^2;

BRANCH_METRIC(5,2:end) = (abs((CHAN_OP(2:end)-QPSK_SYM(1))+(PRED_COEF*(CHAN_OP(1:end-1)-QPSK_SYM(2))))).^2;
BRANCH_METRIC(6,2:end) = (abs((CHAN_OP(2:end)-QPSK_SYM(2))+(PRED_COEF*(CHAN_OP(1:end-1)-QPSK_SYM(2))))).^2;
BRANCH_METRIC(7,2:end) = (abs((CHAN_OP(2:end)-QPSK_SYM(3))+(PRED_COEF*(CHAN_OP(1:end-1)-QPSK_SYM(2))))).^2;
BRANCH_METRIC(8,2:end) = (abs((CHAN_OP(2:end)-QPSK_SYM(4))+(PRED_COEF*(CHAN_OP(1:end-1)-QPSK_SYM(2))))).^2;

BRANCH_METRIC(9,2:end) = (abs((CHAN_OP(2:end)-QPSK_SYM(1))+(PRED_COEF*(CHAN_OP(1:end-1)-QPSK_SYM(3))))).^2;
BRANCH_METRIC(10,2:end) = (abs((CHAN_OP(2:end)-QPSK_SYM(2))+(PRED_COEF*(CHAN_OP(1:end-1)-QPSK_SYM(3))))).^2;
BRANCH_METRIC(11,2:end) = (abs((CHAN_OP(2:end)-QPSK_SYM(3))+(PRED_COEF*(CHAN_OP(1:end-1)-QPSK_SYM(3))))).^2;
BRANCH_METRIC(12,2:end) = (abs((CHAN_OP(2:end)-QPSK_SYM(4))+(PRED_COEF*(CHAN_OP(1:end-1)-QPSK_SYM(3))))).^2;

BRANCH_METRIC(13,2:end) = (abs((CHAN_OP(2:end)-QPSK_SYM(1))+(PRED_COEF*(CHAN_OP(1:end-1)-QPSK_SYM(4))))).^2;
BRANCH_METRIC(14,2:end) = (abs((CHAN_OP(2:end)-QPSK_SYM(2))+(PRED_COEF*(CHAN_OP(1:end-1)-QPSK_SYM(4))))).^2;
BRANCH_METRIC(15,2:end) = (abs((CHAN_OP(2:end)-QPSK_SYM(3))+(PRED_COEF*(CHAN_OP(1:end-1)-QPSK_SYM(4))))).^2;
BRANCH_METRIC(16,2:end) = (abs((CHAN_OP(2:end)-QPSK_SYM(4))+(PRED_COEF*(CHAN_OP(1:end-1)-QPSK_SYM(4))))).^2;

DEC_SYM = VITERBI_ALGORITHM(FRAME_SIZE,DECODING_DELAY,BRANCH_METRIC);    
% DEMAPPING TO QPSK SYMBOLS
DEC_QPSK_SEQ = QPSK_SYM(DEC_SYM);
DEC_DATA = zeros(1,2*(FRAME_SIZE-DECODING_DELAY));
DEC_DATA(1:2:end) = real(DEC_QPSK_SEQ)<0;
DEC_DATA(2:2:end) = imag(DEC_QPSK_SEQ)<0;

% BIT ERRORS IN EACH FRAME (IGNORING 1ST SYMBOL AND LAST TRANSIENT SYMBOLS)
C_BER = C_BER + nnz(DEC_DATA(3:end)-DATA(3:end-2*DECODING_DELAY));
end
% BIT ERROR RATE
 BER = C_BER/((NUM_BIT-2*DECODING_DELAY-2)*NUM_FRAMES)
 toc()
