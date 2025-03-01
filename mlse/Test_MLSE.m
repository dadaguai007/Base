% Test MLSE PAM-4 or OOK
clc;clear;close all;
addpath("D:\PhD\Codebase\");
type='PAM2';
% 参数
sps = 2;
Rs = 510e9;
Fs = Rs*sps;
if strcmp(type,'PAM2')
    % 信号生成
    M=2;
    msgLen = 8000;
    data_2bit=randi([0,1],log2(M),msgLen);
    % 相当于四个电平
    symbols = 2.^(0:log2(M)-1)*data_2bit;
    symbols=symbols.';
    % Mapeia bits para pulsos eletricos
    symbTx = pammod(symbols,M,0,'gray');
    % symbTx = pnorm(symbTx);
    % Pulso
    hsqrt = rcosdesign(0.01,256,sps,'sqrt');

    % Upsampling
    symbolsUp = upsample(symbTx,sps);
    % % pulse shaping
    sig=conv(symbolsUp,hsqrt,'same');

    % 星座点
    const = pammod([0:M-1],M);

    % MLSE 回溯长度
    tblen =  10;

    % 信道
    chanest = [0.986; 0.845; 0.237; 0.12345+0.31i];
    msgFilt = filter(chanest,1,sig);

    eqSym = mlseeq(msgFilt,chanest,const,tblen,'rst',sps);
    eqMsg = pamdemod(eqSym,M);

    tran_sig=downsample(msgFilt,sps);
    eqMsg_tran = pamdemod(tran_sig,M);


    [nerrs_mlse ber_mlse] = biterr(data_2bit.', eqMsg)
    [nerrs ber] = biterr(data_2bit.', eqMsg_tran)

elseif strcmp(type,'PAM4')

    % 信号生成
    M=4;
    msgLen = 8000;
    data_2bit=randi([0,1],log2(M),msgLen);
    % 相当于四个电平
    symbols = 2.^(0:log2(M)-1)*data_2bit;
    symbols=symbols.';
    % Mapeia bits para pulsos eletricos
    symbTx = pammod(symbols,M,0,'gray');
    % 参考
    label=double(symbTx);

    % Pulso
    hsqrt = rcosdesign(0.01,256,sps,'sqrt');

    % Upsampling
    symbolsUp = upsample(symbTx,sps);
    % % pulse shaping
    sig=conv(symbolsUp,hsqrt,'same');

    % 星座点
    const = pammod([0:M-1],M);

    % MLSE 回溯长度
    tblen =  100;

    % 信道
    chanest = [0.986; 0.845; 0.237; 0.12345+0.31i];
    msgFilt = filter(chanest,1,sig);

    eqSym = mlseeq(msgFilt,chanest,const,tblen,'rst',sps);

    Chan_sig=downsample(msgFilt,sps);
    % 重新量化
    A1=[-2 0 2];
    % 参考序列
    [~,tran_sig] = quantiz(real(Chan_sig),A1,[-3,-1,1,3]);

    ncut=1e3;
    [ber_mlse,num_mlse,~] = CalcBER(eqSym(ncut:end),label(ncut:end)); %计算误码率
    fprintf('After Num of Errors = %d, BER = %1.7f\n',num_mlse,ber_mlse);


    [ber,num,~] = CalcBER(tran_sig(ncut:end),label(ncut:end)); %计算误码率
    fprintf('Num of Errors = %d, BER = %1.7f\n',num,ber);

end

% figure;
% plot((sig))