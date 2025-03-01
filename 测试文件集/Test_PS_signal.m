% test ps signal
clc;clear;close all;
N_sym=4e5;
M=16;
C=3.5;
N=log2(M);
% 获得 h 和 v 比特的起始位置
h_bit=N-1;
v_bit=floor(N/2)-1;

% 输出两个概率
[p0, p1] = qam_optimal_p0(C);
vh_bits=randi([0,1],2*N_sym,1);
% p0就是比特0出现的概率,p1为1出现的概率
% 现在的bug就是 假设0出现的概率，需要将概率对调以下才能显示出需要的状态
amp_bits=double((rand((N-2)*N_sym,1)>p1));
figure;
histogram(amp_bits)

% 垂直比特的开始索引
v_bits_ind=(v_bit:N:N_sym*N);
% 水平比特的开始索引
h_bits_ind=(h_bit:N:N_sym*N);
% amp_bits_ind
X=(1:N_sym*N).';
Y=horzcat(v_bits_ind,h_bits_ind).';
% 两个数据集的差集
% 幅度信息的索引
amp_bits_ind=setdiff(X,Y);
% 比特映射长度
bit_stream_mapping=zeros(N_sym*N,1);
% 一半长度
half_size = floor(length(vh_bits) / 2);

% 映射到不同部分
bit_stream_mapping(v_bits_ind) = vh_bits(1:half_size);
bit_stream_mapping(h_bits_ind) = vh_bits(half_size+1:end);
bit_stream_mapping(amp_bits_ind) = amp_bits;
% bit_stream_mapping(H)=amp_bits;
% bit_stream_mapping(amp_bits_ind)=vh_bits;
% 到这里，依旧是比特0占据主导。
bit_stream=reshape(bit_stream_mapping,N,[]).';
% bit_stream=reshape(bit_stream_mapping,N,[]);
jj=num2str(bit_stream);
% 
for i = 1:N_sym

g(i,:)=bin2dec(jj(i,:));
end
figure;
histogram(g)
S=qammod(g,M,'PlotConstellation',true);
% S=qammod(bit_stream_mapping,M,'InputType','bit');
S=awgn(S,18,'measured');
scatterplot(S)
%% ccdm
clc;clear;close all;
N_sym=4e5;
M=16;
C=3.5;
N=log2(M);
% 获得 h 和 v 比特的起始位置
h_bit=N-1;
v_bit=floor(N/2)-1;


[p0, p1] = qam_optimal_p0(C);
%p1 将合适的1进行概率赋值。 p0 合适的0进行概率赋值
result=ccdm_make_cfg_probability(32, 48, p1);
% 将CCDM进行验证
%调制信号长度
N_sym = N_sym + (result.c_sz - mod(N_sym,result.c_sz));
vh_bits=randi([0,1],2*N_sym,1);
% 强度比特信号长度调整
amp_bits_nc=randi([0,1],floor(N_sym*(N-2)*result.d_sz/result.c_sz),1);
%编码
amp_bits= ccdm_encode(amp_bits_nc,result.one_bits,result.d_sz,result.c_sz);

figure;
histogram(amp_bits)
% 垂直比特的开始索引
v_bits_ind=(v_bit:N:N_sym*N);
% 水平比特的开始索引
h_bits_ind=(h_bit:N:N_sym*N);
% amp_bits_ind
X=(1:N_sym*N).';
Y=horzcat(v_bits_ind,h_bits_ind).';
% 两个数据集的差集
% 幅度信息的索引
amp_bits_ind=setdiff(X,Y);
% 比特映射长度
bit_stream_mapping=zeros(N_sym*N,1);
% 一半长度
half_size = floor(length(vh_bits) / 2);

% 映射到不同部分
bit_stream_mapping(v_bits_ind) = vh_bits(1:half_size);
bit_stream_mapping(h_bits_ind) = vh_bits(half_size+1:end);
bit_stream_mapping(amp_bits_ind) = amp_bits;

S=qammod(bit_stream_mapping,M,'InputType','bit');
S=awgn(S,15,'measured');
scatterplot(S)

%% ccdm 解码：
bit_stream_demap = qamdemod(S, M,'OutputType','bit');
%水平比特和垂直比特
v_bits = bit_stream_demap(v_bits_ind);
h_bits = bit_stream_demap(h_bits_ind);
vh_bits = horzcat([v_bits, h_bits]);
% 强度比特
amp_bit = bit_stream_demap(amp_bits_ind);
% ccdm decode
ccdm_amp_bits_nc = ccdm_decode(amp_bit, result.one_bits, result.d_sz, result.c_sz);
%ccdm_amp_bits_nc.'-amp_bits_nc;
ccdm_amp_bits_nc =ccdm_amp_bits_nc.';

[ber,nErr]=CalcBER(ccdm_amp_bits_nc,amp_bits_nc);
fprintf('误码率为：%0.3f\n',ber)