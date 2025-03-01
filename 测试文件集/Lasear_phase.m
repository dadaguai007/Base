
%Laser激光器的建模

%在频域生成FM噪声功率谱密度（PSD）。
%使用IDFT将FM噪声PSD转换为时域。
%通过积分FM噪声来生成PM噪声。
%FM噪声PSD的长度在model.Lpsd中指定，或者在构造函数param.L中指定。
%除非L=1，否则最低频率分量（DC）处的噪声将被强制为零，这意味着除非L=1，否则相位噪声将是零均值。

%对于洛伦兹激光，需要指定linewidth。它会在linewidth/pi水平上创建一个平坦的频率功率谱密度（PSD），对应于FWHM为linewidth的洛伦兹光学光谱。
%对于半导体激光器，需要指定LFLW1GHZ、HFLW、fr、K和alpha。它根据这5个参数创建所需的半导体频率噪声PSD，并输出与所需PSD匹配的时域相位噪声序列。


% 频域中生成频率调制（FM）噪声功率谱密度（PSD）。
% 使用（IDFT）将 FM噪声PSD  转换为时域
% 将FM噪声PSD与白噪声进行卷积。
% 通过积分FM噪声生成相位调制（PM）噪声。

% 各种输入参数
%> carrier frequency  (Hz, can also be set by parameters)
Fc = const.c/1550e-9;
%> Lorentzian linewidth for standard (non-semiconductor) laser
linewidth;
%> Sampling frequency (Hz)
Fs;
%> Symbol rate (/sec)
Rs = 1;
%> Length of the noise signal (samples)
Lnoise;
%> Frequency Modulation noise
FMnoiseCal;
%> Phase Modulation noise
PMnoiseCal;

%选择使用哪种激光器的模型
%洛伦兹型，一个是半导体激光器
if isfield(param, 'linewidth') || ~isfield(param, 'LFLW1GHZ')
    disp('Using Lorentzian mode.');
    param_model.linewidth = paramdefault(param, 'linewidth', 100e3);
    param_model.Lpsd = 1;
else
    robolog('Using SCL mode.');
    param_model = param;
    param_model.Lpsd = paramdefault(param, {'L', 'Lir'}, 2^11);
end


%这段进行选择型号更为合适
% If linewidth is specified, lorenztial lineshape is assumed.
if isfield(param, 'linewidth') && ~isempty(param.linewidth)
    obj.linewidth = param.linewidth;
    obj.zoom = paramdefault(param, 'zoom', 2*obj.linewidth);
else
    obj.LFLW1GHZ=param.LFLW1GHZ;
    obj.HFLW=param.HFLW;
    obj.fr=param.fr;
    obj.K=param.K;
    obj.alpha=param.alpha;
    if obj.LFLW1GHZ > 50
        obj.zoom = paramdefault(param, 'zoom', 1e5*sqrt(obj.LFLW1GHZ));
    else
        obj.zoom = paramdefault(param, 'zoom', 2*obj.HFLW);
    end
end


%相应的频率轴，可以用来进行PSD的计算
% Generates the frequency axis
%生成频率轴，根据对象的绘图类型（'linear' 或 'log'）生成线性或对数刻度的频率轴。
function FMfreq = FMfreq(obj)
switch obj.type
    case 'linear'
        FMfreq=linspace(obj.MinFreq,obj.MaxFreq,obj.Lpsd)';
    case 'log'
        FMfreq = logspace(log10(obj.MinFreq),log10(obj.MaxFreq),obj.Lpsd)';
end
end


%得到pn之后，进行光场的计算
% if ~isnan(pn)
    % Create the complex baseband optical field
    %计算基带光场，相位噪声 pn
    %单位应该使用W
%     field = sqrt(obj.Power).*exp(1j*pn(1:Lnoise));
% 
% end

%产生噪声 phase noise
%> @brief save parameters, generate noise
function [fn,pn] = process(obj,varargin)
%> This function calls genNoise to generate
%> appropriate phase noise and frequency noise sequences.
%>

% The final time sequence will be mirrored in order to avoid
% discontinuities in case of repetition
obj.Lnoise = (obj.Lnoise + 2*obj.model.Lpsd)/2;
% Generate the desired frequency noise PSD
%频率噪声 PSD
[obj.FMnoiseCal, obj.PMnoiseCal] = obj.model.genPSD();
% Generate a gaussian time-domain sequence with desired PSD
[fn, pn] = obj.genNoise();
end

% 非常重要，要得到信号的PSD
%PSD的计算
% Generates the Frequency and Phase power spectral densities
%频率和相位功率谱密度，根据对象的属性生成频率噪声功率谱密度。
% 指定了 linewidth 属性，生成一个简单的洛伦兹型功率谱密度，否则调用 PSDmodel 方法生成半导体相位噪声功率谱密度。
function [FMnoiseCal, PMnoiseCal] = genPSD(obj)
% Build FMnoiseCal
FMfreq = obj.FMfreq();%生成频率轴
if ~isempty(obj.linewidth)
    %linewidth 属性存在且非空，假定采用洛伦兹线型
    FMnoiseCal = obj.linewidth/pi*ones(size(FMfreq));
else
    %linewidth 未指定，那么它会调用 PSDmodel 方法来生成半导体相位噪声功率谱密度
    FMnoiseCal=obj.PSDmodel(obj.FMfreq,obj.LFLW1GHZ,obj.HFLW,obj.K,obj.fr,obj.alpha);
end
if obj.f_cutoff ~= 0 && obj.f_cutoff>=obj.MinFreq
    %检查是否设置了 f_cutoff 属性，并且该值大于等于最小频率 MinFreq
    for a=1:length(FMfreq)
        %低于 f_cutoff 的频率部分的功率谱密度设为零
        if(FMfreq(a)<obj.f_cutoff)
            FMnoiseCal(a)=0;
        end
    end
end
%相位功率谱密度，将频率噪声功率谱密度 FMnoiseCal 除以频率的平方，乘以一个常数 (2*pi)^2。
PMnoiseCal = (2*pi)^2*FMnoiseCal./((2*pi)^2*FMfreq.^2);
end

% Spectral Lineshape Calculators
%用于获取光谱线形，它通过调用 CalcSpectrum 方法来计算频谱，并返回光谱线形及其对应的频率
function [SPEC,F] = getLineShape(obj, fmax)
NoLWPts = obj.Llw;
[SPEC,F]=obj.CalcSpectrum(obj.genPSD,obj.FMfreq,NoLWPts,fmax/2);
end
%获取指定点处的线宽。首先调用 getLineShape 方法获取光谱线形，然后找到与指定点相交的频率点，计算其线宽
function lw = getLinewidth(obj, point)
[SPEC,F] = getLineShape(obj, 15e6);
line = repmat(point, [length(SPEC), 1]);
intersects = intersections(F,line,F,pow2db(SPEC/max(SPEC)));
lw = max(diff(intersects));
end
%三种待调用的函数
%频率噪声功率谱密度
function FMnoise = PSDmodel(FMfreq,LFLW1GHZ,HFLW,K,fr,alpha)
%LFLW1GHZ 除以 π 并乘以 1e9（这是单位转换，将 GHz 转换为 Hz），然后除以频率轴 FMfreq。
% 计算得到了一个与频率轴有关的常数功率谱密度
FMnoise=LFLW1GHZ/pi*1e9./FMfreq;%低频
%添加了一个与频率无关的常数功率谱密度，表示高频段的噪声。
FMnoise=FMnoise+HFLW/pi*(1/(1+alpha^2));%高频
%计算了与频率有关的功率谱密度，表示了一个共振峰，其形状由这些参数决定。
FMnoise=FMnoise+HFLW/pi*(alpha^2/(1+alpha^2))...
    *fr^4./((fr^2-FMfreq.^2).^2+(K/2/pi)^2*fr^4*FMfreq.^2);
end
function[SPEC,F]=CalcSpectrum2(FMnoise,freq,NoLWPts,fmax)
%计算频谱，
% 根据输入的频率噪声 FMnoise、频率轴 freq、频谱点数 NoLWPts 和最大频率 fmax，计算出光谱 SPEC 和相应的频率轴 F。
%             akf = autocorr(FMnoise, NoLWPts);
df = fmax/NoLWPts;
for i=1:NoLWPts
    tau=(i-1)/df;
    akf(i) = exp(-2*(pi*tau)^2*trapz(freq,FMnoise.*sinc(pi*freq*tau).^2));
    %                F(i) = (i-1)*df;
end
SPEC = abs(fftshift(fft(akf)));
F = linspace(-fmax,fmax,NoLWPts);
end
function[SPEC,F]=CalcSpectrum(FMnoise,freq,NoLWPts,fmax)
%也计算频谱，采用不同的方法计算
%             FMnoise=db2pow(FMnoise);
df= 2*2*fmax/NoLWPts;
taumax=1/(2*df);
AKF=zeros(NoLWPts/2+1,1);
F=zeros(NoLWPts/2+1,1);
TAU=zeros(NoLWPts/2+1,1);
dtau=2*taumax/NoLWPts;
for i=1:NoLWPts/2
    tau = taumax-(i-1)*dtau;
    AKF(i)=exp(-pi*tau*PhaseNoiseModel_v1.EffLinewidth(1/tau, FMnoise, freq));
    TAU(i)=-tau;
    F(i)=(i-1)*df;
    %if not(round(tau-0.5)==round(tau-0.5-dtau))
    % round(tau);
    %end;
end
AKF(NoLWPts/2+1)=1;
Tau(NoLWPts/2+1)=0;
F(NoLWPts/2+1)=NoLWPts/2*df;
AKF=[AKF;AKF(end-1:-1:1)];
TAU=[TAU;-TAU(end-1:-1:1)];
SPEC=abs(fftshift(fft(AKF)));
F=[-flipud(F);F(2:end)];
end

function EffLW = EffLinewidth(BaudRate, FMnoise, freq)
%计算有效线宽度
x = freq;
tau=1/BaudRate;
%表示了频率噪声功率谱密度在频率轴上的加权。
y = FMnoise.*(sinc(freq*tau).^2);
%使用 trapz 函数对频率轴上的加权功率谱密度 y 进行积分
EffLW=2*(pi*tau)*trapz(x,y);
end






% Noise calculation methods
%计算了最大频率和最小频率
function val = MaxFreq(obj)
val = obj.Fs/2;
end

function val = MinFreq(obj)
val = obj.Fs/obj.model.Lpsd;
end




%Note:

% 频域中生成频率调制（FM）噪声功率谱密度（PSD）。
% 使用（IDFT）将 FM噪声PSD  转换为时域
% 将FM噪声PSD与白噪声进行卷积。
% 通过积分FM噪声生成相位调制（PM）噪声。
function [fn, pn] = genNoise(obj)
%生成频率和相位噪声的时间序列
% Time and frequency vectors (we double everything)
%频域量
freqs = [-flipud(obj.model.FMfreq);0; obj.model.FMfreq];
%噪声信号的长度除以采样率
time=linspace(-obj.model.Lpsd,obj.model.Lpsd,2*obj.model.Lpsd+1)/obj.Fs;
% Frequency domain filters  频域噪声模型 PSD
H_FN = obj.FMnoiseCal(:);
%频率域滤波器 H_FN
if length(H_FN)>1
    %PSD 包含多个频率点
    %包括正频率、零频率和负频率部分，以确保频域滤波器是实对称的
    %flipud(H_FN(:))：将 H_FN 倒序排列，即将正频率部分反转成负频率部分。
    % [ flipud(H_FN(:)); 0; H_FN(:)]：将倒序的 H_FN 和一个零值（对应零频率）以及原始的 H_FN 组合在一起，得到一个实对称的频率域滤波器。
    % /2：对滤波器的值进行除以2的操作，以确保频率响应保持一致。
    H_FN = [ flipud(H_FN(:)); 0; H_FN(:)]/2;
    %进行开方操作，以得到频域中的振幅响应。
    H_FN = sqrt(H_FN);
    %使用自定义的idft转换为时域
    % Time domain filters
    h_fn = obj.idft(freqs, H_FN, time);
    h_fn = real(h_fn(:));
else
    %一个频率点，表示频率噪声 PSD 仅有一个值，通常对应于洛伦兹线宽。
    h_fn = obj.Fs*sqrt(H_FN/2); %Hacked:
    %std(diff(pn)) = sqrt(2*pi*linewidth/obj.Fs)
    %H_FN = linewidth/pi
    %algebra based on lines 308 to 310
    %gives this equation
end
%时域滤波器 h_fn 功率谱的时域表示
% Frequency and phase noise generation
%白噪声
noise = awgn(zeros(2*obj.Lnoise,1),0)/sqrt(obj.Fs);
%卷积操作使用 h_fn 对白噪声进行滤波，得到频率噪声 fn。
%时域卷积，频域相乘
fn =conv(noise,h_fn,'valid');
%积分操作，相位噪声pn
pn = cumsum(2*pi*fn/obj.Fs);
%对相位噪声进行限制，确保在一定范围内
if obj.limitPn>0
    limit=[-1 1]*obj.limitPn;
    pn = pn-mean(pn);
    pn(pn>max(limit)) = -pn(pn>max(limit))+2*max(limit);
    pn(pn<min(limit)) = -pn(pn<min(limit))+2*min(limit);
end
%长度重新辅助
obj.Lnoise = 2*(obj.Lnoise - obj.model.Lpsd);
end




%自定义的idft
function x=idft(f,X,t)
% Compute IDFT (Inverse Discrete Fourier Transform) at times given
% in t, given frequency terms X taken at frequencies f:
% x(t) = sum { X(k) * e**(2*pi*j*t*f(k)) }
% k

shape = size(t);
f = f(:); % Format 'f' into a column vector
X = X(:); % Format 'X' into a column vector
t = t(:); % Format 't' into a column vector

df=diff(f);
df(length(df)/2)=0;
df=([0;df]+[df;0])/2;

%df(end)=0;
%df=f(2)-f(1);
%X(1)=X(1)/2;
%X(end)=X(end)/2;
x(1:length(t))=0;
x=x';
fspan=2*f(end);
for k=1:length(t)
    % It's just this simple:
    %x(k) = exp( 2*pi*1i * round(t(k)*fspan)/fspan*f')*(X.*df) ;
    x(k) = exp( 2*pi*1i * t(k)*f')*(X.*df) ;
end
x = reshape(x,shape);
end