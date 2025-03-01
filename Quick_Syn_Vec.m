% [DeWaveform,P,OptSampPhase]=Quick_Syn_Vec(RecWaveform,label,SampleT_Osci,SampleT_AWG,OptSampPhase)
% Origin version was got from Jane
% log: 2014.4.2 halfLen & quadLen are added to replace original fixed value
% log: 20170619 Vectorize the correlation calculation part to have 75%
%      running time reduction

function [DeWaveform,P,OptSampPhase,MaxCorrIndex]=Quick_Syn_Vec(RecWaveform,label,SampleT_Osci,SampleT_AWG,OptSampPhase)
%RecWaveform：接收到的波形数据
%label：一个参考波形
%SampleT_Osci：示波器采样时间间隔。
%SampleT_AWG：任意波形发生器（AWG）的采样时间间隔。
%OptSampPhase：可选的最优采样相位。
if nargin < 5
    OptSampPhase = [];
end

% find the optimal sampe phase
if isempty(OptSampPhase)
    seq_len = 1e5;
    % seq_len = ceil(length(RecWaveform)/2);
    SeqForSync = RecWaveform(1:seq_len);
    
    maxPhaseError = 0.005;
    if 1/maxPhaseError - floor(1/maxPhaseError) ~= 0
        error('Wrong maxPhaseError value!');
    end
    % 时间序列
    tSeqBefSamp=(0:numel(SeqForSync)-1).'*SampleT_Osci;
    %搜索的相位
    SampPhaseVec=(0:maxPhaseError:(1-maxPhaseError)).';
    % 重采样后的时间点，搜索最小相位值 与 信号点数间隔时间  的乘积
    % 对SeqForSync按不同相位重采样,得到AftSampWaveform
    % 实现了对信号的重采样,通过改变时间轴,可以获得不同相位下的采样信号。
    
    % AWG下信号的时间索引
    tSeqAftSamp = 0:maxPhaseError*SampleT_AWG:tSeqBefSamp(end);
    AftSampWaveform = interp1(tSeqBefSamp,SeqForSync,tSeqAftSamp,'spline');

    %将重采样后的信号分割成多个符号波形，并对这些符号波形进行处理。这样做可以提高后续处理的效率，例如在数字通信系统中，每个符号可以单独进行解码和错误检测。
    % 通过这种方式，信号处理可以更加灵活和高效，同时也有助于提高信号的同步质量和性能。
    
    %根据maxPhaseError将AftSampWaveform划分为nSymbols个符号波形,每个符号波形对应一个采样相位。
    nSymbols = floor(length(AftSampWaveform)*maxPhaseError);
    AftSampWaveform = reshape(AftSampWaveform(1:nSymbols/maxPhaseError),1/maxPhaseError,[]).';
    % use "xcorr is equivalent to conv with a vector's flip"
    % 计算AftSampWaveform与本地参考信号label的卷积CorrMat。
    % 互相关是相当于数组倒置的卷积（互相关需要减去直流）
    CorrMat = convn(AftSampWaveform,flipud(label));
    CorrResult = max(abs(CorrMat)).';
    
    [~, MaxCorrIndex]=max(CorrResult);
%     拟合求取最佳采样相位而实现的二次曲线拟合。

    if MaxCorrIndex==1
        polyForQ=polyfit([SampPhaseVec(end);SampPhaseVec(1:2)+1],[CorrResult(end);CorrResult(1:2)],2);
        % warning('Optimum resampling phase near phase Zero');
    elseif MaxCorrIndex==numel(CorrResult)
        polyForQ=polyfit([SampPhaseVec(end-1:end);SampPhaseVec(1)+1],[CorrResult(end-1:end);CorrResult(1)],2);
        % warning('Optimum resampling phase near phase One');
    else
        polyForQ = polyfit(SampPhaseVec(MaxCorrIndex-1:MaxCorrIndex+1),CorrResult(MaxCorrIndex-1:MaxCorrIndex+1), 2);
    end

    % 拟合出CorrResult的峰值点相位,作为最佳采样相位OptSampPhase。
    %  最大互相关对应的相位OptSampPhase,这就是最佳采样相位
    %   一个是多项式微分，一个是求解多项式；多项式计算
    OptSampPhase=roots(polyder(polyForQ));
    OptCorr=polyval(polyForQ,OptSampPhase);
    OptSampPhase=OptSampPhase-floor(OptSampPhase);
    %对 OptSampPhase 进行了向下取整操作，确保它是一个整数。
    % 这是必要的，因为采样相位通常表示为一个离散的整数，而不是一个实数
end
% 用OptSampPhase对整个RecWaveform信号进行重采样,得到重采样后的信号AftSampWaveform。
% 计算AftSampWaveform与label的互相关Corr,找到主峰P,这表示找到了同步位置。
% 得到同步后的信号DeWaveform
tSeqBefSamp=(0:numel(RecWaveform)-1)*SampleT_Osci;
tSeqBefSamp=tSeqBefSamp.';
tSeqAftSamp=OptSampPhase*SampleT_AWG:SampleT_AWG:tSeqBefSamp(end);
tSeqAftSamp=tSeqAftSamp.';
AftSampWaveform=interp1(tSeqBefSamp,RecWaveform,tSeqAftSamp,'spline');

% 同步
Corr=abs(xcorr(label,AftSampWaveform-mean(AftSampWaveform)));
%在一个【】中赋值两个索引，都赋值为零
Corr([1:1000,length(AftSampWaveform)-1000+1:length(AftSampWaveform)+1000]) = 0;
figure(200);
plot(Corr)
close(200);
peak = max(Corr);
P = find(Corr>0.9*peak);
%翻转
P = flipud(numel(AftSampWaveform)-P);
P(P<1) = [];
%移除了 P 中所有小于1的索引。这是为了确保所有的同步位置都是有效的，因为索引小于1的位置上不会有有效的采样点
DeWaveform=AftSampWaveform;

end
