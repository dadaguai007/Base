classdef SigClass < handle
%SigClass 定义了信号类，输入信噪比和相关参数后可以得到对应的SER、BER等属性
%description
%   SigClass v1.0
%   本代码定义了一种信号类
%   其中包含了信号的信噪比、调制格式、采样率、信号波特率、光波长、光偏振态 等参数
%   属性提取可以通过Dot完成，如Sig.BER
%   定义信号后可以完成SNR、OSNR间的转换，也可同时生成理论BER、SER
%   目前代码中只可以定义QAM、PAM、PSK三种信号
%   本代码的计算公式来源于《数字通信（第五版）》中文版教材中4.2节-4.3节
%   如有错误或不妥之处，欢迎指正！
%   首次编辑于2023.10.18                         作者：翟孟晟
%   请依次输入：OSNR/SNR/OSNR_Line/SNR_Line（输入选项），信噪比，调制模式，阶数，*sps，*波特率，*波长，*偏振(*为可缺省项）
%   使用实例：SigClass('SNR_Line',SNR_Line,'QAM',4)  或  SigClass('SNR',SNR_dB,'QAM',4)
%   若 Sig=SigClass('SNR_Line',SNR_Line,'QAM',4),则 BER=Sig.BER;  SER=Sig.SER;

    properties      %属性设置
        SNR                 %信噪比
        OSNR                %光信噪比
        EsNo                %符号信噪比
        EbNo                %比特信噪比
        M                   %调制阶数
        k                   %每个符号中的比特数
        Model               %调制格式
        SER                 %符号错误率
        BER                 %比特错误率
        Rb                  %比特率
        Rs                  %波特率
        Fs                  %采样率
        sps=1               %=Fs/Rs,默认采样率=波特率
        lambda=1550e-9      %光波长，默认波长1550nm
        p=1                 %偏振态，默认偏振模式为单偏振
    end
    methods         %内置方法
        function obj = SigClass(inputmod,SNR,Model,M,sps,Baud,lambda,p)
        %阐明输入数值含义、调制格式、信号属性后方可计算
        if inputmod=="OSNR"||inputmod=="SNR"         %检验输入为SNR或OSNR(dB)
            switch nargin
                case {4,5}   %输入为SNR
                    if inputmod=="SNR"
                        if nargin==5
                            obj.sps = sps;
                        end
                        obj.M = M;
                        obj.k=log2(obj.M);
                        obj.Model=Model;
                        obj.SNR = dB2Line(SNR);
                        obj.EsNo=obj.SNR*obj.sps;
                        obj.EbNo=obj.EsNo/obj.k;
                        Checksig(obj)
                    elseif inputmod=="OSNR"
                        error('参数输入有误，请检查信号波特率是否已输入');
                    else
                        error('参数输入有误，请检查信噪比输入');
                    end
                case {6,7,8} %输入为OSNR
                    obj.M = M;
                    obj.k=log2(M);
                    obj.Model=Model;

                      switch nargin
                        case 6
                        obj.Rs=Baud;
                        obj.Rb=Baud*obj.k;
                        obj.Fs=obj.Rs*obj.sps;
                        case 7
                        obj.Rs=Baud;
                        obj.Rb=Baud*obj.k;
                        obj.Fs=obj.Rs*obj.sps;
                        obj.lambda=lambda;
                        case 8
                        obj.Rs=Baud;
                        obj.Rb=Baud*obj.k;
                        obj.Fs=obj.Rs*obj.sps;
                        obj.lambda=lambda;                   
                        obj.p=p;
                      end
                      c0=299792458;
                      d=1e-10;
                      Bref=c0/obj.lambda^2*d;
                    if inputmod=="OSNR"
                        obj.OSNR = dB2Line(SNR);%输入为OSNR
                        obj.SNR=2*Bref/(obj.p*obj.Fs)*obj.OSNR;%SNR与OSNR转换
                    elseif inputmod=="SNR"
                        obj.SNR = dB2Line(SNR);
                        obj.OSNR = obj.SNR./(2*Bref/(obj.p*obj.Fs));
                    end
                    obj.EsNo=obj.SNR*sps;
                    obj.EbNo=obj.SNR*sps/obj.k;
                    Checksig(obj)
                otherwise
                    error('输入参数有误，请重新输入')
            end
            %线性转dB
            obj.SNR=Line2dB(obj.SNR);
            obj.OSNR=Line2dB(obj.OSNR);
            obj.EsNo=Line2dB(obj.EsNo);
            obj.EbNo=Line2dB(obj.EbNo);
        elseif  inputmod=="OSNR_Line"||inputmod=="SNR_Line"%检验输入为SNR或OSNR(Line)
            switch nargin
                case {4,5}   %输入为SNR
                    if inputmod=="SNR_Line"
                        if nargin==5
                            obj.sps = sps;
                        end
                        obj.M = M;
                        obj.k=log2(obj.M);
                        obj.Model=Model;
                        obj.SNR = SNR;
                        obj.EsNo=obj.SNR*obj.sps;
                        obj.EbNo=obj.EsNo/obj.k;
                        Checksig(obj)
                    elseif inputmod=="OSNR_Line"
                        error('参数输入有误，请检查信号波特率是否已输入');
                    else
                        error('参数输入有误，请检查信噪比输入');
                    end
                case {6,7,8} %输入为OSNR
                    obj.M = M;
                    obj.k=log2(M);
                    obj.Model=Model;
                      switch nargin %输入具体情况
                        case 6
                        obj.Rs=Baud;
                        obj.Rb=Baud*obj.k;
                        obj.Fs=obj.Rs*obj.sps;
                        case 7
                        obj.Rs=Baud;
                        obj.Rb=Baud*obj.k;
                        obj.Fs=obj.Rs*obj.sps;
                        obj.lambda=lambda;
                        case 8
                        obj.Rs=Baud;
                        obj.Rb=Baud*obj.k;
                        obj.Fs=obj.Rs*obj.sps;
                        obj.lambda=lambda;                   
                        obj.p=p;
                      end
                    c0=299792458;
                    d=1e-10;
                    Bref=c0/obj.lambda^2*d;
                    if inputmod=="OSNR_Line"
                        obj.OSNR = SNR;%输入为OSNR
                        obj.SNR=2*Bref./(obj.p*obj.Fs)*obj.OSNR;%SNR与OSNR转换
                    elseif inputmod=="SNR_Line"
                        obj.SNR = SNR;%输入为SNR
                        obj.OSNR = obj.SNR./(2*Bref/(obj.p*obj.Fs));%SNR与OSNR转换
                    end
                    obj.EsNo=obj.SNR*sps;
                    obj.EbNo=obj.SNR*sps/obj.k;
                    Checksig(obj)
                otherwise
                    error('输入参数有误，请重新输入')
           end
            
        else
            error('参数输入有误，请检查信噪比输入');
        end     %if_Check end
        end     %function SigClass end
%% 调制模式分流Checksig
        function Checksig(obj)
            %调制模式分流
            switch(obj.Model)
                case {'QAM'}
                    SER_QAM(obj)
                case {'PSK'}
                    SER_PSK(obj)
                case {'PAM'}
                    SER_PAM(obj)
            end
            BER_T(obj)
        end     %function Checksig end
%% QAM相关计算SER_QAM
        function SER_QAM(obj)
            % QAM相关计算
            P_half=2*(1-1/sqrt(obj.M)).*qfunc(sqrt(3*log2(obj.M).*obj.EbNo/(obj.M-1)));
            obj.SER=2*P_half.*(1-P_half/2);
        end
%% PAM相关计算SER_PAM
        function SER_PAM(obj)
            % PAM相关计算
            obj.SER=2*(1-1/obj.M)*qfunc(sqrt(6*obj.k.*obj.EbNo/(obj.M^2-1)));
        end
%% PSK相关计算SER_PSK
        function SER_PSK(obj)
            % PSK相关计算
            dv=0.01;d_theta=pi/1000;
            d1=-pi:d_theta:pi;
            V1=0:dv:1000;
            format long
            P_theta=nan(length(d1));
%             syms V thetar
            for i=1:length(d1)
                P_V_theta=@(V,thetar)exp(-(V-sqrt(obj.EsNo*2)*cos(thetar)).^2/2).*V;
                P_V_theta1=P_V_theta(V1,d1(i));
                P=trapz(V1,P_V_theta1);
                P_theta(i)=exp(-obj.EsNo*sin(d1(i)).^2)/pi/2.*P;
            end
            d2=-pi/obj.M:d_theta:pi/obj.M;
            st=find(abs(d1-d2(1))==min(abs(d1-d2(1))));nd=find(abs(d1-d2(end))==min(abs(d1-d2(end))));
            obj.SER=1-sum(P_theta(st:nd).*d_theta);
        end
%% SER、BER转化   
        function BER_T(obj)
            % SER、BER转化
            obj.BER=1-2.^(log2(1-obj.SER)/log2(obj.M));
        end
    end
end
%% 线性与dB互转
function dB=Line2dB(m)
    % 线性转dB
    dB=10.*log10(m);
    return
end
function Line=dB2Line(m)
    % dB转线性
    Line=10.^(m/10);
    return
end