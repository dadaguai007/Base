classdef SigClass < handle
%SigClass �������ź��࣬��������Ⱥ���ز�������Եõ���Ӧ��SER��BER������
%description
%   SigClass v1.0
%   �����붨����һ���ź���
%   ���а������źŵ�����ȡ����Ƹ�ʽ�������ʡ��źŲ����ʡ��Ⲩ������ƫ��̬ �Ȳ���
%   ������ȡ����ͨ��Dot��ɣ���Sig.BER
%   �����źź�������SNR��OSNR���ת����Ҳ��ͬʱ��������BER��SER
%   Ŀǰ������ֻ���Զ���QAM��PAM��PSK�����ź�
%   ������ļ��㹫ʽ��Դ�ڡ�����ͨ�ţ�����棩�����İ�̲���4.2��-4.3��
%   ���д������֮������ӭָ����
%   �״α༭��2023.10.18                         ���ߣ�������
%   ���������룺OSNR/SNR/OSNR_Line/SNR_Line������ѡ�������ȣ�����ģʽ��������*sps��*�����ʣ�*������*ƫ��(*Ϊ��ȱʡ�
%   ʹ��ʵ����SigClass('SNR_Line',SNR_Line,'QAM',4)  ��  SigClass('SNR',SNR_dB,'QAM',4)
%   �� Sig=SigClass('SNR_Line',SNR_Line,'QAM',4),�� BER=Sig.BER;  SER=Sig.SER;

    properties      %��������
        SNR                 %�����
        OSNR                %�������
        EsNo                %���������
        EbNo                %���������
        M                   %���ƽ���
        k                   %ÿ�������еı�����
        Model               %���Ƹ�ʽ
        SER                 %���Ŵ�����
        BER                 %���ش�����
        Rb                  %������
        Rs                  %������
        Fs                  %������
        sps=1               %=Fs/Rs,Ĭ�ϲ�����=������
        lambda=1550e-9      %�Ⲩ����Ĭ�ϲ���1550nm
        p=1                 %ƫ��̬��Ĭ��ƫ��ģʽΪ��ƫ��
    end
    methods         %���÷���
        function obj = SigClass(inputmod,SNR,Model,M,sps,Baud,lambda,p)
        %����������ֵ���塢���Ƹ�ʽ���ź����Ժ󷽿ɼ���
        if inputmod=="OSNR"||inputmod=="SNR"         %��������ΪSNR��OSNR(dB)
            switch nargin
                case {4,5}   %����ΪSNR
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
                        error('�����������������źŲ������Ƿ�������');
                    else
                        error('�������������������������');
                    end
                case {6,7,8} %����ΪOSNR
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
                        obj.OSNR = dB2Line(SNR);%����ΪOSNR
                        obj.SNR=2*Bref/(obj.p*obj.Fs)*obj.OSNR;%SNR��OSNRת��
                    elseif inputmod=="SNR"
                        obj.SNR = dB2Line(SNR);
                        obj.OSNR = obj.SNR./(2*Bref/(obj.p*obj.Fs));
                    end
                    obj.EsNo=obj.SNR*sps;
                    obj.EbNo=obj.SNR*sps/obj.k;
                    Checksig(obj)
                otherwise
                    error('���������������������')
            end
            %����תdB
            obj.SNR=Line2dB(obj.SNR);
            obj.OSNR=Line2dB(obj.OSNR);
            obj.EsNo=Line2dB(obj.EsNo);
            obj.EbNo=Line2dB(obj.EbNo);
        elseif  inputmod=="OSNR_Line"||inputmod=="SNR_Line"%��������ΪSNR��OSNR(Line)
            switch nargin
                case {4,5}   %����ΪSNR
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
                        error('�����������������źŲ������Ƿ�������');
                    else
                        error('�������������������������');
                    end
                case {6,7,8} %����ΪOSNR
                    obj.M = M;
                    obj.k=log2(M);
                    obj.Model=Model;
                      switch nargin %����������
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
                        obj.OSNR = SNR;%����ΪOSNR
                        obj.SNR=2*Bref./(obj.p*obj.Fs)*obj.OSNR;%SNR��OSNRת��
                    elseif inputmod=="SNR_Line"
                        obj.SNR = SNR;%����ΪSNR
                        obj.OSNR = obj.SNR./(2*Bref/(obj.p*obj.Fs));%SNR��OSNRת��
                    end
                    obj.EsNo=obj.SNR*sps;
                    obj.EbNo=obj.SNR*sps/obj.k;
                    Checksig(obj)
                otherwise
                    error('���������������������')
           end
            
        else
            error('�������������������������');
        end     %if_Check end
        end     %function SigClass end
%% ����ģʽ����Checksig
        function Checksig(obj)
            %����ģʽ����
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
%% QAM��ؼ���SER_QAM
        function SER_QAM(obj)
            % QAM��ؼ���
            P_half=2*(1-1/sqrt(obj.M)).*qfunc(sqrt(3*log2(obj.M).*obj.EbNo/(obj.M-1)));
            obj.SER=2*P_half.*(1-P_half/2);
        end
%% PAM��ؼ���SER_PAM
        function SER_PAM(obj)
            % PAM��ؼ���
            obj.SER=2*(1-1/obj.M)*qfunc(sqrt(6*obj.k.*obj.EbNo/(obj.M^2-1)));
        end
%% PSK��ؼ���SER_PSK
        function SER_PSK(obj)
            % PSK��ؼ���
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
%% SER��BERת��   
        function BER_T(obj)
            % SER��BERת��
            obj.BER=1-2.^(log2(1-obj.SER)/log2(obj.M));
        end
    end
end
%% ������dB��ת
function dB=Line2dB(m)
    % ����תdB
    dB=10.*log10(m);
    return
end
function Line=dB2Line(m)
    % dBת����
    Line=10.^(m/10);
    return
end