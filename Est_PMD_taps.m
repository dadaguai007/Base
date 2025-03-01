function [Ncd, Npmd] = Est_PMD_taps(self, Rs, ros, lambda)
%% Estimated number of taps in DSP required to compensate for CD and PMD in coherent detection link
% Based on Ip, E., & Kahn, J. M. (2007). Digital equalization of chromatic dispersion and polarization mode dispersion.
% Journal of Lightwave Technology, 25(8), 2033-2043.
Ncd = 2*pi*abs(self.beta2(lambda))*self.L*Rs^2*ros; % eq (35)
Npmd = self.tauDGD*ros*Rs; % eq (36)
end

% 这段代码是一个函数，用于估计在相干检测链路中补偿色散（CD）和偏振模式色散（PMD）所需的数字信号处理（DSP）抽头数。这个估计是基于Ip和Kahn在2007年发表的一篇论文《Journal of Lightwave Technology》中的公式。
% 函数使用三个输入参数：
% - `Rs`：符号率（symbols per second），即每秒传输的符号数量。
% - `ros`：每符号的采样率（samples per symbol），即每个符号期间的采样次数。
% - `lambda`：传输信号的波长（以米为单位）。
% 函数返回两个输出：
% - `Ncd`：补偿色散所需的抽头数。
% - `Npmd`：补偿偏振模式色散所需的抽头数。
% 函数内部的工作原理如下：
% 1. `Ncd`的计算基于公式(35)在论文中，该公式用于估计补偿色散所需的抽头数。计算公式为：
%    \[ N_{cd} = 2\pi \left| \beta_2(\lambda) \right| L R_s^2 r_{os} \]
%    其中：
%    - `self.beta2(lambda)`是群速度色散（GVD）系数，它是一个与波长相关的函数。
%    - `self.L`是光纤链路的长度。
%    - `Rs`是符号率。
%    - `ros`是每符号的采样率。
% 2. `Npmd`的计算基于公式(36)，用于估计补偿偏振模式色散所需的抽头数。计算公式为：
%    \[ N_{pmd} = \tau_{DGD} r_{os} R_s \]
%    其中：
%    - `self.tauDGD`是差分群延迟（DGD），它是偏振模式色散的一个度量。
%    - `ros`是每符号的采样率。
%    - `Rs`是符号率。
% 这两个公式的输出`Ncd`和`Npmd`分别给出了在DSP中用于补偿色散和偏振模式色散所需的抽头数。这些抽头数可以用于设计数字均衡器，以在接收端补偿光纤信道引入的色散效应。
