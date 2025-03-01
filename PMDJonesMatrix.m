function [H,N] = PMDJonesMatrix(omega,param)
% self include: L,PMD_span_length,meanDGD
% meanDGDps = 0.1; % mean DGD  (ps/sqrt(km))
% PMD_span_length = 1e3;  % Section length for simulating PMD (m)
tauDGD = param.meanDGDps*1e-12*sqrt(param.L/1e3);

%总长度除去每段长度
Nsect = ceil(param.L/param.PMD_span_length);
%每段的微分群时延
dtau = tauDGD/sqrt(Nsect);


%创建一个空的Jones矩阵M，维度为2x2xlength(omega)，即两个极化状态，长度为omega
M = zeros(2,2,length(omega));
%通过随机旋转矩阵来初始化第一个omega值的Jones矩阵M(:,:,1)，初始偏振态
M(:, :, 1) = randomRotationMatrix(); % rotation to an arbritary polarizataion state
N = {};
H=eye(2);
H = repmat(H, [1, 1, length(omega)]);
for k = 1:Nsect
    U = randomRotationMatrix();

    for m = 2:length(omega)
        Dw = [exp(1j*dtau*omega(m)/2), 0; 0, exp(-1j*dtau*omega(m)/2)]; % Birefringence matrix
        M(:,:,m) = U*Dw*U';
    end
    N = [N, {M}];
    for m=2:length(omega)
    H(:,:,m)=H(:,:,m)*M(:,:,m);
    end
end


    function U = randomRotationMatrix()
        %函数生成随机的旋转矩阵（U）
        %三参量模型
        phi = rand(1, 3)*2*pi;
        U1 = [exp(-1j*phi(1)/2), 0; 0 exp(1j*phi(1)/2)];
        U2 = [cos(phi(2)/2) -1j*sin(phi(2)/2); -1j*sin(phi(2)/2) cos(phi(2)/2)];
        U3 = [cos(phi(3)) -sin(phi(3)); sin(phi(3)) cos(phi(3))];

        U = U1*U2*U3;
    end
end
% % 创建一个2x2的单位矩阵
% unitMatrix = eye(2);
% 
% % 指定复制的次数，这里是在每个维度上复制3次
% repeatedMatrix = repmat(unitMatrix, [1, 1, 3]);
% 
% % 显示结果
% disp(repeatedMatrix);
