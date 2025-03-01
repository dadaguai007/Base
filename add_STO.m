function y_STO=add_STO(y,nSTO,Type,sps)
%应该分为两种情况，一种基频，一种已经进行多采样
%注意，需要将信号进行上采样后，再进行STO的处理
%步骤如下：
% sps=10;
%输入向量必须为行向量
if isrow(y)
    % 输入向量已经是行向量，无需更改
    y = y;
else
    % 输入向量不是行向量，转换为行向量
    y = y';
end


if nargin < 4
    switch lower(Type)
        case 'zero'
            if nSTO >= 0
                %补零的方式
                y_STO = [y(nSTO+1:end), zeros(1, nSTO)];
            else
                y_STO = [zeros(1, -nSTO)', y(1:end+nSTO)];
            end
        case 'cycle'
            if isreal(y)
                % 实数信号平移
                y_STO = circshift(y, nSTO);  % 实部平移
            else
                %复数信号的平移
                %nSTO为负值时，信号向右平移
                real_part_shifted = circshift(real(y), nSTO);  % 实部平移
                imaginary_part_shifted = circshift(imag(y), nSTO);  % 虚部平移
                y_STO = complex(real_part_shifted, imaginary_part_shifted);  % 平移后的复数信号
            end
    end

else


    l=1:1:length(y);
    l2=1:1/sps:length(y);
    z = interp1(l,y,l2,"spline"); % 对接收信号进行插值
    switch lower(Type)
        case 'zero'
            if nSTO >= 0
                %补零的方式
                y_STO = [z(nSTO+1:end), zeros(1, nSTO)];
            else
                y_STO = [zeros(1, -nSTO), z(1:end+nSTO)];
            end
        case 'cycle'
            if isreal(z)
                % 实数信号平移
                y_STO = circshift(z, nSTO);  % 实部平移
            else
                %复数信号的平移
                %nSTO为负值时，信号向右平移
                real_part_shifted = circshift(real(z), nSTO);  % 实部平移
                imaginary_part_shifted = circshift(imag(z), nSTO);  % 虚部平移
                y_STO = complex(real_part_shifted, imaginary_part_shifted);  % 平移后的复数信号
            end
    end
    %下采样输出
    y_STO=downsample(y_STO,sps);
end
end