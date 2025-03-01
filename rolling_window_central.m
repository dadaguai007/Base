
function data_out = rolling_window_central(data, size, wrap)

%    将一个一维数组重新排列成一个有重叠帧的二维数组【多个较小的数组】，每个数组向量的大小为size
%    第n个数组的中心元素都为第n个元素，size必须为单数
%     Examples: 0-10, size is 3.


if nargin < 3
    wrap = 'false';
end

if strcmp(wrap,'true')
    dt = size - floor(size/2);
    shape=[length(data),size];
    %水平堆叠数组
    data=horzcat(data(end-dt+1:end),data,data(1:dt));
    %vertcat 垂直堆叠数组
elseif strcmp(wrap,'false')
    shape=[length(data)-size+1,size];
end
out=ones(shape);

for i=1:length(out)
    out(i,:)=data(1:size);
    data = circshift(data, [0, -1]);
end
% 将第一行移动至最后一行
data_out=[out(2:end,:);out(1,:)];
end

