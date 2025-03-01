
function out = rolling_window(data, size, wrap)
% it is uesd for the LUT
%    将一个一维数组重新排列成一个有重叠帧的二维数组【多个较小的数组】，每个数组向量的大小为size
%     Examples: 0-10, size is 3.
%     >>> rolling_window(0:10, 3)
%     ([[0, 1, 2],
%       [1, 2, 3],
%       [2, 3, 4],
%       [3, 4, 5],
%       [4, 5, 6],
%       [5, 6, 7],
%       [6, 7, 8],
%       [7, 8, 9]])

if nargin < 3
    wrap = 'false';
end

if strcmp(wrap,'true')
    dt = size - 1;
    shape=[length(data),size];
    %水平堆叠数组
    data=horzcat(data,data(1:dt));
    %vertcat 垂直堆叠数组
elseif strcmp(wrap,'false')
    shape=[length(data)-size+1,size];
end
out=ones(shape);

for i=1:length(out)
    out(i,:)=data(1:size);
    data = circshift(data, [0, -1]);
end

end

