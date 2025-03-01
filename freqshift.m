function xshift = freqshift(x,t,fshift)
% 频率上的移位，相当于相位调制
if size(x,2)>size(x,1)
    x=x';
    transpose = true;
end
% xshift = x.*repmat(exp(1j*2*pi*fshift*t),size(x,1),1);
% 进行相位调制
xshift = x.*exp(1j*2*pi*fshift*t).';
if transpose
    xshift = xshift.';
end

end
