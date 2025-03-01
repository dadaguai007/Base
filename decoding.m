function[y] = decoding(sig,chose)
switch lower(chose)
    case 'binary'
        y1 = ((sig)+3)/2;
        y2 = dec2bin(y1);
        y3 = boolean(y2-'0');
        y = reshape(y3',[1,2*length(y3)]);
        y = y(:);
    case  'gray'
%                 y1 = ((sig)+3)/2;
%                 y2 = dec2bin(y1);
%                 y3 = boolean(y2-'0');	%将字符转成逻辑量
%                 for i =1:length(y3)
%                     a(i,:) = y3(i,1);  %首位保留
%                     b(i,:) = xor(y3(i,1),y3(i,2));  %其他位与上一位异或
%                 end
%                 c_1 = num2str(a);
%                 c_2 = num2str(b);
%                 y4 = [c_1,c_2];
%                 y5 = boolean(y4-'0');
%                 y = reshape(y5',[1,2*length(y5)]);
%                 y = y(:);

        %%PAM4
        for i = 1:length(sig)
            if sig(i) == -3;
                y(i*2-1) = 0;
                y(i*2) = 0;
            elseif sig(i) == -1;
                y(i*2-1) = 1;
                y(i*2) = 0;
                %原始对应的是01
            elseif sig(i) == 1;
                y(i*2-1) = 1;
                y(i*2) = 1;
            else sig(i) == 3;
                y(i*2-1) = 0;
                y(i*2) = 1;
                %原始对应的是10
            end
        end
        y = y(:);


end
