%ER 该代码没有问题，注意symbTx的生成即可，需要的是bit数
function   rexdB=Calc_ER(E,M,symbTx,SpS)
% Calculate extinction ratio and intensity levels position
% Note:symbTx must be noted，it should be the symbols not be modulator。
% for example ： OOK should be the 0,1 。PAM4 should be the 0,1,2,3
P = abs(E).^2;
%down
Psamp = P(1:SpS:end);
%每个电平对应一个功率
Pl = zeros(1, M);
for k = 1:M
    Pl(k) = mean(Psamp(symbTx == k-1));
end
%升序排列
Pl = sort(Pl, 'ascend');
%ER
rexdB = 10*log10(min(Pl)/max(Pl));
fprintf('Estimated extinction ratio = %.2f dB\n', rexdB)
%归一化
Pl = Pl - Pl(1);
Pl = (M-1)*Pl/Pl(end);
fprintf('Normalized optical levels: %s\n', sprintf('%.3f\t', Pl));
end