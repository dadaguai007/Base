function Output=MLSE_viterbi(input,tap)
% input=[1,3,2.5,2,1,8];
% tap=0.1;
symbol=[-3 -1 1 3];%h.Constellation;
M=length(symbol);
L=length(input);
min_length=zeros(M,L);

min_sym=zeros(M,L);
for ii=1:L
    if ii==1
        for jj=1:M
            %参考与输入之间的差值
            d=symbol(jj)-input(1);
            %模平方
            min_length(jj,1)=d*conj(d);
        end
        continue;
    end
    for jj=1:M
        temp=0;
        flag=0;
        for kk=1:M
            d=symbol(jj)+tap*symbol(kk)-input(ii);
%             d=symbol(jj)-(input(ii)-tap*input(ii-1));
            e=d*conj(d)+min_length(kk,ii-1);
            if flag==0
                temp=e;
                flag=kk;
            else
                if e<temp
                    temp=e;
                    flag=kk;
                end
            end
        end
        min_length(jj,ii)=temp;
        min_sym(jj,ii)=flag;
    end
end
Output=zeros(1,L);
flag=0;
for ii=L:-1:1
    if flag==0
        [a,flag]=min(min_length(:,ii));
    else
        flag=min_sym(flag,ii+1);
    end
    Output(ii)=symbol(flag);
end
end