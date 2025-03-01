function Output=MLSE_viterbi2(input,tap)
symbol=[-1 1];%h.Constellation;
M=length(symbol);
L=length(input);
min_length=zeros(M^2,L);
min_sym=zeros(M^2,L);
for ii=1:L
    if ii==1
        for jj=1:M
            d=symbol(jj)-input(1);
            min_length(jj,1)=d*conj(d);
        end
        continue;
    end
    if ii==2
        for jj=1:M
            temp=0;
            flag=0;
            for kk=1:M
                d=symbol(jj)+tap(1)*symbol(kk)-input(2);
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
            min_length(jj:M:end,ii)=temp*ones(M,1);
            min_sym(jj:M:end,ii)=flag*ones(M,1);
        end
        continue;
    end
    
    for jj=1:M
        for kk=1:M
            temp=0;
            flag=0;
            for ll=1:M
                d=symbol(jj)+tap(1)*symbol(kk)+tap(2)*symbol(ll)-input(ii);
                e=d*conj(d)+min_length(kk+M*(ll-1),ii-1);
                if flag==0
                    temp=e;
                    flag=ll;
                else
                    if e<temp
                        temp=e;
                        flag=ll;
                    end
                end                   
            end
            min_length(jj+M*(kk-1),ii)=temp;
            min_sym(jj+M*(kk-1),ii)=flag;
        end
    end
end
Output=zeros(1,L);
flag=0;
node1=0;
node2=0;
for ii=L:-1:1
    if ii==L
        [~,flag]=min(min_length(:,ii));
        node1=mod(flag,M);
        if node1==0
            node1=M;
        end
    elseif ii==L-1
        node2=ceil(flag/M);        
    else
        if mod(L-ii,2)==0
            node1=min_sym(node1+(node2-1)*M,ii+2);
        else
            node2=min_sym((node1-1)*M+node2,ii+2);
        end
    end
    if mod(L-ii,2)==0
        Output(ii)=symbol(node1);
    else
        Output(ii)=symbol(node2);
    end
end
end