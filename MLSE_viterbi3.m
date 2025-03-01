function Output=MLSE_viterbi3(input,tap,h,delay)
symbol=h.Constellation;
M=length(symbol);
L=length(input);
min_length=zeros(M^3,L);
min_sym=zeros(M^3,L);
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
            min_length(jj:M:end,ii)=temp*ones(M^2,1);
            min_sym(jj:M:end,ii)=flag*ones(M^2,1);
        end
        continue;
    end
    if ii==3
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
                min_length(jj+M*(kk-1):M^2:end,ii)=temp*ones(M,1);
                min_sym(jj+M*(kk-1):M^2:end,ii)=flag*ones(M,1);
            end
        end
    end
    
    for jj=1:M
        for kk=1:M
            for ll=1:M
                temp=0;
                flag=0;
                for mm=1:M
                    d=symbol(jj)+tap(1)*symbol(kk)+tap(2)*symbol(ll)+tap(3)*symbol(mm)-input(ii);
                    e=d*conj(d)+min_length(kk+M*(ll-1)+M^2*(mm-1),ii-1);
                    if flag==0
                        temp=e;
                        flag=mm;
                    else
                        if e<temp
                            temp=e;
                            flag=mm;
                        end
                    end                   
                end
                min_length(jj+M*(kk-1)+M^2*(ll-1),ii)=temp;
                min_sym(jj+M*(kk-1)+M^2*(ll-1),ii)=flag;
            end
        end
    end
    
end
Output=zeros(1,L);
flag=0;
node1=0;
node2=0;
node3=0;
for ii=L:-1:1
    if ii==L
        [~,flag]=min(min_length(:,ii));
        node1=mod(flag,M);
        if node1==0
            node1=M;
        end
    elseif ii==L-1
        node2=mod(flag,M^2);
        if node2==0
            node2=M^2;
        end
        node2=ceil(node2/M);
    elseif ii==L-2
        node3=ceil(flag/M^2);
    else
        if mod(L-ii,3)==0
            node1=min_sym(node1+(node2-1)*M+(node3-1)*M^2,ii+3);
        elseif mod(L-ii,3)==1
            node2=min_sym(node2+(node3-1)*M+(node1-1)*M^2,ii+3);
        else
            node3=min_sym(node3+(node1-1)*M+(node2-1)*M^2,ii+3);
        end
    end
    if mod(L-ii,3)==0
        Output(ii)=symbol(node1);
    elseif mod(L-ii,3)==1
        Output(ii)=symbol(node2);
    else
        Output(ii)=symbol(node3);
    end
end
end