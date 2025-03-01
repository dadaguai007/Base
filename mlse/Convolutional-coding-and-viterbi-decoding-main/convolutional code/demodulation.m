function [data_out] = demodulation(data_in,mod)
    switch(mod)
        case 'BPSK'
            L = length(data_in);
            data_out = zeros(1,L);

        for k=1:L
            if (data_in(k)>0)
                data_out(k)= 0;
            else
                data_out(k)= 1;

            end
        end

        case 'QPSK'
            data_inI = (1 - sign(real(data_in))) / 2;
            data_inQ = (1 - sign(imag(data_in))) / 2;
            L=length(data_inI);
            data_out=zeros(1,2*L);
            for k=1:L
                data_out(2*k-1:2*k)=[data_inI(k) data_inQ(k)];
                data_out(1,2*k-1)=data_inI(k);
                data_out(1,2*k)=data_inQ(k);
            end
 
        case 'AMPM'
            const = [1 + -1i,-3 + 3i,1 +  3i,-3 - 1i,3 - 3i,-1 + 1i,+3 + 1i,-1 - 3i]/sqrt(10);
            L = length(data_in);
            y_term = repmat(data_in,8,1);
            const_term = repmat(transpose(const),1,L);
            dis = abs(y_term - const_term);
            [~,index] = min(dis);
            bits_hat = de2bi((index-1)','left-msb');
            data_out = reshape(bits_hat',1,3*L);
            data_out = data_out(1:1e5);

    end

end