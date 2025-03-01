function [data_out] = modulation(data_in,mod)
    switch(mod)
        case 'BPSK'
            data_out = zeros(1,length(data_in));
            for k=1:length(data_in)
                if (data_in(k)==1)
                    data_out(k)= -1;
                else
                    data_out(k)= 1;
                end

            end

        case 'QPSK'
            Kmod=1/sqrt(2);
            L = length(data_in)/2;
            data_outI = zeros(1,L);
            data_outQ = zeros(1,L);

            for k=1:L
                data_outI(k)=Kmod * (1 - 2*data_in(2*k-1));
                data_outQ(k)=Kmod * (1 - 2*data_in(2*k));
            end
            data_out = data_outI+1i*data_outQ; 
 
        case 'AMPM'
            const = [1 + -1i,-3 + 3i,1 +  3i,-3 - 1i,3 - 3i,-1 + 1i,+3 + 1i,-1 - 3i]/sqrt(10);
            index = bi2de((buffer(data_in,3))','left-msb');
            data_out = const(index+1);


    end
end