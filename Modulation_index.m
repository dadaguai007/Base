function m=Modulation_index(E,Vpi,type)
% 计算调制器的调制深度
if strcmp(type,'ofdm')
    Vpp=2*3*std(real(E));

    m=Vpp/Vpi;

else
    Vpp=max(E)-min(E);

    m=Vpp/Vpi;
end

end