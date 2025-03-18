function codedData = Upper_Lower_code_pam(data)
    if mod(length(data), 2) ~= 0
        error('The length of the input data must be even for Alamouti coding.');
    end
    data = repelem(data, 2);
    codedData = zeros(2, size(data, 2) );
    for i = 1:size(data, 2)
        s1 = data(i);
        codedData(1, i) = s1;        % first time slot
        codedData(2, i) = (-1).^i*s1;  % second time slot
    end
end
