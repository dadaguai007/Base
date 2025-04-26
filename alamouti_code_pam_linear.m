% 复现：A better than alamouti OSTBC for MIMO backscatter communications

function codedData = alamouti_code_pam_linear(data)
    if mod(length(data), 2) ~= 0
        error('The length of the input data must be even for Alamouti coding.');
    end
    reverse_martix=[1,1;-1,1];
    data = reshape(data, 2, []);
    codedData = zeros(2, size(data, 2) * 2);
    for i = 1:size(data, 2)
        s1 = data(1, i);
        s2 = data(2, i);

        codedData(:, 2*i-1) = [s1; s2];        % first time slot
        codedData(:, 2*i) = [-(s2); (s1)];  % second time slot
        codedData(:, 2*i-1:2*i)=reverse_martix* codedData(:, 2*i-1:2*i);
    end
end
