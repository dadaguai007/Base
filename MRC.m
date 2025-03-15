function [ber_simulation] = MRC(SNR, A, N, NRx)

% Rx_MRC(SNR3, A, N, 2); 示例
    absolute_amplitude = sqrt(2) * A;
    M = 2; % BPSK
    
    % Energy Definitions
    Eb = (absolute_amplitude^2) / 2; 
    Es = Eb;

    % Voltage of the sent signals
    voltages = zeros(1,M);
    for i = 1:M
        theta = (i-1) *  180;
      
        voltages(i) = -Es * cosd(theta);
       
    end

    % Generate Random Symbols
    tx = randsrc(1, N, voltages);

    % Linear SNR
    linear_SNR = 10.^(SNR / 10);
    
    % Analytical and simulated error rate arrays
    ber_simulation = zeros(1, length(SNR));

    % Fading, Energy and rx Arrays
    h = zeros(NRx, N);
    energy_of_received_symbol1 = zeros(NRx);
    N0 = zeros(NRx);
    rx = zeros(NRx, N);
    decision = zeros(1,N);

    for i = 1:length(SNR)
        for j = 1: NRx
            % Rayleigh Fading Channels
            h(j,:) = (1 / sqrt(2)) .* (randn(1,N) + 1i * randn(1,N));

            % Power Calculations
            energy_of_received_symbol1(j) = mean(abs(h(j,:).^2)) * Es;
            N0(j) = energy_of_received_symbol1(j) / linear_SNR(i);

            rx(j,:) = h(j,:) .* tx + sqrt(N0(1) / 2) * (randn(1, N) + 1i * randn(1, N));
        end

        % Decision Circuit
        rx_sum = sum(conj(h) .* rx, 1);
        decision(rx_sum > 0) = voltages(2);
        decision(rx_sum < 0) = voltages(1);
        
        % Ber Calculation
        ber_simulation(i) = sum(decision ~= tx)/ N;
    end
end