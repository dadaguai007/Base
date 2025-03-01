function [E_rec, phase] = phase_partition_16qam(E, SpS)
    % 16-QAM blind phase recovery using QPSK partitioning.

    % A blind phase estimator for 16-QAM signals based on partitioning the signal
    % into 3 rings, which are then phase estimated using traditional V-V phase
    % estimation after Fatadin et al [1].


    % References
    % ----------
    % [1] I. Fatadin, D. Ives, and S. Savory, “Laser linewidth tolerance
    %    for 16-QAM coherent optical systems using QPSK partitioning,”
    %    Photonics Technol. Lett. IEEE, vol. 22, no. 9, pp. 631–633, May 2010.
    
    % E should be 2*N  （两个模式）
    E2d = E; % Ensure E is 2D
    dphi = pi / 4 + atan(1 / 3);
    modes = size(E2d, 1);
    ph_out = [];
    L=size(E2d,1);
    for j = 1:modes
        % Partition QPSK signal into qpsk constellation and non-qpsk const
        [c1_m, c2_m] = partition_16qam(E2d(j, :));
        
        Sx = zeros(L,1);
        Sx(c2_m) = (E2d(j, c2_m) .* exp(1i * dphi)).^4;
        
        So = zeros(L,1);
        So(c2_m) = (E2d(j, c2_m) .* exp(1i * -dphi)).^4;
        
        S1 = zeros(L,1);
        S1(c1_m) = (E2d(j, c1_m)).^4;
        

        phi_est = zeros(L,1);
        
        for i = 1:SpS:numel(E2d(j, :))
            S1_sum = sum(S1(i:i + SpS - 1));
            Sx_tmp = min([S1_sum - Sx(i:i + SpS - 1), S1_sum - So(i:i + SpS - 1)], [], 2);
            phi_est(i:i + SpS - 1) = angle(S1_sum + sum(Sx_tmp));
        end
        
        ph_out = [ph_out, unwrap(phi_est) / 4 - pi / 4];
    end
%     ph_out 是一个矩阵，这里或许会有大的bug
    phi_out = reshape(ph_out, length(E2d),[]);
    
    if size(E, 1) == 1
        E_rec = (E2d .* exp(-1i * phi_est)).';
        phase = phi_out.';
    else
        E_rec = E2d .* exp(-1i * phi_est);
        phase = phi_out;
    end
end
