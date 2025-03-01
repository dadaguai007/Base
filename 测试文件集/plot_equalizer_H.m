% plot equalizer H

% 绘图均衡前后的频率响应
% W is the coefficient
% time domain
eq.h = W; % [Ntaps x Nfilters]
% eq.Wdc = Wdc; % ajustable dc bias
% Frequency response is calculated for only one of the filters
% frequency domain
eq.Hff = @(f) freqz(eq.h(:, 1), 1, 2*pi*f).*exp(1j*2*pi*f*grpdelay(eq.h(:, 1), 1, 1)); % removed group delay
eq.MSE = abs(e).^2;

% Plot MSE
figure(400), clf
subplot(221), hold on, box on
MSE = abs(e).^2;
plot(MSE);
plot(smooth(MSE, 201), 'r', 'Linewidth', 2);
a = axis;
plot(eq.Ndiscard(1)*[1 1], a(3:4), ':k')
plot((Nsymb-eq.Ndiscard(2))*[1 1], a(3:4), ':k')
legend('MSE', 'Smooth MSE', 'Valid window')
xlabel('Iteration')
ylabel('MSE')
title('Equalizer MSE')


fnorm = linspace(-0.5, 0.5);
for k = 1:size(W, 2)
    subplot(222), hold on, box on
    stem(-(eq.Ntaps-1)/2:(eq.Ntaps-1)/2, abs(W(:, k)))
    %     plot([-(eq.Ntaps-1)/2 (eq.Ntaps-1)/2], Wdc(k)*[1 1])

    % frequency domain
    Hw = freqz(W(:, k), 1, 2*pi*fnorm);
    subplot(223), hold on, box on
    plot(Rs*eq.ros*fnorm/1e9, abs(Hw).^2)

    % phase
    subplot(224), hold on, box on
    plot(Rs*eq.ros*fnorm/1e9, unwrap(angle(Hw)))
end
subplot(222)
xlabel('n')
ylabel('|W(n)|')
title('Equalizer taps')
legend('Equalizer taps', 'DC bias')

subplot(223)
a = axis;
axis([-Rs*eq.ros/2e9 Rs*eq.ros/2e9 a(3:4)]);
xlabel('Frequency (GHz)')
ylabel('|H_W(f)|^2')

subplot(224)
a = axis;
axis([-Rs*eq.ros/2e9 Rs*eq.ros/2e9 a(3:4)]);
xlabel('Frequency (GHz)')
ylabel('arg(H_W(f))')
drawnow
