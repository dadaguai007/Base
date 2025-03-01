%% pnle3_lms
function pnle3 = PNLE3_LMS(S, sig, train)
sig = sig(:);
step = S.step;
t_ff = S.t_ff1;
t_fb1 = S.t_fb1;
t_fb2 = S.t_fb2;
c_fe1 = zeros(t_ff,1);
c_fe2 = zeros(t_fb1,1);
c_fe3 = zeros(t_fb2,1);
x_fe1 = zeros(t_ff,1);

for i = 1:N/sps-1
    % linear data
    x_fe1 = cat(1,x_fe1(sps+1:end),sig(sps*i-sps+1:1:sps*i));
    if t_ff-sps*i <= 0
        % nonlinear data
        x_fe2 = x_fe1(1:t_fb1).^2;
        x_fe3 = x_fe1(1:t_fb2).^3;
        % linear data * linear coefficients + linear data * linear coefficients
        y(i,:) = x_fe1.'*c_fe1 + x_fe2.'*c_fe2 + x_fe3.'*c_fe3;
        % update coefficients
        if i>1
            e(i,:) = train(i-1)-y(i);
            c_fe1 = c_fe1+step*e(i)*x_fe1;
            c_fe2 = c_fe2+step*e(i)*x_fe2;
            c_fe3 = c_fe3+step*e(i)*x_fe3;
        end
    end
end
e = e(:);
y = y(:);
pnle3=y;