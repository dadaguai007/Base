function Plot_error_weight(w,e,squaredError)
figure
stem(w)
title("均衡器抽头")
figure
plot(abs(e))
xlabel("迭代次数")
ylabel("误差")

% % 画图
% %抽头对比图
% figure
% stem(w_lms)
% hold on
% stem(w_rls)
% legend("FFE-LMS","FFE-RLS")
% title("均衡器抽头")
%
% %训练的对比图
% figure
% plot(abs(e_lms))
% hold on
% plot(abs(e_rls))
% legend("FFE-LMS","FFE-RLS")
% xlabel("迭代次数")
% ylabel("误差")



figure;
semilogy(squaredError)
xlabel("迭代次数")
ylabel("误差")

L = 100;
MA_SE = conv(squaredError,ones(1, L)/L, 'same');
figure;
semilogy(MA_SE)
xlabel("迭代次数")
ylabel("误差")

end
