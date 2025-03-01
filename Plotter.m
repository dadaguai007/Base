%% figure settings
function Plotter(titlename,xname,yname)
% legend('FDE CP=1','FDE CP=3','FDE CP=5','FDE CP=10', 'TDE Taps=10');
xlabel(xname);
ylabel(yname);
title(titlename);
grid on;
axis square;
set(gca,'FontSize',14);
ylim tight
end

