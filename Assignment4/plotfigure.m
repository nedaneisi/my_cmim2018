
function plotfigure(T,Analytical_Solution,Xn,Tit)
Figure_setup; hold on
plot(T,Analytical_Solution)
plot(T,Xn,'r*')
% plot(T,Xn,'r--')
grid on
legend('analytical',Tit)
% set(lgd,'Location','BestOutside')
xlabel('Times ,s')
ylabel('Displacement ,m')
% title(Tit)
fname = 'C:\Data\neda\my matalb reaserch\computional mech_assingment\Assingment Neda\Assignment4\Figures';
filename=Tit
saveas(gca, fullfile(fname, filename), 'emf');
saveas(gca, fullfile(fname, filename), 'fig');

end