
function plotfigure3mass(T,X,Tit)
Figure_setup; hold on
plot(T,X(1,:))
plot(T,X(2,:))
plot(T,X(3,:))
legend('m1','m2', 'm3')
grid on
xlabel('Times ,s')
ylabel('Displacement ,m')
% title(Tit)
fname = 'C:\Data\neda\my matalb reaserch\computional mech_assingment\Assingment Neda\Assignment4\Figures';
filename=Tit
saveas(gca, fullfile(fname, filename), 'emf');
saveas(gca, fullfile(fname, filename), 'fig');
end