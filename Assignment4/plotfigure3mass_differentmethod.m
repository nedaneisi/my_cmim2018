
function plotfigure3mass_differentmethod(T,Xn_3MASS,Xn_RK43MASS,X,y_45,y_15s,k,Tit)
Figure_setup; hold on
% k=1 shows mass 1 k=2 mass2 k3=mass 3
plot(T,Xn_3MASS(k,:))
plot(T,Xn_RK43MASS(k,:),'b--')
plot(T,X(k,:),'r:')
plot(T,y_45(:,k)')
plot(T,y_15s(:,k)')

legend('SIE','RK4', 'Modal solution', 'ode 45', 'ode 15s')
grid on
xlabel('Times ,s')
ylabel('Displacement ,m')
fname = 'C:\Data\neda\my matalb reaserch\computional mech_assingment\Assingment Neda\Assignment4\Figures';
filename=Tit
saveas(gca, fullfile(fname, filename), 'emf');
saveas(gca, fullfile(fname, filename), 'fig');

end