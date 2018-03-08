Figure_setup; hold on
plot(T,yn(1,:))
plot(T,Xn,'r*') % position- SIE and analyitical
plot(T,Analytical_Solution) % position- SIE and analyitical
legend('forward euler','SIE','Analytical solution')
% title('forward euler SIE Analytical solution')
fname = 'C:\Data\neda\my matalb reaserch\computional mech_assingment\Assingment Neda\Assignment4\Figures';
filename='forwardeuler_SIE_Analytical '
saveas(gca, fullfile(fname, filename), 'emf');

% Total energy Forward Euler
potential_energy_feuler=0.5*k.*yn(1,:).^2;
kinetic_energy_feuler=0.5*m.*yn(2,:).^2;
Figure_setup; hold on
plot(T,kinetic_energy_feuler)
plot(T,potential_energy_feuler)
plot(T,kinetic_energy_feuler+potential_energy_feuler)
xlabel('Times ,s')
ylabel('Energy')
grid on
legend('potential energy','Kinetic energy','total energy')
% title('Energy of system, Forward Euler ')
filename='Energy_Forward Euler' 
saveas(gca, fullfile(fname, filename), 'emf');


% % Total energy SIE 
Figure_setup; hold on
plot(T,kinetic_energy)
plot(T,potential_energy)
plot(T,kinetic_energy+potential_energy)
xlabel('Times ,s')
ylabel('Energy')
grid on
legend('potential energy','Kinetic energy','total energy')
% title('Energy of system, SIE ')
filename='Energy_SIE' 
saveas(gca, fullfile(fname, filename), 'emf');



Figure_setup; hold on
plot(T,Analytical_Solution,'--.')
plot(T,Xn)
plot(T,Xn_RK4(1,:))
legend('analytical solution','Euler expilicit', 'Runge Kutta method')
grid on
xlabel('Times ,s')
ylabel('Displacement ,m')
% title('Position for SIF,RK4 and analytical method')
filename='PositionSIF_RK4_analytical' 
saveas(gca, fullfile(fname, filename), 'emf');
