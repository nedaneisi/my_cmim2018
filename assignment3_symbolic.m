
clc
close all
clear all

syms  f0 omega_n omega  zeta gamma   A_s B_s real
assumeAlso(zeta > 0);
assumeAlso(omega_n > 0);
assumeAlso(omega > 0);


syms x_p(t)  x(t) X(t) theta(t)

x_p=A_s*cos(omega*t)+B_s*sin(omega*t);
Dx = diff(x_p,t);
DDx = diff(Dx,t);

eq=diff(x, t,2)+2*zeta*omega_n*diff(x, t)+omega_n^2*x-f0*cos(omega*t);
eq=expand(subs(subs(subs(eq,diff(x, t,2),DDx),diff(x, t),Dx),x,x_p));

% collect(eq1, 'sin' )

eq=collect(eq, {'sin' 'cos'});

eq1=subs(eq,sin(omega*t),0);
eq2=subs(eq,cos(omega*t),0);

eq1=subs(eq1,cos(omega*t),1);
eq2=subs(eq2,sin(omega*t),1);
% 
S=solve([eq1==0,eq2==0], [A_s B_s]);

S.A_s;
S.B_s;

X=simplify(sqrt(S.A_s^2+S.B_s^2));
theta=atan2(S.B_s,S.A_s);
disp(X)
disp(theta)

%% plot result
k=1;m=1; OMEGA_n=sqrt(k/m); % naturl freq for the input of the probllem 
F=1;

zeta0=[0.00001 0.1 0.2 0.3 0.5 1];% damping ratio
gamma=(0:0.01:3); % freq. ratio
omega0=gamma*OMEGA_n;

 
%% plot result 
f_width=450;
f_height=300;
figure ('Color','w','Position',[200,300,f_width,f_height]); 

set(gcf,'color','w');
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 12)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 12)
set(groot,'defaultLineLineWidth',2)


figure;% amplitude
hold on
for i=1:length(zeta0)
 plot(gamma, double(subs(subs(subs(subs(X,omega_n, OMEGA_n),zeta,zeta0(i)),f0,F),[omega], {omega0} )),'k')
  Legend{i}=strcat('{\zeta=}', num2str(zeta0(i)));
end
xlabel('Frequency ratio {\gamma}={\omega}/{\omega}_n')
ylabel('Amplitude X k/F_0')
 ylim([0 5])
% legend(Legend);
grid on
% title('steady state responce due to unbalance mass')
fname = '\\maa1.cc.lut.fi\home\h2004\Desktop\Assignment 3';
filename=' vibration 1dof harmonic'
saveas(gca, fullfile(fname, filename), 'emf');
saveas(gca, fullfile(fname, filename), 'fig');



figure;% phase
hold on
for i=1:length(zeta0)
plot(gamma, (180/pi)*double(subs(subs(subs(subs(theta,omega_n, OMEGA_n),zeta,zeta0(i)),f0,F),[omega], {omega0} )),'k')
Legend{i}=strcat('{\zeta=}', num2str(zeta0(i)));
end
xlabel('Frequency ratio {\gamma}={\omega}/{\omega}_n')
ylabel('phase angle (degree)  {\theta} ')
legend(Legend);
ylim([0 180])
grid on
title('bode plot for various damping level')
fname = '\\maa1.cc.lut.fi\home\h2004\Desktop\Assignment 3';
filename=' phase angle 1dof harmonic'
saveas(gca, fullfile(fname, filename), 'emf');
saveas(gca, fullfile(fname, filename), 'fig');



