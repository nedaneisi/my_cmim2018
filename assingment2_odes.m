%% Computional mechanic course_assignment2

clear all;close all;clc
%% Direct solution
% direct solution is based on state space method. it is required to form the state space matrix, where A=[ x:xdot]
% Input data
M=[10]; % mass matrix Kg
K=[10];%Stiffness matrix N/m
% zeta=0.1;
% C=2*zeta*sqrt(K*M);
 C=[1]; % Damping matrix N/mm
omega_f=0:0.1:10; % range of excitation frequencies
F0=2; % external force , in the state space method the external force will be divied into sin and cos components.
%M*X"+C*X'+KX=F0*sin(omega_f*t)+F0*cos(omega_f*t)
F=[F0 0 ]'; %External force vector [x1,x2,x1dot,x2dot]
t=0:0.01:10;% simulation time


% initial conditions for numerical time intergarion for specific naturalfreq
x0=0.02;
xdot0=0;
y0=[x0(1) xdot0(1)]';
omega_f0=0.6;% natural frequency of the exictation force


%To find a noninteger value, use a tolerance value based on your data.
index = find(abs(omega_f-omega_f0) < 0.0001);% in the last section the amplitude X is obtained for differnet range of omega_f. 



% eigen value natural frequency of undamped system
[FiiN d]=eig(K,M); %FiiN=eigenvectors, d=eigenvalues
omega_n=sqrt(d) % natural freq rad/s
freq=omega_n/(2*pi) %natural freq rad/s Hz
for i=1:length(omega_f)
A=[K-omega_f(i)^2*M -omega_f(i)*C
omega_f(i)*C K-omega_f(i)^2*M];
ab(i,:)=A^(-1)*F; % ab is the matrix for the coefficient of direct solution which is describe in the pdf file attached to this assingment
end

X_direct(:,:)=sqrt(ab(:,1).^2+ab(:,2).^2);% amplitudes
phi=atan2(ab(:,2),ab(:,1))'; % atan2(y,x)  the angle phi is phase angle

% in next section the amplitude is compared for sample frequeny of
% omega_f0=5. Therefore the index is used to see it is for which array of matrix phi and X__direct
x_times_direct=ab(index,1)*sin(omega_f0*t)+ab(index,2)*cos(omega_f0*t);


zeta=C/(2*M*omega_n);
omega_d=omega_n*sqrt(1-zeta^2);
XX=(F0/M)/sqrt(  (omega_n^2-omega_f(index)^2)^2  +  (2*zeta*omega_n*omega_f(index))^2);
teta=atan2(omega_n^2-omega_f(index)^2,2*zeta*omega_n*omega_f(index));
phiphi=atan2(  xdot0+(x0-XX*cos(teta))*zeta*omega_n-omega_f(index)*XX*sin(teta), omega_d*(x0-XX*cos(teta)));
CC=(x0-XX*cos(teta))/sin(phiphi);
x_time_harmonic=CC*exp(-zeta*omega_n*t).*sin(omega_d*t+phiphi)+XX*cos(omega_f(index)*t-teta);

%% numerical time integration


err=1e-4;
options = odeset('RelTol',err,'MaxStep',1e-4,'Stats','off');
[t,y_ode45]=ode45(@eom,t,y0,options,M,C,K,F0,omega_f0);% ode45 Mtlab solver for time integration
[t,y_ode15s]=ode15s(@eom,t,y0,options,M,C,K,F0,omega_f0);% ode45 Mtlab solver for time integration
[t,y_ode23]=ode23(@eom,t,y0,options,M,C,K,F0,omega_f0);% ode45 Mtlab solver for time integration
[t,y_ode113]=ode113(@eom,t,y0,options,M,C,K,F0,omega_f0);% ode45 Mtlab solver for time integration


	
%% exact solution by dsolve Matlab
syms y(tt)
ode = M*diff(y,tt,2)+C*diff(y,tt)+K*y-F0*sin(omega_f0*tt) == 0;
Dy = diff(y);
cond1 = y(0) == x0;
cond2 = Dy(0) == xdot0;
cond = [cond1 cond2 ];
ySol(tt) = dsolve(ode,cond);
tt=t;
for i=1:length(t)
Y_exact_matlab(i)=double(subs(ySol,t(i)));
end
figure
 plot(Y_exact_matlab)
%% plot results
f_width=450;
f_height=350;
set(gcf,'color','w');
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 10)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 10)
set(groot,'defaultLineLineWidth',1)
figure ('Color','w','Position',[200,300,f_width,f_height]);
subplot(2,1,1)
plot(omega_f,X_direct(:))
grid on
title('Vibrtion amplitude-Direct solution')
xlabel('Excitation frequency [rad/s]')
ylabel('Displacement [m]')
subplot(2,1,2)
plot(omega_f,phi(:))
grid on
title('phase-direct solution')
xlabel('Excitation frequency [rad/s]')
ylabel('phase [rad]')


figure ('Color','w','Position',[200,300,f_width,f_height]); hold on
plot(t,x_times_direct,'k')
plot(t,x_time_harmonic,'r')
grid on
title('Vibrtion amplitude direct solution and harmonic solution')
xlabel('Time [s]')
ylabel('Displacement [m]')




figure ('Color','w','Position',[200,300,f_width,f_height]); hold on
plot(t,x_times_direct)
plot(t,x_time_harmonic)
plot(t,y_ode45(:,1))
plot(t,y_ode15s(:,1))
plot(t,y_ode23(:,1))
plot(t,y_ode113(:,1))
grid on
title('Vibrtion amplitude for various numerical integrations')
xlabel('Time [s]')
ylabel('Displacement [m]')
legend('direct solution','harmonic','ode 45','ode15s','ode23','ode113')



figure ('Color','w','Position',[200,300,f_width,f_height]); hold on
plot(t,y_ode45(:,1)-Y_exact_matlab','r')
plot(t,y_ode15s(:,1)-Y_exact_matlab','y')
plot(t,y_ode23(:,1)-Y_exact_matlab','g')
plot(t,y_ode113(:,1)-Y_exact_matlab','c')
grid on
title('error for various numerical integrations')
xlabel('Time [s]')
ylabel('Displacement [m]')
legend('ode 45','ode15s','ode23','ode113')



% time integartion part
function xdot=eom(t,y0,M,C,K,F0,omega_f0)
Ft=F0*sin(omega_f0*t);
AA=[ 0 eye(size(M));-inv(M)*K -inv(M)*C];
xdot=AA*y0+[0 ;inv(M)*Ft];
end