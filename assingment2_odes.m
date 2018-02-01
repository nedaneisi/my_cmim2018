%% Computional mechanic course_assignment2

clear all;close all;clc
%% Direct solution
% direct solution is based on state space method. it is required to form
the state space matrix, where A=[ x:xdot]
% Input data
M=[10]; % mass matrix Kg
K=[10];%Stiffness matrix N/m
C=[1]; % Damping matrix N/mm
omega_f=0:0.1:10; % range of excitation frequencies
F0=2; % external force , in the state space method the external force will be
divied into sin and cos components.
%M*X"+C*X'+KX=F0*sin(omega_f*t)+F0*cos(omega_f*t)
F=[F0 0 ]'; %External force vector [x1,x2,x1dot,x2dot]
t=0:0.01:10;% simulation time
% eigen value natural frequency of undamped system
[FiiN d]=eig(K,M); %FiiN=eigenvectors, d=eigenvalues
omega_n=sqrt(d) % natural freq rad/s
freq=omega_n/(2*pi) %natural freq rad/s Hz
for i=1:length(omega_f)
A=[K-omega_f(i)^2*M -omega_f(i)*C
omega_f(i)*C K-omega_f(i)^2*M];
ab(i,:)=A^(-1)*F; % ab is the matrix for the coefficient of direct
solution which is describe in the pdf file attached to this assingment
end
size(ab)
X_direct(:,:)=sqrt(ab(:,1).^2+ab(:,2).^2);% amplitudes
X_direct=sqrt(a_i^2+b_i^2)
%% numerical time integration
% initial conditions for numerical time intergarion for specific naturalfreq
x0=0;
xdot0=0;
y0=[x0(1) xdot0(1)]';
omega_f0=5;% natural frequency of the exictation force
err=1e-1;
options = odeset('RelTol',err,'MaxStep',1e-3,'Stats','off');
[t,y_ode45]=ode45(@eom,t,y0,options,M,C,K,F0,omega_f0);% ode45 Mtlab solver
for time integration
[t,y_ode15s]=ode15s(@eom,t,y0,options,M,C,K,F0,omega_f0);% ode45 Mtlab solver
for time integration
[t,y_ode23]=ode23(@eom,t,y0,options,M,C,K,F0,omega_f0);% ode45 Mtlab solver
for time integration
[t,y_ode113]=ode113(@eom,t,y0,options,M,C,K,F0,omega_f0);% ode45 Mtlab solver
for time integration
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
plot(omega_f,X_direct(:))
grid on
title('Vibrtion amplitude-Direct solution')
xlabel('Excitation frequency [rad/s]')
ylabel('Displacement [m]')
%
figure ('Color','w','Position',[200,300,f_width,f_height]); hold on
plot(t,y_ode45(:,1),'--k')
plot(t,y_ode15s(:,1),'r')
plot(t,y_ode23(:,1),'g')
plot(t,y_ode113(:,1),'c')
grid on
title('Vibrtion amplitude for various numerical integrations')
xlabel('Time [s]')
ylabel('Displacement [m]')
legend('ode 45','ode15s','ode23','ode113')
% time integartion part
function xdot=eom(t,y0,M,C,K,F0,omega_f0)
Ft=F0*sin(omega_f0*t);
AA=[ 0 eye(size(M));-inv(M)*K -inv(M)*C];
xdot=AA*y0+[0 ;inv(M)*Ft];
end