close all
clear all
clc


%% indata part 3,4,5
xmin= 0; xmax=20; delta = 0.05; % stimaluation time and time step
x0=1; v0=0; % initial conidtion
m=1;k=1; % mass and spring 
c=0; % damping
%% analytical solution
syms x(t) 
Dx=diff(x,t);
DDx=diff(Dx,t);

eq=m*DDx+k*x;
Analytical_solution= dsolve(eq==0, [Dx(0) == v0, x(0) == x0] )
% ezplot(Analytical_solution, [0, xmax]) ; % plot of symbolic function for specific range

%% using forward euler for second order numerical integration
fun = @(t, x) (-k/m)*x;
[T,Xn,Vn]=Forward_euler_2order_SIE(fun,xmin,xmax,delta,x0,v0);

potential_energy=0.5*k.*Xn.^2;
kinetic_energy=0.5*m.*Vn.^2;

Analytical_Solution=subs(Analytical_solution,t,T);
Analytical_Solution_V=subs(diff(Analytical_solution),t,T);


A0=[0 1; -m\k  -m\c]; % general form of state space matix t
fun1=@(t, x) A0*x; % input function for forward euler method
[ T, yn] =Forward_euler(fun1,xmin,xmax,delta,x0,v0);

%% using Runge_Kutta for second order numerical integration
fun1 = @(t, x) 0.5*(-k/m)*x; % dv/dt
fun2 = @(t, v) v; % dx/dt


[ T, Xn_RK4] =Runge_Kutta1(k,m,xmin,xmax,delta,x0,v0);

%% plot result first part
Figure_setup

Tit0='Forward Euler-position';
Tit00='Forward Euler-velocity';
Tit='SIE-position';
Tit1='SIE-Velocity';
Tit2='RK4-position';
plotfigure(T,Analytical_Solution,yn(1,:),Tit0) % position- forward euler and analyitical
plotfigure(T,Analytical_Solution_V,yn(2,:),Tit00) % velocity- forward euler and analyitical
plotfigure(T,Analytical_Solution,Xn,Tit) % position- SIE and analyitical
plotfigure(T,Analytical_Solution_V,Vn,Tit1)% velocity- SIE and analyitical
plotfigure(T,Analytical_Solution,Xn_RK4(1,:),Tit2)% position- RK4 and analyitical

plotfigure_firstsection_part2

%% Assingmet 4 3 mass
m1=2; m2=1.5; m3=1; 
k1=20000; k2=15000; k3=10000; 

M=[m1 0 0 ; 0 m2 0 ; 0 0 m3];
K=[k1+k2 -k2 0;  -k2 k2+k3 -k3; 0 -k3 k3];

%% using forward euler and runge gutta for MULTIDOFS
FUN = @(t, x) (-K/M)*x; % dv/dt
xmin=0; xmax=1;delta=0.0001;

[ T, Xn_3MASS] =Forward_euler_2order_SIE1(K,M,xmin,xmax,delta,x0,v0);
% [ T,  Xn_3MASS,vn_3MASS] =Forward_euler_2order_SIE2(FUN,M,xmin,xmax,delta,x0,v0);
[ T, Xn_RK43MASS] =Runge_Kutta1(K,M,xmin,xmax,delta,x0,v0);
 %% modal solution
x0=[x0;x0;x0];
xdot0=[v0;v0;v0];
t=xmin:delta:xmax;


[u d ]=eig(K,M);
omega=sqrt(d)


L=chol(M,'lower'); %(M^(1/2))
k_tilde=inv(L)*K*inv(L)';% mass normalized stiffness matrix

% solve symmetric eigenvalue problem
[P landa ]=eig(k_tilde);
XX=P'*P; % check that mode shape matrix are orthogonal
vv=P(:,1).*P(:,2);
Km=P'*k_tilde*P; % Modal stiffness matrix
S=inv(L)*P ; % S= mode shape matrix , P = matrix of eigenvectors

% convert initial conidtion to modal coordinate
r0=inv(S)*x0  ;
rdot0=inv(S)*xdot0;

Omega=sqrt(landa);
r=zeros(length(M),length(t));

for i=1:length(M)
    r(i,:)=r0(i)*cos(Omega(i,i).*t) + (rdot0(i)/Omega(i,i))*sin( Omega(i,i).*t);
end
% convert modal corrdinate back to physical coordinate
X=S*r;

% analytical solution
for i=1:length(M)
    x_analytical(i,:)=x0(i)*cos(omega(i,i).*t) + (xdot0(i)/omega(i,i))*sin( omega(i,i).*t);
end


%% use ode 45 for multidegree 
% in this exercise we dont have damping in the system, therefore, damping matrix=zero
C=zeros(size(M));

% vector of initial conditions
y0=[x0 ; xdot0]';

% integration paramaters
err=1e-5;
options = odeset('RelTol',err,'MaxStep',1e-3,'Stats','on');

% numerical integration using ode45
[t,y_45]=ode45(@myode,t,y0,options,M,K,C);
[t,y_15s]=ode15s(@myode,t,y0,options,M,K,C);


% plot result 3mass
Tit3='SIE';
Tit4='RK4';
Tit5='modal responce';
Tit6='physical responce';
Tit7='ode45';
Tit8='ode15s';

plotfigure3mass(T,Xn_3MASS,Tit3) % euler method 
plotfigure3mass(T,Xn_RK43MASS,Tit4)% RK4
plotfigure3mass(T,r,Tit5)% modal repsonce-modal solution
plotfigure3mass(T,X,Tit6)% physical responce-modal solution
plotfigure3mass(T,y_45',Tit7)% ode45
plotfigure3mass(T,y_15s',Tit8)% ode15s


Tit9='position mass 1 -3dof';
Tit10='position mass 2 -3dof';
Tit11='position mass 3 -3dof';

plotfigure3mass_differentmethod(T,Xn_3MASS,Xn_RK43MASS,X,y_45,y_15s,1,Tit9) % compare different method for m1 
plotfigure3mass_differentmethod(T,Xn_3MASS,Xn_RK43MASS,X,y_45,y_15s,2,Tit10)% compare different method for m2 
plotfigure3mass_differentmethod(T,Xn_3MASS,Xn_RK43MASS,X,y_45,y_15s,3,Tit11)% compare different method for m3 

