

function  [ T, yn ] =Forward_euler_2order_SIE1(K,M,xmin,xmax,delta,x0,v0)

DOF=length(M);
C=zeros(size(M));
A=[zeros(size(M)) eye(size(M)); -M\K  -M\C]; % general form of state space matix
FUN1=@(t, x) A*x;


T=(xmin:delta:xmax);
Xn=zeros(DOF,length(T));
Vn=zeros(DOF,length(T));


Xn(:,1) = x0; % initial conditions
Vn(:,1) = v0; % initial conditions
yn=[Xn;Vn];


for i=2:length(T)
    yn(:,i)=yn(:,i-1)+FUN1( T(i-1),yn(:,i-1))*delta;
end


end
