

function  [ T, Xn,Vn ] =Forward_euler_2order_SIE2(FUN,M,xmin,xmax,delta,x0,v0)

DOF=length(M);

T=(xmin:delta:xmax);
Xn=zeros(DOF,length(T));
Vn=zeros(DOF,length(T));


Xn(:,1) = x0; % initial conditions
Vn(:,1) = v0; % initial conditions

for i=2:length(T)
    Vn(:,i)=Vn(:,i-1)+FUN( T(i-1),Xn(:,i-1))*delta;
    Xn(:,i)=Xn(:,i-1)+Vn(:,i)*delta;
end
end
