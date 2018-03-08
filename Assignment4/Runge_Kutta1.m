function  [ T, yn ] =Runge_Kutta1(K,M,xmin,xmax,delta,x0,v0)
 
DOF=length(M);

C=zeros(size(M));

% if there is ax=b x=inv(a)*b but for  both execution time and numerical accuracy use x=a\b    mathwork, Grzeg comment! 
% be careful: for Rungkutta use the method given in inman book, page 89. numerical simulation and time responce. You can not use the approach  given  for semi explicit method,  When I tried it I get big error. 

A=[zeros(size(M)) eye(size(M)); -M\K  -M\C]; % general form of state space matix
FUN1=@(t, x) A*x;

T=(xmin:delta:xmax);
Xn=zeros(DOF,length(T));
Vn=zeros(DOF,length(T));
 
Xn(:,1) = x0; % initial conditions
Vn(:,1) = v0; % initial conditions
yn=[Xn;Vn];


for i=2:length(T)
    k1=FUN1( T(i-1),yn(:,i-1))*delta;
    k2=FUN1( T(i-1)+delta/2 ,yn(:,i-1) + k1/2 )*delta;
    k3=FUN1( T(i-1)+delta/2 ,yn(:,i-1) + k2/2 )*delta;
    k4=FUN1( T(i-1)+delta ,yn(:,i-1) + k3 )*delta;
  
    yn(:,i)=yn(:,i-1)+k1/6+k2/3+k3/3+k4/6;
   
end

end