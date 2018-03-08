%% Euler method for second order system 
function  [ T, Xn,Vn ] =Forward_euler_2order_SIE(fun,xmin,xmax,delta,x0,v0)

T=(xmin:delta:xmax);
Xn(1) = x0; % initial conditions
Vn(1) = v0; % initial conditions



for i=2:length(T)
    Vn(i)=Vn(i-1)+fun( T(i-1),Xn(i-1))*delta;
    Xn(i)=Xn(i-1)+Vn(i-1)*delta;
    
end
end
