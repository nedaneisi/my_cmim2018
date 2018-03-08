%% Euler method for second order system 
function  [ T, yn] =Forward_euler(fun,xmin,xmax,delta,x0,v0)

T=(xmin:delta:xmax);
Xn=zeros(1,length(T));
Vn=zeros(1,length(T));


Xn(:,1) = x0; % initial conditions
Vn(:,1) = v0; % initial conditions

yn=[Xn;Vn];

for i=2:length(T)
    yn(:,i)=yn(:,i-1)+fun( T(i-1),yn(:,i-1))*delta;
    
end
end
