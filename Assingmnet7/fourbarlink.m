


function [Phi2,Phi3,Phidot2,Phidot3, Phidotdot2,Phidotdot3]=fourbarlink(Inp,x0,xdot0,xdotdot0);



X_result= fsolve(@(x) constrains(x,Inp),x0) ;
Phi2=X_result(1);
Phi3=X_result(2);

V_result= fsolve(@(Xdot) constrains_V(Xdot,Inp,Phi2,Phi3),xdot0) ;
Phidot2=V_result(1);
Phidot3=V_result(2);

a_result= fsolve(@(Xdotdot) constrains_a(Xdotdot,Inp,Phi2,Phi3,Phidot2,Phidot3),xdotdot0) ;
Phidotdot2=a_result(1);
Phidotdot3=a_result(2);

function F=constrains(X,Inp)
F=[Inp.l1*cos(Inp.phi1)+Inp.l2*cos(X(1))-Inp.l3*cos(X(2))-Inp.d1;
Inp.l1*sin(Inp.phi1)+Inp.l2*sin(X(1))-Inp.l3*sin(X(2))-Inp.d2];
end

function F=constrains_V(Xdot,Inp,Phi2,Phi3)
F =[Inp.l3*sin(Phi3)*Xdot(2) - Inp.l2*sin(Phi2)*Xdot(1) - Inp.l1*sin(Inp.phi1)*Inp.phidot1;
 Inp.l1*cos(Inp.phi1)*Inp.phidot1 + Inp.l2*cos(Phi2)*Xdot(1) - Inp.l3*cos(Phi3)*Xdot(2)];
end

function F=constrains_a(Xdotdot,Inp,Phi2,Phi3,Phidot2,Phidot3)
F =[
 
Inp.l3*sin(Phi3)*Xdotdot(2) - Inp.l2*sin(Phi2)*Xdotdot(1) - Inp.l1*sin(Inp.phi1)*Inp.Phidotdot1 - Inp.l1*cos(Inp.phi1)*Inp.phidot1^2 - Inp.l2*cos(Phi2)*Phidot2^2 + Inp.l3*cos(Phi3)*Phidot3^2;
 
Inp.l1*cos(Inp.phi1)*Inp.Phidotdot1 + Inp.l2*cos(Phi2)*Xdotdot(1) - Inp.l3*cos(Phi3)*Xdotdot(2) - Inp.l1*sin(Inp.phi1)*Inp.phidot1^2 - Inp.l2*sin(Phi2)*Phidot2^2 + Inp.l3*sin(Phi3)*Phidot3^2];
 
end
end