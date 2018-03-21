
close all
clear all
clc

Inp.l1=0.2,
Inp.l2=0.4;
Inp.l3=0.3;
Inp.d1=0.35;
Inp.d2=0.1;

Inp.phi1=2.36;
Inp.phidot1=6.28;
Inp.Phidotdot1 =0;
x0 = [0 ; 2];
xdot0 = [0 ; 2];
xdotdot0 = [0 ; 2];

% PHI1=[2.36,2.52,2.67,2.83,8.49,8.64]; % different range for the angle phi1

 PHI1=[2.36,2.52,2.67,2.83,3,3.5,4,5,6,7,8,8.49,8.64]; % different range for the angle phi1

for i=1:length(PHI1)
    Inp.phi1=PHI1(i);
    
    [ Phi2(i),Phi3(i),Phidot2(i),Phidot3(i), Phidotdot2(i),Phidotdot3(i)]=fourbarlink(Inp,x0,xdot0,xdotdot0);
end

figure; hold on
x1=-Inp.l1*cos(PHI1);
y1=Inp.l1*sin(PHI1);
x2=Inp.d2-Inp.l3*cos(Phi3);
y2=Inp.d1+Inp.l3*sin(Phi3);
X0=zeros(size(x1));
Y0=zeros(size(x1));
xend=zeros(size(x1));
yend=zeros(size(x1));

X0(:)=0; Y0(:)=0;
xend(:)=Inp.d1;     yend(:)=Inp.d2;
XX=[X0;x1; x2;xend];
YY=[Y0;y1; y2;yend];

for i=1:length(PHI1)
    plot(XX(:,i),YY(:,i))
    
    legendInfo{i} = ['{\phi}1 = ' num2str(PHI1(i))]; % or whatever is appropriate
    
    
end
legend(legendInfo)

%
%
Result=[PHI1' Phi2' Phi3'  Phidot2' Phidot3' Phidotdot2' Phidotdot3']

figure; hold on
plot(Phidot2)
plot(Phidot3)
title('velcotiy of rotation angle ')
legend('phi2','phi3')

figure; hold on
plot(Phidotdot2)
plot(Phidotdot3)
title('acceleration of rotation angle')
legend('phi2','phi3')
