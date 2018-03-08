
function xdot=myode(t,x,M,K,C)

% if t < 1 % this was an example of how to apply force in spceific interval
% machnie dynamic exercise 6
%     f=1;
% else
%     f=0;
% end
f=0;
F=[f;f;f];

A=[zeros(size(M)) eye(size(M)); -M\K  -M\C];
xdot=A*x+[zeros(length(M(:,1)),1);M\F];
    
end
