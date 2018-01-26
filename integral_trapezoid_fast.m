function  result=integral_trapezoid_fast( fun, low_limit,no_splits,h )

result = 0;
i = 1:no_splits;
a=low_limit + (i-1)*h;
% b=low_limit + (i)*h;
 b=[a(2:end) low_limit + (no_splits)*h];
%   result = result+0.5*h*sum(feval(fun,a))+0.5*h*sum(feval(fun,b)); 
  result = result+0.5*h*sum(fun(a))+0.5*h*sum(fun(b));
end

