

clear all 
clc

%   Function that calculate the integral of 1D continous function using
%   well-known trapezoidal rule. Naive approach.
%   fun - handle of a function to integrate,
%   low_limit - lower limit of an integral,
%   up_limit - upper limit of an integral,
%   no_splits - number of trapezoids.

%% integral_trapezoid given in the assignment
% The aim of the exercise is to modifiy the code given the the computional method in a way that speed up the caclulation 10 times. The matlab cod was added to the repositiry in the git hub
% %The function sin(x) over the intervall of 0-pi is tested'¤  
tic

fun = @(x) sin(x);
low_limit=0;
up_limit=pi;
no_splits=100;
% result_main=integral_trapezoid( fun, low_limit, up_limit, no_splits )

%% modifed function to reduce the time
% % result_fast=integral_trapezoid_fast( fun, low_limit,no_splits,h ) % when the commnads directly inserted into main file the calcluation time is reduced

h = (up_limit - low_limit) / no_splits;
result = 0;
i = 1:no_splits;
a=low_limit + (i-1)*h;
% b=low_limit + (i)*h;
b=[a(2:end) low_limit + (no_splits)*h];
result_fast = result+0.5*h*sum(fun(a))+0.5*h*sum(fun(b));
toc

return
%% 2D integral
N = 20;
xmin=0;
xmax=2;
ymin=0;
ymax=1;

x = linspace(xmin,xmax,N);
y = linspace(ymin,ymax,N);
dx=(xmax-xmin)/(N-1);
dy=(ymax-ymin)/(N-1);

[x,y] = meshgrid(x,y);
fun = x.^2.*y;
fun(2:end-1,:) = fun(2:end-1,:)*2;
fun(:,2:end-1) = fun(:,2:end-1)*2;
result_2d= sum(fun(:))*dx*dy/4

trapz


fun = @(x,y) x.^2.*y;
result_2d = integral2(fun,xmin,xmax,ymin,ymax)