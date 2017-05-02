clear all;
clc;

ax = 0; ay = 0;
bx = 2*pi; by = 2*pi;
Lambda = 0.5;
N = 100;
M = N;

x = linspace(ax,bx,N);
y = linspace(ay,by,M);
% x = ax:delx:bx;
% y = ay:dely:by;

delx = x(2)-x(1); dely = y(2)-y(1);
delX = 1/(delx^2); delY = 1/(dely^2);

[X,Y] = meshgrid(x,y);

fb = y.*((by-y).^2);
gb = ((by-y).^2).*cos(pi*y/by);
F = sin(pi*((X-ax)/(bx-ax))).*cos((pi*0.5)*(2*(Y-ay)/(by-ay)+1));

% figure 
% colormap('jet')
% surf(x,y,F)

%Boundary Conditions
uax = fb;
ubx = gb;
uay = fb(1) + ((x-ax)/bx-ax).*(gb(1)-fb(1));
dudy = 0;

% figure 
% plot(x,uay,y,uax,y,ubx)