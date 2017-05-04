clear all;
clc;

ax = 0; ay = 0;
bx = 2*pi; by = 2*pi;
Lambda = 0.5;
N = 100;

x = linspace(ax,bx,N+2);
y = linspace(ay,by,N+2);
delx = x(2)-x(1);  
h = delx;
[X,Y] = meshgrid(x,y);
fb = y.*((by-y).^2);
gb = ((by-y).^2).*cos(pi*y/by);
F = sin(pi*((X-ax)/(bx-ax))).*cos((pi*0.5)*(2*(Y-ay)/(by-ay)+1));

w = 1; %Relaxation variable
targeterror = 10^-5;

uax = fb;
ubx = gb;
uay = fb(1) + ((x-ax)/bx-ax).*(gb(1)-fb(1));
dudy = 0;
nflag = 'N';
eflag = 'D';
sflag = 'D';
wflag = 'D';

[ u, e, iter, h ] = HelmholtzSolver( Lambda, N, h, targeterror, F, dudy, nflag, ubx, eflag, uay, sflag, uax, wflag, w );

figure 
colormap('jet')
surf(X,Y,u)
title(strcat(num2str(iter),' Iterations, Relative Error = ', num2str(e),'%'))
xlabel('X Axis')
ylabel('Y Axis')
zlabel('u(x,y)')
