%clear all;
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
w = 1; %Relaxation variable
targeterror = 10^-10;

%First Manufactured Solution
%Boundary conditions + exact sol for manufactured problem
sbc = 0;
nbc = (by^2)*x.^2; 
wbc = 0;
ebc = (bx^2)*y.^2;
sflag = 'D';
nflag = 'D';
eflag = 'D';
wflag = 'D';
F = 2*(X.^2 + Y.^2) - Lambda*(X.^2).*(Y.^2);
v1exact = (X.^2).*(Y.^2);
tic
[ v1, e1, v1iter ] = HelmholtzSolver( Lambda, N, h, targeterror, F, nbc, nflag, ebc, eflag, sbc, sflag, wbc, wflag, w );
v1time = toc
%Plot Manufactured sol
figure 
colormap('jet')
surf(X,Y,v1)
title(strcat(num2str(iter),' Iterations, Relative Error = ', num2str(e),'%'))
xlabel('X Axis')
ylabel('Y Axis')
zlabel('v1(x,y)')
%Plot Exact solution of manufactured problem
figure
colormap(jet);
surf(X,Y,v1exact);
title('Exact Solution for v1(x,y)');
xlabel('X Axis')
ylabel('Y Axis')

%Second Manufactured Solution
%BCs + exact 
sbc = 1; 
nbc = cos(by*x); 
wbc = 1;
ebc = cos(bx*y);
sflag = 'D';
nflag = 'D';
eflag = 'D';
wflag = 'D';
F = -cos(X.*Y).*((X.^2) + (Y.^2)) - Lambda*cos(X.*Y);
v2exact = cos(X.*Y);
tic
[ v2, e2, v2iter ] = HelmholtzSolver( Lambda, N, h, targeterror, F, nbc, nflag, ebc, eflag, sbc, sflag, wbc, wflag, w );
v2time = toc
%Plot Manufactured sol
figure 
colormap('jet')
surf(X,Y,v2)
title(strcat(num2str(iter),' Iterations, Relative Error = ', num2str(e),'%'))
xlabel('X Axis')
ylabel('Y Axis')
zlabel('v2(x,y)')
%Plot Exact solution of manufactured problem
figure
colormap(jet);
surf(X,Y,v2exact);
title('Exact Solution for v2(x,y)');
xlabel('X Axis')
ylabel('Y Axis')



%%

%Boundary conditions for main problem
fb = y.*((by-y).^2);
gb = ((by-y).^2).*cos(pi*y/by);
F = sin(pi*((X-ax)/(bx-ax))).*cos((pi*0.5)*(2*(Y-ay)/(by-ay)+1));

uax = fb;
ubx = gb;
uay = fb(1) + ((x-ax)/bx-ax).*(gb(1)-fb(1));
dudy = 0;
nflag = 'N';
eflag = 'D';
sflag = 'D';
wflag = 'D';

%Call Solver
tic
[ u, e, iter ] = HelmholtzSolver( Lambda, N, h, targeterror, F, dudy, nflag, ubx, eflag, uay, sflag, uax, wflag, w );
utime = toc

%Plot
figure 
colormap('jet')
surf(X,Y,u)
title(strcat(num2str(iter),' Iterations, Relative Error = ', num2str(e),'%'))
xlabel('X Axis')
ylabel('Y Axis')
zlabel('u(x,y)')
