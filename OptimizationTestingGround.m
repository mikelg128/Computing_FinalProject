N = 100; Nmax = 2;


ax = 0; ay = 0;
bx = 2*pi; by = 2*pi;
Lambda = 0.5;

%in = 1;
x = linspace(ax,bx,N+2);
y = linspace(ay,by,N+2);
delx = x(2)-x(1);  
h = delx;
[X,Y] = meshgrid(x,y);
w = 1.9; %Relaxation variable
targeterror = 10^-5;


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
time = toc
e 
iter

