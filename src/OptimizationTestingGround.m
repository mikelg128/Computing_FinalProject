N = 100; Nmax = 2;


ax = 0; ay = 0;
bx = 2*pi; by = 2*pi;
Lambda = 0.5;

%in = 1;
x = linspace(ax,bx,N);
y = linspace(ay,by,N);
delx = x(2)-x(1);  
h = delx;
[X,Y] = meshgrid(x,y);
w = 1.9; %Relaxation variable
targeterror = 10^-1;
in = 1
e1 = zeros(1,6);
e2=e1;e3=e1;iter1=e1;iter2=e1;iter3=e1;
time1=e1;time2=e1;time3=e1;
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

while targeterror>=10^-6
err(in)=targeterror;
%Call Solver
tic
[ u1, e1(in), iter1(in) ] = HelmholtzSolverJ( Lambda, N, h, targeterror, F, dudy, nflag, ubx, eflag, uay, sflag, uax, wflag );
time1(in) = toc
e1 
iter1
tic
[ u2, e2(in), iter2(in) ] = HelmholtzSolver( Lambda, N, h, targeterror, F, dudy, nflag, ubx, eflag, uay, sflag, uax, wflag );
time2(in) = toc
e2 
iter2
tic
[ u3, e3(in), iter3(in) ] = HelmholtzSolverSOR( Lambda, N, h, targeterror, F, dudy, nflag, ubx, eflag, uay, sflag, uax, wflag, w );
time3(in) = toc
e3 
iter3

targeterror=targeterror/10;
in = in +1;
end

