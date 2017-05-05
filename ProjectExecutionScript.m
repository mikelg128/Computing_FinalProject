%clear all;
clc;
ax = 0; ay = 0;
bx = 2*pi; by = 2*pi;
Lambda = 0.5;
w = 1; %Relaxation variable
targeterror = 10^-4;
k=5; Nmax = 50;
in = 1;
f1=0;f2=0;f3=0;

for N = k:5:Nmax
tic

x = linspace(ax,bx,N+2);
y = linspace(ay,by,N+2);
delx = x(2)-x(1);  
h = delx;
[X,Y] = meshgrid(x,y);


%First Manufactured Solution-----------------------------------------------
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
v1time(in) = toc;


%Absolute Error Calculations
[ei, I] = max(abs(v1exact - v1));
[e_abs_Linf, J] = max(ei); %L infinity abs error
e_abs_L1 = (1/(N^2))*sum(sum(abs(v1exact-v1))); clear ei
e_abs_L2 = sqrt((1/(N^2))*sum(sum(abs(v1exact-v1))));
e_abs_array_Linf_1(in) = e_abs_Linf;
e_abs_array_L1_1(in) = e_abs_L1;
e_abs_array_L2_1(in) = e_abs_L2;
Narray1(in)=N;

if N >= 100 && f1 == 0
    f1=1;
    %Plot Manufactured sol (v1)
    figure 
    colormap('jet')
    surf(X,Y,v1)
    title(strcat(num2str(v1iter),' Iterations, Relative Error = ', num2str(e1),'%'))
    xlabel('X Axis')
    ylabel('Y Axis')
    zlabel(strcat('v1(x,y) at ',num2str(N),' Axis points'))
    %Plot Exact solution of manufactured problem
    figure
    colormap(jet);
    surf(X,Y,v1exact);
    title('Exact Solution for v1(x,y)');
    xlabel('X Axis')
    ylabel('Y Axis')
end



%Second Manufactured Solution----------------------------------------------
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
v2time(in) = toc;
%Absolute Error Calculations
[ei, I] = max(abs(v2exact - v2));
[e_abs_Linf, J] = max(ei); %L infinity abs error
e_abs_L1 = (1/(N^2))*sum(sum(abs(v2exact-v2))); clear ei
e_abs_L2 = sqrt((1/(N^2))*sum(sum(abs(v2exact-v2))));
e_abs_array_Linf_2(in) = e_abs_Linf;
e_abs_array_L1_2(in) = e_abs_L1;
e_abs_array_L2_2(in) = e_abs_L2;

if N >= 100 && f2 == 0
    f2=1;
    %Plot Manufactured sol (v2)
    figure 
    colormap('jet')
    surf(X,Y,v2)
    title(strcat(num2str(v2iter),' Iterations, Relative Error = ', num2str(e2),'%'))
    xlabel('X Axis')
    ylabel('Y Axis')
    zlabel(strcat('v2(x,y) at ',num2str(N),' Axis points'))
    %Plot Exact solution of manufactured problem
    figure
    colormap(jet);
    surf(X,Y,v2exact);
    title('Exact Solution for v2(x,y)');
    xlabel('X Axis')
    ylabel('Y Axis')
end


%%
%Main Solution-------------------------------------------------------------
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
utime(in) = toc

if N >= 100 && f3 == 0
    f3=1;
    %Plot
    figure 
    colormap('jet')
    surf(X,Y,u)
    title(strcat(num2str(iter),' Iterations, Relative Error = ', num2str(e),'%'))
    xlabel('X Axis')
    ylabel('Y Axis')
    zlabel('u(x,y)')
end
%N=N*2
in = in+1;
toc
end

%Convergence Plots
figure
semilogx(Narray1,e_abs_array_Linf_1,Narray1,e_abs_array_L1_1,Narray1,e_abs_array_L2_1);
title('Grid Convergence Analysis for v1(x,y)')
legend('L Infinity Absolute Error','L1 Absolute Error','L2 Absolute Error')
xlabel('Number of axis points');
ylabel('Absolute Error');
figure
loglog(Narray1, v1time)
title('Time to compute v1(x,y)')
xlabel('Time in Seconds')
ylabel('Number of axis points')

figure
semilogx(Narray1,e_abs_array_Linf_2,Narray1,e_abs_array_L1_2,Narray1,e_abs_array_L2_2);
title('Grid Convergence Analysis for v2(x,y)')
legend('L Infinity Absolute Error','L1 Absolute Error','L2 Absolute Error')
xlabel('Number of axis points');
ylabel('Absolute Error');
figure
loglog(Narray1, v2time)
title('Time to compute v2(x,y)')
xlabel('Time in Seconds')
ylabel('Number of axis points')

figure
loglog(Narray1, utime)
title('Time to compute u(x,y)')
xlabel('Time in Seconds')
ylabel('Number of axis points')
