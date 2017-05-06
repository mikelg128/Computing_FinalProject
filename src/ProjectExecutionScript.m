clear all;
clc;
ax = 0; ay = 0;
bx = 2*pi; by = 2*pi;
Lambda = 0.5;
w = 1.5; %Relaxation variable
targeterror = 10^-4;
k=5; Nmax = 300;
in = 1;
f1=0;f2=0;f3=0;
man1=1;man2=1;main=1;



for N = k:20:Nmax
tic
N
Narray(in)=N;
x = linspace(ax,bx,N);
y = linspace(ay,by,N);
delx = x(2)-x(1);  
h = delx;
[X,Y] = meshgrid(x,y);

if man1 == 1
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
[ v1, e1, v1iter ] = HelmholtzSolverJ( Lambda, N, h, targeterror, F, nbc, nflag, ebc, eflag, sbc, sflag, wbc, wflag);
%[ v1, e1, v1iter ] = HelmholtzSolver( Lambda, N, h, targeterror, F, nbc, nflag, ebc, eflag, sbc, sflag, wbc, wflag);
%[ v1, e1, v1iter ] = HelmholtzSolverSOR( Lambda, N, h, targeterror, F, nbc, nflag, ebc, eflag, sbc, sflag, wbc, wflag, w);
v1time(in) = toc;


%Absolute Error Calculations
[ei, Iv1] = max(abs(v1 - v1exact));
[e_abs_Linf, J] = max(ei) %L infinity abs error
e_abs_L1 = (1/(N^2))*sum(sum(abs(v1-v1exact))) 
clear ei
e_abs_L2 = sqrt((1/(N^2))*sum(sum(abs(v1-v1exact).^2)))
e_abs_array_Linf_1(:,in) = [e_abs_Linf;N];
e_abs_array_L1_1(:,in) = [e_abs_L1;N];
e_abs_array_L2_1(:,in) = [e_abs_L2;N];


if N >= 100 && f1 == 0
    f1=1;
    %Plot Manufactured sol (v1)
    figure(1) 
    subplot(1,2,1)
    colormap('jet')
    contourf(X,Y,v1)
    title(strcat(num2str(v1iter),' Iterations, Relative Error = ', num2str(e1),'%'))
    xlabel('X Axis')
    ylabel('Y Axis')
    zlabel(strcat('v1(x,y) at ',num2str(N),' Axis points'))
    %Plot Exact solution of manufactured problem
    figure(1)
    subplot(1,2,2)
    colormap(jet);
    contourf(X,Y,v1exact);
    title('Exact Solution for v1(x,y)');
    xlabel('X Axis')
    ylabel('Y Axis')
end
end

if man2 == 1
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
[ v2, e2, v2iter ] = HelmholtzSolverJ( Lambda, N, h, targeterror, F, nbc, nflag, ebc, eflag, sbc, sflag, wbc, wflag);
%[ v2, e2, v2iter ] = HelmholtzSolver( Lambda, N, h, targeterror, F, nbc, nflag, ebc, eflag, sbc, sflag, wbc, wflag);
%[ v2, e2, v2iter ] = HelmholtzSolverSOR( Lambda, N, h, targeterror, F, nbc, nflag, ebc, eflag, sbc, sflag, wbc, wflag,w);

v2time(in) = toc;
%Absolute Error Calculations
[ei, Iv1] = max(abs(v2 - v2exact));
[e_abs_Linf, J] = max(ei); %L infinity abs error
e_abs_L1 = (1/(N^2))*sum(sum(abs(v2-v2exact))); clear ei
e_abs_L2 = sqrt((1/(N^2))*sum(sum(abs(v2-v2exact).^2)));
e_abs_array_Linf_2(:,in) = [e_abs_Linf;N];
e_abs_array_L1_2(:,in) = [e_abs_L1;N];
e_abs_array_L2_2(:,in) = [e_abs_L2;N];

if N >= 100 && f2 == 0
    f2=1;
    %Plot Manufactured sol (v2)
    figure(2)
    subplot(1,2,1)
    colormap('jet')
    contourf(X,Y,v2)
    title(strcat(num2str(v2iter),' Iterations, Relative Error = ', num2str(e2),'%'))
    xlabel('X Axis')
    ylabel('Y Axis')
    zlabel(strcat('v2(x,y) at ',num2str(N),' Axis points'))
    %Plot Exact solution of manufactured problem
    figure(2)
    subplot(1,2,2)
    colormap(jet);
    contourf(X,Y,v2exact);
    title('Exact Solution for v2(x,y)');
    xlabel('X Axis')
    ylabel('Y Axis')
end
end

%%
if main == 1
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
[ u, e, iter ] = HelmholtzSolverJ( Lambda, N, h, targeterror, F, dudy, nflag, ubx, eflag, uay, sflag, uax, wflag);
%[ u, e, iter ] = HelmholtzSolver( Lambda, N, h, targeterror, F, dudy, nflag, ubx, eflag, uay, sflag, uax, wflag);
%[ u, e, iter ] = HelmholtzSolverSOR( Lambda, N, h, targeterror, F, dudy, nflag, ubx, eflag, uay, sflag, uax, wflag,w);
utime(in) = toc

if N >= 100 && f3 == 0
    f3=1;
    %Plot
    figure(3) 
    colormap('jet')
    contourf(X,Y,u)
    title(strcat(num2str(iter),' Iterations, Relative Error = ', num2str(e),'%'))
    xlabel('X Axis')
    ylabel('Y Axis')
    zlabel('u(x,y)')
    figure(10)
    colormap('jet')
    surf(X,Y,u)
    title(strcat(num2str(iter),' Iterations, Relative Error = ', num2str(e),'%'))
    xlabel('X Axis')
    ylabel('Y Axis')
    zlabel('u(x,y)')
end
end

in = in+1;
toc
end

%Convergence Plots
if man1 == 1
figure(4)
subplot(2,1,1)
plot(e_abs_array_Linf_1(2,:),e_abs_array_Linf_1(1,:),e_abs_array_L1_1(2,:),e_abs_array_L1_1(1,:),e_abs_array_L2_1(2,:),e_abs_array_L2_1(1,:));
title(strcat('Grid Convergence Analysis for v1(x,y), Target Error: ', num2str(targeterror)))
legend('L Infinity Absolute Error','L1 Absolute Error','L2 Absolute Error')
xlabel('Number of axis points');
ylabel('Absolute Error');
[minv1, Iv1] = min(e_abs_array_Linf_1(1,:));
text(e_abs_array_Linf_1(2,Iv1),e_abs_array_Linf_1(1,Iv1),strcat('\leftarrow Minimum Error at N = ',num2str(e_abs_array_Linf_1(2,Iv1))))
figure(4)
subplot(2,1,2)
loglog(Narray, v1time)
title('Time to compute v1(x,y)')
xlabel('Time in Seconds')
ylabel('Number of axis points')
end
if man2 == 1
figure(5)
subplot(2,1,1)
semilogx(e_abs_array_Linf_2(2,:),e_abs_array_Linf_2(1,:),e_abs_array_L1_2(2,:),e_abs_array_L1_2(1,:),e_abs_array_L2_2(2,:),e_abs_array_L2_2(1,:));
title('Grid Convergence Analysis for v2(x,y)')
title(strcat('Grid Convergence Analysis for v2(x,y), Target Error: ', num2str(targeterror)))
legend('L Infinity Absolute Error','L1 Absolute Error','L2 Absolute Error')
xlabel('Number of axis points');
ylabel('Absolute Error');
[minv2, Iv2] = min(e_abs_array_Linf_2(1,:));
text(e_abs_array_Linf_2(2,Iv2),e_abs_array_Linf_2(1,Iv2),strcat('\leftarrow Minimum Error at N = ',num2str(e_abs_array_Linf_2(2,Iv2))))

figure(5)
subplot(2,1,2)
loglog(Narray, v2time)
title('Time to compute v2(x,y)')
xlabel('Time in Seconds')
ylabel('Number of axis points')
end
if main == 1
figure(6)
loglog(Narray, utime)
title('Time to compute u(x,y)')
xlabel('Time in Seconds')
ylabel('Number of axis points')
end





