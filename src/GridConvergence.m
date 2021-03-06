clear all

ax = 0; ay = 0;
bx = 2*pi; by = 2*pi;
Lambda = 0.5;
w = 1; %Relaxation variable
in = 1;

iterations = zeros(1,n);
err = zeros(1,n);
targeterror = 10^-4;

elapsed=0;
N = 2
while(N <= 2^8)
    e=1;
    x = linspace(ax,bx,N+2);
    y = linspace(ay,by,N+2);

    delx = x(2)-x(1); dely = y(2)-y(1);
    delX = 1/(delx^2); delY = 1/(dely^2);
    h = delx;
    [X,Y] = meshgrid(x,y);
    uexact = u;


    e
    [e1, I] = max(abs(uexact - u));
    [e_abs_Linf, J] = max(e1); %L infinity abs error
    e_abs_L1 = (1/(N^3))*sum(sum(abs(uexact-u)));
    e_abs_L2 = sqrt((1/(N^3))*sum(sum(abs(uexact-u))));
    e_abs_array_Linf(in2) = e_abs_Linf;
    e_abs_array_L1(in2) = e_abs_L1;
    e_abs_array_L2(in2) = e_abs_L2;
    Narray(in2)=N;
    harray(in2)=h;
    N=N*2
    elapsed = toc
    time(in2) = elapsed;
    in2 = in2+1;
end

figure
semilogx(Narray,e_abs_array_Linf,Narray,e_abs_array_L1,Narray,e_abs_array_L2);
title('Grid Convergence Analysis')
legend('L Infinity Absolute Error','L1 Absolute Error','L2 Absolute Error')
xlabel('Number of axis points');
ylabel('Absolute Error');
figure
semilogx(harray,e_abs_array_Linf,harray,e_abs_array_L1,harray,e_abs_array_L2);
title('Grid Convergence Analysis')
legend('L Infinity Absolute Error','L1 Absolute Error','L2 Absolute Error')
xlabel('Step Size');
ylabel('Absolute Error');
figure
loglog(Narray, time)
title('Time to compute')
xlabel('Time in Seconds')
ylabel('Number of axis points')

