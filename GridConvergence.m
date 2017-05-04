clear all

ax = 0; ay = 0;
bx = 2*pi; by = 2*pi;
Lambda = 0.5;
N = 100;
M = N;
cur=0;
ex=0;
e = 1;
w = 1; %Relaxation variable
iter = 1;
in = 1;
n = 4;
iterations = zeros(1,n);
err = zeros(1,n);
targeterror = 10^-15;

tic
N = 2
in2=1;
while(N <= 10^n)
    x = linspace(ax,bx,N+2);
    y = linspace(ay,by,N+2);
    % x = ax:delx:bx;
    % y = ay:dely:by;

    delx = x(2)-x(1); dely = y(2)-y(1);
    delX = 1/(delx^2); delY = 1/(dely^2);
    h = delx;
    delta = 1/(Lambda*h^2 + 4);
    [X,Y] = meshgrid(x,y);
    u = zeros(N+2,N+2);
    uprev = u;
%     fb = y.*((by-y).^2);
%     gb = ((by-y).^2).*cos(pi*y/by);
%     F = sin(pi*((X-ax)/(bx-ax))).*cos((pi*0.5)*(2*(Y-ay)/(by-ay)+1));
%     uax = fb;
%     ubx = gb;
%     uay = fb(1) + ((x-ax)/bx-ax).*(gb(1)-fb(1));
%     dudy = 0;
    u(1,:) = 0;
    u(N+2,:) = (by^2)*x.^2; 
    u(:,1) = 0;
    u(:,N+2) = (bx^2)*y.^2;
    F = 2*(X.^2 + Y.^2) + Lambda*(X.^2).*(Y.^2);
    uexact = (X.^2).*(Y.^2);
    

%     u(1,:) = uay;
%     %u(N+2,:) -> Neumann Boundary Condition
%     u(N+2,1) = u(N,1);
%     u(N+2,N+2) = u(N,N+2);
%     u(:,1) = uax;
%     u(:,N+2) = ubx;

    ustart = u;
    du = u;
    if abs(1/delta) < 4
        error('Matrix is not diagonally dominant');
    end
    while(e>=targeterror)

        uprev = u;
        %du = u;
        if iter == 10^in
            e
            toc
            iterations(in) = 10^in;
            err(in) = e;
            in = in + 1

            tic
        end
        for j = 2:N+1
            for i = 2:1:N+1
%                 if i == N+1
%                     u(i,j) = delta*(F(i,j)*h^2+(u(i,j-1)+2*u(i-1,j)+u(i,j+1)));
%     %                 elseif i == N
%     %                     u(i,j) = delta*(F(i,j)*h^2+(u(i,j-1)+u(i-1,j)+u(i+1,j)+u(i,j+1)));
%     %                     k=i+1;
%     %                     u(k,j) = delta*(F(k,j)*h^2+(u(k,j-1)+2*u(k-1,j)+u(k,j+1)));
%                 else
                    u(i,j) = delta*(F(i,j)*h^2+(u(i,j-1)+u(i-1,j)+u(i+1,j)+u(i,j+1)));

    %                     k=i+1;
    %                     u(k,j) = delta*(F(k,j)*h^2+(u(k,j-1)+u(k-1,j)+u(k+1,j)+u(k,j+1)));
%                 end
                u(i,j) = w*u(i,j) + (1-w)*uprev(i,j);
            end
%             u(N+2,j) = u(N,j);
        end
        iter = iter + 1;
        e = max(max(abs((uprev - u)./u)))*100;
    end
    e_abs_Linf = max(max(abs(uexact - u))); %L infinity abs error
    e_abs_L1 = (1/(N^3))*sum(sum(abs(uexact-u)));
    e_abs_L2 = sqrt((1/(N^3))*sum(sum(abs(uexact-u))));
    e_abs_array_Linf(in2) = e_abs_Linf;
    e_abs_array_L1(in2) = e_abs_L1;
    e_abs_array_L2(in2) = e_abs_L2;
    Narray(in2)=N;
    harray(in2)=h;
    in2 = in2+1;
    
%     figure
%     mesh(X,Y,u)
%     title(strcat(num2str(N),' Axis Points, Absolute Error = ', num2str(e_abs)))
%     xlabel('X-Axis')
%     ylabel('Y-Axis')
    N=N*2
end
toc
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
