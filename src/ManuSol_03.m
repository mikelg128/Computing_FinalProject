clear all

ax = 0; ay = 0;
bx = 2*pi; by = 2*pi;
Lambda = 0.5;
N = 100;
M = N;
cur=0;
e = 1;
w = 1; %Relaxation variable
iter = 1;
in = 1;
n = 4;
iterations = zeros(1,n);
err = zeros(1,n);
targeterror = 10^-10;

x = linspace(ax,bx,N+2);
y = linspace(ay,by,M+2);
% x = ax:delx:bx;
% y = ay:dely:by;

delx = x(2)-x(1); dely = y(2)-y(1);
delX = 1/(delx^2); delY = 1/(dely^2);
h = delx;
delta = 1/(Lambda*h^2 + 4);

[X,Y] = meshgrid(x,y);

v = zeros(N+2,N+2);
vprev = v;
vexact = v;
v(1,:) = 1; %<-- Neumann Condition
v(N+2,:) = cos(by*x); 
v(:,1) = 1;
v(:,N+2) = cos(bx*y);
G = cos(X.*Y).*((X.^2) + (Y.^2)) + Lambda*cos(X.*Y);
vexact = cos(X.*Y);
tic
while(e>=targeterror)
    vprev = v;
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
%             if i == 1
%                 v(i,j) = delta*(G(i,j)*h^2+(v(i,j-1)+2*v(i+1,j)+v(i,j+1)));
%             else
                v(i,j) = delta*(G(i,j)*h^2+(v(i,j-1)+v(i-1,j)+v(i+1,j)+v(i,j+1)));
%             end
            %v(i,j) = w*v(i,j) + (1-w)*vprev(i,j);
        end
%         v(1,j) = v(3,j);
    end
    iter = iter + 1;
    e = max(max(abs((vprev - v)./v)))*100;
end
[e1, I] = max(abs(vexact - v));
[e_abs_Linf, J] = max(e1); %L infinity abs error
e_abs_Linf
e_abs_L1 = (1/(N^3))*sum(sum(abs(vexact-v)))
e_abs_L2 = sqrt((1/(N^3))*sum(sum(abs(vexact-v))))
figure 
colormap('jet')
surf(X,Y,G)
title('Forcing Function G');
xlabel('X Axis')
ylabel('Y Axis')
figure
colormap(jet);
surf(X,Y,vexact);
title('Exact Solution for v(x,y)');
xlabel('X Axis')
ylabel('Y Axis')
figure 
colormap('jet')
surf(X,Y,v)
title(strcat(num2str(iter),' Iterations, Relative Error = ', num2str(e),'%'))
xlabel('X Axis')
ylabel('Y Axis')
zlabel('u(x,y)')

%%
%--------------------------Grid Convergence--------------------------------