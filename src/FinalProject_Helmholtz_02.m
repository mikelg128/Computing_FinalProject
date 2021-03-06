clear all;
clc;

ax = 0; ay = 0;
bx = 2*pi; by = 2*pi;
Lambda = 0.5;
N = 100;
M = N;
cur=0;
ex=0;

x = linspace(ax,bx,N+2);
y = linspace(ay,by,M+2);
% x = ax:delx:bx;
% y = ay:dely:by;

delx = x(2)-x(1); dely = y(2)-y(1);
delX = 1/(delx^2); delY = 1/(dely^2);
h = delx;
delta = 1/(Lambda*h^2 + 4);

[X,Y] = meshgrid(x,y);

fb = y.*((by-y).^2);
gb = ((by-y).^2).*cos(pi*y/by);
F = sin(pi*((X-ax)/(bx-ax))).*cos((pi*0.5)*(2*(Y-ay)/(by-ay)+1));






% figure 
% colormap('jet')
% surf(X,Y,F)
%%
%-----------------------Calculations for Solution--------------------------
e = 1;
w = 1; %Relaxation variable
iter = 1;
in = 1;
n = 4;
iterations = zeros(1,n);
err = zeros(1,n);
targeterror = 10^-5;

%%
%-----------------------Boundary Conditions--------------------------------
uax = fb;
ubx = gb;
uay = fb(1) + ((x-ax)/bx-ax).*(gb(1)-fb(1));
dudy = 0;
nflag = 'N';
eflag = 'D';
sflag = 'D';
wflag = 'D';

[ u, e, iter, h ] = HelmholtzSolver( Lambda, N, h, targeterror, F, dudy, nflag, ubx, eflag, uay, sflag, uax, wflag, w );

% 
% % figure 
% % plot(x,uay,y,uax,y,ubx)
% 
% u = zeros(N+2,N+2);
% uprev = u;
% 
% u(1,:) = uay;
% %u(N+2,:) -> Neumann Boundary Condition
% u(N+2,1) = u(N,1);
% u(N+2,N+2) = u(N,N+2);
% u(:,1) = uax;
% u(:,N+2) = ubx;
% 
% ustart = u;
% du = u;
% if abs(1/delta) < 4
%     error('Matrix is not diagonally dominant');
% end


%%
%------------------------Grid Independance Study---------------------------


%%
%-----------------------Targeted Error Calculation-------------------------

% tic
% while(e>=targeterror)
%    
%     uprev = u;
%     %du = u;
%     if iter == 10^in
%         e
%         toc
%         iterations(in) = 10^in;
%         err(in) = e;
%         in = in + 1
%         
%         tic
%     end
%     for j = 2:N+1
%         for i = 2:1:N+1
%             if i == N+1
%                 u(i,j) = delta*(F(i,j)*h^2+(u(i,j-1)+2*u(i-1,j)+u(i,j+1)));
% %                 elseif i == N
% %                     u(i,j) = delta*(F(i,j)*h^2+(u(i,j-1)+u(i-1,j)+u(i+1,j)+u(i,j+1)));
% %                     k=i+1;
% %                     u(k,j) = delta*(F(k,j)*h^2+(u(k,j-1)+2*u(k-1,j)+u(k,j+1)));
%             else
%                 u(i,j) = delta*(F(i,j)*h^2+(u(i,j-1)+u(i-1,j)+u(i+1,j)+u(i,j+1)));
% 
% %                     k=i+1;
% %                     u(k,j) = delta*(F(k,j)*h^2+(u(k,j-1)+u(k-1,j)+u(k+1,j)+u(k,j+1)));
%             end
%             u(i,j) = w*u(i,j) + (1-w)*uprev(i,j);
%         end
%         u(N+2,j) = u(N,j);
%     end
%     iter = iter + 1;
%     e = max(max(abs((uprev - u)./u)))*100;
% end
% toc
% iter
%%
figure 
colormap('jet')
surf(X,Y,u)
title(strcat(num2str(iter),' Iterations, Relative Error = ', num2str(e),'%'))
xlabel('X Axis')
ylabel('Y Axis')
zlabel('u(x,y)')
%%
%---------------------------Error Analysis---------------------------------

% tic
% while(iter<=10^n)
%     u = ustart;
%     du = ustart;
%     in = in+1
%     iterations(in) = iter;
%     tic
%     while(cur<iter)
%         
%         uprev = u;
%         %du = u;
%         
%         for j = 2:N+1
%             for i = 2:1:N+1
%                 if i == N+1
%                     u(i,j) = delta*(F(i,j)*h^2+(u(i,j-1)+2*u(i-1,j)+u(i,j+1)));
% %                 elseif i == N
% %                     u(i,j) = delta*(F(i,j)*h^2+(u(i,j-1)+u(i-1,j)+u(i+1,j)+u(i,j+1)));
% %                     k=i+1;
% %                     u(k,j) = delta*(F(k,j)*h^2+(u(k,j-1)+2*u(k-1,j)+u(k,j+1)));
%                 else
%                     u(i,j) = delta*(F(i,j)*h^2+(u(i,j-1)+u(i-1,j)+u(i+1,j)+u(i,j+1)));
%                     
% %                     k=i+1;
% %                     u(k,j) = delta*(F(k,j)*h^2+(u(k,j-1)+u(k-1,j)+u(k+1,j)+u(k,j+1)));
%                 end
%                 u(i,j) = w*u(i,j) + (1-w)*uprev(i,j);
%             end
%             u(N+2,j) = u(N,j);
%         end
%         
%         cur = cur+1;
%     end
%     toc
%     cur = 0;
%     e = max(max(abs((uprev - u)./u)))*100;
%     err(in) = e;
%     figure 
%     colormap('jet')
%     surf(X,Y,u)
%     title(strcat(num2str(iter),' Iterations, Relative Error = ', num2str(e),'%'))
%     xlabel('X Axis')
%     ylabel('Y Axis')
%     zlabel('u(x,y)')
%     iter = iter*10;
%     
%     
% end
% toc



% figure
% semilogx(iterations, err)
% title('Relative Error vs. Number of Iterations')
% xlabel('Iterations')
% ylabel('% Relative Error')
% axis([10,max(iterations),-100,max(err)])




