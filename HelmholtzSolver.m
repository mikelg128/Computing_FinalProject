function [ u, e, iter ] = HelmholtzSolver( Lambda, N, h, etarget, F, nbc, nflag, ebc, eflag, sbc, sflag, wbc, wflag, w )
%HelmholtzSolver takes in a forcing function, boundary conditions, and
%other parameters needed to solve the 2-Dimensional Helmholtz Equation:
%Laplacian(u) - Lambda*u = F
%   This solver uses the Gauss-Seidel Iteritive Method to find a solution
%   to the Helmholtz PDE
%Input Parameters:
%   Lambda:
%   N:
%   etarget:
%   F:
%   NBC-WBC: 
%   w:
%Output:
%   u:
%   e:
%   iter:
%   h:

    u = zeros(N+2,N+2);
    uprev = u;
    delta = 1/(Lambda*h^2 + 4);
    e = 1;
    iter = 1;
    if abs(1/delta) < 4
        error('Matrix is not diagonally dominant');
    end
    neue =0;
    neun =0;
    neuw =0;
    neus =0;

    if eflag == 'D'
        u(:,N+2) = ebc; 
    elseif eflag == 'N'
        neue = -ebc*2*h;
    end
    if nflag == 'D'
        u(N+2,:) = nbc; 
    elseif nflag == 'N'
        neun = -nbc*2*h;
    end
    if wflag == 'D'
        u(:,1) = wbc; 
    elseif wflag == 'N'
        neuw = -wbc*2*h;
    end
    if sflag == 'D'
        u(1,:) = sbc; 
    elseif sflag == 'N'
        neus = -sbc*2*h;
    end
    F = delta*F*h^2;
   
    while(e>=etarget)
        uprev = u;
        
        for j = 2:N+1
            if eflag == 'N' && j == N+1
                    u(:,j+1) = neue + u(:,j-1);
            end
            if wflag == 'N' && j == 2
                    u(:,j-1) = neuw + u(:,j+1);
            end
            for i = 2:2:N+1                
                if sflag == 'N' && i == 2
                    u(i-1,j) = neus + u(i+1,j);
                end
                if i == N+1
                    if nflag == 'N' && i == N+1
                        u(i+1,j) = neun + u(i-1,j);
                    end
                    u(i,j) = delta*(u(i,j-1)+u(i-1,j)+u(i+1,j)+u(i,j+1))-F(i,j);
                else
                    if nflag == 'N' && i == N
                        u(i+2,j) = neun + u(i,j);
                    end
                    u(i,j) = delta*(u(i,j-1)+u(i-1,j)+u(i+1,j)+u(i,j+1))-F(i,j);                    
                    u(i+1,j) = delta*(u(i+1,j-1)+u(i,j)+u(i+2,j)+u(i+1,j+1))-(F(i+1,j));
                end
                
                %u(i,j) = w*u(i,j) + (1-w)*uprev(i,j);
            end
        end
        
        iter = iter + 1;
        e = max(max(abs((uprev - u)./u)));
    end
    
end

