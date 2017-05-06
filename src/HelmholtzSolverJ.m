function [ u, e, iter ] = HelmholtzSolverJ( Lambda, N, h, etarget, F, nbc, nflag, ebc, eflag, sbc, sflag, wbc, wflag )
%HelmholtzSolver takes in a forcing function, boundary conditions, and
%other parameters needed to solve the 2-Dimensional Helmholtz Equation:
%Laplacian(u) - Lambda*u = F
%   This solver uses the Jacobi Iteritive Method to find a solution
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

    u = zeros(N,N);
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
        u(:,N) = ebc; 
    elseif eflag == 'N'
        neue = -ebc*2*h;
    end
    if nflag == 'D'
        u(N,:) = nbc; 
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
        
        for j = 2:N-1
            if eflag == 'N' && j == N-1
                    u(:,j+1) = neue + uprev(:,j-1);
            end
            if wflag == 'N' && j == 2
                    u(:,j-1) = neuw + uprev(:,j+1);
            end
            for i = 2:2:N-1                
                if sflag == 'N' && i == 2
                    u(i-1,j) = neus + uprev(i+1,j);
                end
                if i == N-1
                    if nflag == 'N' && i == N-1
                        u(i+1,j) = neun + uprev(i-1,j);
                    end
                    u(i,j) = delta*(uprev(i,j-1)+uprev(i-1,j)+uprev(i+1,j)+uprev(i,j+1))-F(i,j);
                else
                    if nflag == 'N' && i == N-2
                        u(i+2,j) = neun + uprev(i,j);
                    end
                    u(i,j) = delta*(uprev(i,j-1)+uprev(i-1,j)+uprev(i+1,j)+uprev(i,j+1))-F(i,j);                    
                    u(i+1,j) = delta*(uprev(i+1,j-1)+uprev(i,j)+uprev(i+2,j)+uprev(i+1,j+1))-(F(i+1,j));
                end
                
                
            end
        end
        
        iter = iter + 1;
        e = max(max(abs((uprev - u)./u)));
    end
    
end
