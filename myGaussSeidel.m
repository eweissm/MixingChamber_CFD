function phi = myGaussSeidel(phi, f, h, ninter)
    
    %get mesh size from input var (-2 to remove ghosts cells)
    [M,N] = size(phi);
    N = N-2;
    M=M-2;
    
    %step 1 --> define solution domain
    %not including ghost cells:
    % x = linspace(ys-h/2, (ys-h/2)+(N+2)*h, N+2)
    % x = linspace(xs-h/2, (xs-h/2)+(M+2)*h, M+2)

    %step 2--> Define mesh --> cell centered 
    %grid is M+2 x N+2 with corresponding i = 1...M+2, j = 1...N+2

    %step 3,4 --> skipping for now

    %step 5 --> apply BCs
    phi = bcGS(phi);

    %precalculate optimal omega
    rho_pj = 0.5*(cos(pi/M)+cos(pi/N));

    %Using same code as SOR algorithm, but with omega = 1
    omega =1;
    
    %Apply SOR algorithm to iteratively solve equation A*phi = b where A
    %and b are given from finite difference method
    
    for k = 1:ninter 
        for j = 2:N+1 %loop over interior j indices
            for i = 2:M+1 %loop over interior i indices
                %update Phi with SOR scheme
                phi(i,j)=phi(i,j)+omega*((phi(i,j-1)+phi(i-1,j)+phi(i+1,j)+phi(i,j+1))/4-h^2*f(i,j)/4-phi(i,j));
            end
        end
         %make sure that BC are still maintained
         phi = bcGS(phi);
    end
   
end