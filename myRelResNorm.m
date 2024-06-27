function rr = myRelResNorm(phi, f, h)
%calculates infinity norm of relative residual relative to f_inf. This is
%for a cell centered mesh

    %get 
    [M, N] = size(phi);
    N = N - 2;
    M = M - 2;

    %calculate the resisduals using fucntion coded in problem 2
    r = myResidual(phi, f, h);
    
    %calulate infintiy norm of r and f
    r_inf_norm = max(max(abs(r(2:M+1 , 2:N+1))));

    f_inf_norm= max(max(abs(f(2:M+1 , 2:N+1))));

    %calculate relative residual norm
    rr = r_inf_norm/f_inf_norm;

end