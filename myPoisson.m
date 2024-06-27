function [phi, Linf, iter] = myPoisson(phi, f, h, nIterMax, epsilon)
%function to solve elliptical Poisson Equation using GaussSiedel Multigrid
%method on cell centered mesh


    phi = myMultigrid(phi, f, h); %perform at least 1 multigrid iteration
    Linf = myRelResNorm(phi, f, h); %calculate relative residual norm

    iter = 1; %we've performed 1 iteration

    while(iter < nIterMax &&  Linf>epsilon)
        phi = myMultigrid(phi, f, h); %perform multigrid iteration
        Linf = myRelResNorm(phi, f, h); %calculate relative residual norm
        iter = iter+1; %increment iteration counter
    end
end