function phi = myMultigrid(phi, f, h)
%multigrid recursive function for cell centered mesh
%essentially the same code from lecture notes

    %set some constants:
    n = 1;
    nc = 3;
    
    %find the size of M and N
    M = size(phi,1)-2;
    N = size(phi,2)-2;

    %perform n iters of GS algorithm
    phi = myGaussSeidel(phi,f,h,n);

    %if we are not on the coarsest possible mesh:
    if mod(M,2)==0 && mod(N,2)==0

        rh = myResidual(phi,f,h); %get phi residuals
        r2h = myRestrict(rh);   %restrict fine mesh
        e2h = zeros(M/2+2,N/2+2); %declare e2h cell centered mesh
        e2h = myMultigrid(e2h,r2h,2*h); %recursively run myMultigrid
        eh = myProlong(e2h); %prolong coarse mesh
        phi = phi + eh; % add error to phi to improve phi 
        phi = bcGS(phi); %reapply BCs
        phi = myGaussSeidel(phi,f,h,n); %run GS again for n iters

    else
        %perform nc iters of GS algorithm at coarsest mesh
        phi = myGaussSeidel(phi,f,h,nc);

    end

end