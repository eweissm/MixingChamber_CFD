function [Hu,Hv] = hyperbolic_uv_2D(u,v)
    global h
    %Solves for convective terms in the 2d- viscous burgers' eq on a
    %staggered mesh (see mod 7 for details) when moved to the RHS of the equation

    %get M and N from v
    M = size(v,1)-2; % in x dir
    N = size(v,2)-1;  %y dir

    %lets declare Hu and Hv arrays for staggered mesh
    Hu = zeros(M+1, N+2);
    Hv = zeros(M+2, N+1);

    %calculate Hu on interior points only
    for j = 2:N+1
        for i = 2:M
            
            duu_dx =-( (0.5*(u(i+1,j)+u(i,j)))^2 - (0.5*(u(i,j)+u(i-1,j)))^2 ) / h;
            duv_dy =-( (0.5*(u(i,j)+u(i,j+1)))*(0.5*(v(i,j)+v(i+1,j))) - (0.5*(u(i,j-1)+u(i,j)))*(0.5*(v(i,j-1)+v(i+1,j-1))) ) / h;

            Hu(i,j) = duu_dx + duv_dy;
        end
    end

    %calculate Hv on interior points only
    for j = 2:N
        for i = 2:M+1
            
            dvv_dy =-((0.5*(v(i,j) + v(i,j+1)))^2 - (0.5*(v(i,j)+v(i,j-1)))^2 ) / h;
            duv_dx =-( (0.5*(v(i,j)+v(i+1,j)))*(0.5*(u(i,j)+u(i,j+1))) - (0.5*(v(i,j)+v(i-1,j)))*(0.5*(u(i-1,j)+u(i-1,j+1))) ) / h;

            Hv(i,j) = dvv_dy + duv_dx;
        end
    end

end