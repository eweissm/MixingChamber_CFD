function [u, v]= projectV(u, v, phi, dt)
    global t h;

     %get size of M and N
    M = size(u,1)-1; 
    N = size(u,2)-2; 
    
    for j = 2:N+1
        for i = 2:M
            %calculate u at n+1
            u(i,j) = u(i,j)-dt*(phi(i+1,j) - phi(i,j))/h; 
        end
    end 

    for j = 2:N
        for i = 2:M+1
            %calculate v at n+1
            v(i,j) = v(i,j)-dt*(phi(i,j+1) - phi(i,j))/h; 
        end
    end 
    
    %apply BCs
    u = bcGhost_u(u, t);
    v = bcGhost_v(v, t);
end