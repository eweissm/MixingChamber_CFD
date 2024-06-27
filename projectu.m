function u = projectu(u, phi, dt)
    global t h;

     %get size of u
    M = size(u,1)-1; %node based in x dir
    N = size(u,2)-2; %cell centered in y dir
    
    for j = 2:N+1
        for i = 2:M
            %calculate u at n+1
            u(i,j) = u(i,j)-dt*(phi(i+1,j) - phi(i,j))/h; 

        end
    end 
    
    %apply BCs
    u = bcGhost_u(u, t);
end