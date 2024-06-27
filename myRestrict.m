function r2h = myRestrict(rh)
%Restrict function for cell centered mesh

    %get fine mesh dimensions
    M = size(rh,1)-2;
    N = size(rh,2)-2;

    %determine courser mesh size
    M_2h = M/2;
    N_2h = N/2;

    %delcare r2h with ghost cells
    r2h = zeros(M_2h+2, N_2h+2);

    for j = 2:N_2h+1
        for i = 2:M_2h+1
 
            % average fine mesh points to the coarse mesh
            r2h(i,j) = (1/4) * (rh(2*i-2,2*j-2) +rh(2*i-1,2*j-2) +rh(2*i-2,2*j-1)+ rh(2*i-1,2*j-1));
        end
    end
end