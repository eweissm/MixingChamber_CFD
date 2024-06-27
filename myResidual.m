function r = myResidual(phi, f, h)

    %get mesh size from input var (-2 to remove ghosts cells)
    [M,N] = size(phi);
    N = N - 2;
    M = M - 2;
    
    % declare / initialized r array
    r = zeros(size(phi));

    for j = 2:N+1 %loop over interior j cells
        for i = 2:M+1 %loop over interior i cells
            %use a second order central difference method to calculate residuals
            %r_ij = f_ij - (phi''_x_ij + phi''_y_ij)

            d2phidx2_ij= (phi(i+1,j) - 2* phi(i,j) +phi(i-1,j)) / (h^2) ; 

            d2phidy2_ij= (phi(i,j+1) - 2* phi(i,j) +phi(i,j-1)) / (h^2) ;

            r(i,j) = f(i,j) - (d2phidx2_ij + d2phidy2_ij);
        end
    end

end