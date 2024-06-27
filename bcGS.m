function phi = bcGS(phi)

      %get size of phi
    M = size(phi,1)-2;
    N = size(phi,2)-2; 


    %apply zero neumann condition to left bc
    phi(:,1) = phi(:,2);

    %apply zero neumann condition to right bc
    phi(:,N+2) = phi(:, N+1);

    %apply zero neumann condition top bc
    phi(M+2,:) = phi(M+1,:);

    %apply zero neumann condition to bottom bc
    phi(1,:) = phi(2,:);

end