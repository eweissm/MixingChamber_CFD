function eh = myProlong(e2h)
%Prolongation function for cell centered mesh

    %get fine mesh dimensions
    M_2h = size(e2h,1)-2;
    N_2h = size(e2h,2)-2;

    %determine courser mesh size
    M_h = M_2h*2;
    N_h = N_2h*2;

    %delcare eh with ghost cells
    eh = zeros(M_h+2, N_h+2);

    for j = 2:N_2h+1
        for i = 2:M_2h+1

            %apply prolongation scheme for cell centered mesh
            eh(2*i-2, 2*j-2) = e2h(i,j);   
            eh(2*i-1, 2*j-2) = e2h(i,j); 
            eh(2*i-2, 2*j-1) = e2h(i,j); 
            eh(2*i-1, 2*j-1) = e2h(i,j); 
        
        end
    end
    
    eh = bcGS(eh);

end