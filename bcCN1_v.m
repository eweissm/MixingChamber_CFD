function [a, b, c, d] = bcCN1_v(a, b, c, d, t)

    global yf;

    %get size of v
    M = size(a,1)-2; %x dir
    N = size(a,2)-1; %y dir

    %for each row
    for j = 2:N
        %apply left BC
        if((yf(j) >= 0.25) && (yf(j) <= 0.75) ) %outlet condition
            b(2, j) = b(2, j) + a(2, j);
            a(2, j) = 0;
            %d and c unchanged
        else %wall condition
            b(2, j) = b(2, j) - a(2, j);
            a(2, j) = 0;
            %d and c unchanged
        end
        
        %apply right BC
        b(M+1, j) = b(M+1, j) - c(M+1, j);
        c(M+1, j) = 0; 
        % d and a unchanged
       
    end
end