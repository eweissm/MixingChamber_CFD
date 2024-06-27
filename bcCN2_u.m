function [a, b, c, d] = bcCN2_u(a, b, c, d, t)

    global xf;

    %get size of u
    M = size(a,1)-1; %node based in x dir
    N = size(a,2)-2; %cell centered in y dir

    %for each col
    for i = 2:M

        %apply top BC
        if((xf(i) > 1.5) && (xf(i) < 2) ) %inlet condition
            z2 = 2*(xf(i)-1.5);
            g2 = 6*z2*(1-z2);
                
            d(i,2) = d(i,2) - 2*g2*a(i,2);
            b(i,2) = b(i,2) - a(i,2);
            a(i,2) = 0;
            %c remains unchanged
        else %wall condition
            b(i,2) = b(i,2) - a(i,2);
            a(i,2) = 0;
            %c and d remain unchanged
        end
        

        %apply bottom BC
        if((xf(i) > 2) && (xf(i) < 2.5) ) %inlet condition
            z3 = 2*(xf(i)-2);
            g3 = sqrt(3)*(3/2)*6*z3*(1-z3);

            d(i,N+1) = d(i,N+1) - 2*g3*c(i,N+1);
            b(i,N+1) = b(i,N+1) - c(i,N+1);
            c(i,N+1) = 0;
            %a remains unchanged

        else %wall condition
            b(i,N+1) = b(i,N+1) - c(i,N+1);
            c(i,N+1) = 0;
            %a and d remain unchanged
        end
    end
end