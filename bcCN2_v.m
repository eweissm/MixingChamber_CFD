function [a, b, c, d] = bcCN2_v(a, b, c, d, t)

    global xc;

    %get size of v
    M = size(a,1)-2; %x dir
    N = size(a,2)-1; %y dir

    %for each col
    for i = 2:M+1

        %apply bottom BC
        if((xc(i) >= 1.5) && (xc(i) <= 2) ) %inlet condition
            z2 = 2*(xc(i)-1.5);
            g2 = sqrt(3)*6*z2*(1-z2);
                
            d(i,2) = d(i,2) - g2* a(i,2);
            a(i,2) = 0;

            %b and c remains unchanged
        else %wall condition
            a(i,2) = 0;
            %b, c and d remain unchanged
        end
        

        %apply top BC
        if((xc (i) >= 2) && (xc(i) <= 2.5) ) %inlet condition
            z3 = 2*(xc(i)-2);
            g3 = -9*z3*(1-z3);

            d(i,N) = d(i,N) - g3*c(i,N);
            c(i,N) = 0;
            %a and b remains unchanged

        else %wall condition
            c(i,N) = 0;
            %b, a and d remain unchanged
        end
    end
end