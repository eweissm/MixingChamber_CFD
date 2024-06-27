function f =  calc_adYdx_WENO(Y, a)

    global h;

    %get size of Y
    M = size(Y,1)-6;

    %declare f
    f = zeros( M+6,1);
    dp = zeros( M+6,1);
    dmdp = zeros( M+6,1);
    
    % get dp and dmdp
    dp(1:M+5) = Y(2:M+6) - Y(1:M+5);
    dmdp(2:M+5) = Y(3:M+6) - 2*Y(2:M+5) + Y(1:M+4);

    
    %calculate f
    for i = 4:M+3
        
        %get WENO value --> will be added to dphidx
        if(a(i)>=0)
            WENO_val = psiWENO(dmdp(i-2)/h, dmdp(i-1)/h,dmdp(i)/h,dmdp(i+1)/h);
        else
            WENO_val = -psiWENO(dmdp(i+2)/h, dmdp(i+1)/h,dmdp(i)/h,dmdp(i-1)/h);
        end

        dphidx = (1/(12*h)) *(-dp(i-2)+ 7*dp(i-1) +7*dp(i) - dp(i+1)) - WENO_val;

        f(i) = a(i)*dphidx;

    end
    

end