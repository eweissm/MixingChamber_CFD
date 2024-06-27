function [u, v] = correctOutletAlt(u,v)
    global h xc ucorr;
    
    %from problem definition --> outlet on bottom face for 1<x<2.5
    OutletXLimits = [1, 2.5];
    
    %get M and N
    M = size(v,1)-2; % in x dir
    N = size(v,2)-1;  %y dir

    % calculate q_dot* 
    q_dot_star = sum(u(1,2:N+1)*h) - sum(u(M+1,2:N+1)*h) + sum(v(2:M+1,1)*h) - sum(v(2:M+1,N+1)*h);
    
    %calculate u_corr
    ucorr = q_dot_star/ (OutletXLimits(2) - OutletXLimits(1));
    
    %correction will be applied to v at the outlet bc
    for i = 2:M+1
        if(xc(i)>OutletXLimits(1) && xc(i)<OutletXLimits(2))
            v(i,1) = v(i,1)-ucorr;
        end
    end

end
