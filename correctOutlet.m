function [u, v] = correctOutlet(u,v)
    global h yc ucorr;
    
    %from problem definition --> outlet on bottom face for 1<x<2.5
    OutletYLimits = [.25, .75];
    
    %get M and N from v
    M = size(v,1)-2; % in x dir
    N = size(v,2)-1;  %y dir

    % calculate q_dot* 
    q_dot_star = h*(sum(u(1,2:N+1)) - sum(u(M+1,2:N+1)) + sum(v(2:M+1,1)) - sum(v(2:M+1,N+1)));
    
    %calculate u_corr
    ucorr = q_dot_star/ (OutletYLimits(2) - OutletYLimits(1));
    
    %correction will be applied to u at the outlet bc
    for j = 2:N+1
        if(yc(j)>OutletYLimits(1) && yc(j)<OutletYLimits(2))
            u(1,j) = u(1,j)-ucorr;
        end
    end

end
