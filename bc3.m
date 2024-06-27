function [Y] = bc3(Y, t)

    %get size of Y
    M = size(Y,1)-6;

    %apply left BC
    

%for questions 10
% g = (21+5*sqrt(2)+20*cos(5*t*pi))/50;
    % Y(1) = 2*g - Y(6);
    % Y(2) = 2*g - Y(5);
    % Y(3) = 2*g - Y(4);
    % 
    % 
    % %apply right bc
    % Y(M+4) = Y(M+3);
    % Y(M+5) = Y(M+2);
    % Y(M+6) = Y(M+1);

    %for questions 11 --periodic BC  
    
    Y(1) = Y(M+1);
    Y(2) = Y(M+2);
    Y(3) = Y(M+3);


    %apply right bc
    Y(M+4) = Y(4);
    Y(M+5) = Y(5);
    Y(M+6) = Y(6);

end