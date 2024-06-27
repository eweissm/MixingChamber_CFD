function Sol = mySolveTriDiag(a, b, c, d)
    % This functions uses gaussian elimination to solve tri-diagonal matices 
    % See Slides 42 to 60 of module 2 to see detailed derivation of method
    
    
    % Check that all vectors are same length by attempting to concatonate all
    % vectors --> this will throw an error if vectors are not equal lengths
    try
        [a b c d];
    catch 
        error("Vectors must all be same Length")
    end
    
    %find length of p
    p = length(a);
    
    %perform first elimination step
    for i = 2:p
        b(i) = b(i) - c(i-1)*a(i)/b(i-1);
        d(i) = d(i) - d(i-1)*a(i)/b(i-1);
    end
    
    %perform back substitution step
    d(p) = d(p)/b(p);
    for i = (p-1):-1:1
        d(i) = (d(i) - c(i)*d(i+1))/b(i);
    end
    
    Sol = d;
end