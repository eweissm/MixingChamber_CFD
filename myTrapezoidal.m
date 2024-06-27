function I = myTrapezoidal(x, y)
    I = 0;
    for i = 2:length(x)
        %perform composite trapazoidal integration
        I = I + ( 0.5 * (x(i) - x(i-1)) * (y(i-1)+y(i)) );
    end

end