function [dt, outputFlag] = calcDt(t, outputTime, a)

    if exist('a', 'var') %if we pass an a argument, then calcDt if for a hyperbolics function

        global h CFL;
        
        %initialize output flag as 0
        outputFlag = 0;
        
        %calculate max stable dt step
        dt_max =h/max(abs(a));
        
        %apply security factor
        dt = dt_max*CFL;
        
        %check to see if we overshot our output time
        
        if (t < outputTime) && (t + dt >= outputTime)   
            outputFlag = 1;
            dt = outputTime-t;
        end

    else %otherwise its a parabolic
         global h CFL Re Sc;
    
        %initialize output flag as 0
        outputFlag = 0;
        
        %calculate max stable dt step
        dt_max =min(0.25*h^2*Re, 0.25*h^2*Re*Sc);
    
        %apply security factor
        dt = dt_max*CFL;
    
        %check to see if we overshot our output time
    
        if (t < outputTime) && (t + dt >= outputTime)   
            outputFlag = 1;
            dt = outputTime-t;
        end
    end
end