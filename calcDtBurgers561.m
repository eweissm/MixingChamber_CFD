function [dt, outputFlag] =  calcDtBurgers561(t,  outputTime, u, v)
    global h CFL
      
        %initialize output flag as 0
        outputFlag = 0;
        
        %calculate max stable dt step for u and v
        dt_max_u =h/(2*max(max(abs(u))) + max(max(abs(v))));
        dt_max_v =h/(max(max(abs(u))) + 2*max(max(abs(v))));

        dt_max = min(dt_max_u, dt_max_v);
        
        %apply security factor
        dt = dt_max*CFL;
        
        %check to see if we overshot our output time
        
        if (t < outputTime) && (t + dt >= outputTime)   
            outputFlag = 1;
            dt = outputTime-t;
        end

end