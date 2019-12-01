
function output = gpops_Endpoint(input)

    cost  = input.phase.integral;
    output.objective = cost;

    phase_num = input.auxdata.phase_num;
    order = input.auxdata.order;
    
    t0 = zeros(phase_num,1);
    tf = zeros(phase_num,1);
    x0 = zeros(phase_num,order);
    xf = zeros(phase_num,order);
    
    for i = 1:phase_num-1

        t0(i) = input.phase(i).initialtime;
        tf(i) = input.phase(i).finaltime;
        x0(i,:) = input.phase(i).initialstate;
        xf(i,:) = input.phase(i).finalstate;
        
    end
    
    for i = 1:phase_num-1
        
        output.eventgroup(i).event = [x0(i+1,:)-xf(i,:), t0(i+1)-tf(i)];

    end
    
    
    
end
