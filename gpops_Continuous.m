function phaseout = gpops_Continuous( input )

    phase_num = input.auxdata.phase_num;
    xf = input.auxdata.xf;
    dt = input.auxdata.dt;
    
    for j = 1:phase_num
        
        x_1 = input.phase(j).state(:,1);
        x_2 = input.phase(j).state(:,2);
        
        p = input.phase(j).parameter;
        
        u_1 = p(:,3*j-2);
        u_2 = p(:,3*j-1);
        u_3 = p(:,3*j);
        
        
        if j == 1
            u_1_next = 0;
            u_2_next = 0;
            u_3_next = 0;
        else
            u_1_next = p(:,3*j-5);
            u_2_next = p(:,3*j-4);
            u_3_next = p(:,3*j-3);
        
        end
        
        
        
        x = [x_1, x_2]';
        u = [u_1, u_2, u_3]';
        
        vf2 = vf_koopman(x,u);

        dynamics = vf2';

        phaseout(j).dynamics = dynamics;
        
        phaseout(j).integrand = 100* (x_1 - xf(1)).^2 + 100*(x_2 - xf(2)).^2 + 1*(u_1.^2 + u_2.^2 +u_3.^2);
        
        slope_u_1 = (u_1_next - u_1)/dt;
        slope_u_2 = (u_2_next - u_2)/dt;
        slope_u_3 = (u_3_next - u_3)/dt;
        
        phaseout(j).path = [slope_u_1,slope_u_2,slope_u_3];
    end
    
end


