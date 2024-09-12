function dx = ecm_model1RC_global(t,t_ode, x, i, tau1)  %x is the current  
    if mod(t_ode,1000)==0
        t_ode
    end
    current_tau1 = interp1(t,tau1,t_ode);
    dx = [-1/(current_tau1) * x(1) + 1/(current_tau1) * i];
end