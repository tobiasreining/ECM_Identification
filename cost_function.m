function error = cost_function(params, t, i, v_measured, tau_limit,spline_par_limit,I_points)
    charging_cond=i<0;
    discharging_cond=i>=0;
    N=length(I_points);
    lambda_error=0.001;
    spline_par=params(1:N)*(spline_par_limit(2)-spline_par_limit(1))+spline_par_limit(1);
    tau1 = params(N+1) * (tau_limit(2) - tau_limit(1)) + tau_limit(1);
    K_I = interp1(I_points, spline_par, abs(i), 'linear', 'extrap');
    K_I_charging=K_I(charging_cond);
    K_I_discharging=K_I(discharging_cond);
    x0 = [0]; % Initial state for current x
    % Solve the ODE for x: (x=[iR1])
    x = ode4(@(t_ode, x) ecm_model(t_ode, x, interp1(t, i, t_ode), tau1), t, x0);
    
    I_matrix_charging = [i(charging_cond)', x(charging_cond)];
    I_matrix_discharging = [i(discharging_cond)', x(discharging_cond)];
    R_charging = lsqnonneg(I_matrix_charging, v_measured(charging_cond)');
    R_discharging = lsqnonneg(I_matrix_discharging, v_measured(discharging_cond)');
    R0_charging = R_charging(1);
    R0_discharging = R_discharging(1);
    R1_charging=R_charging(2);
    R1_discharging=R_discharging(2);
    v_model_charging = K_I_charging.*(R0_charging * i(charging_cond)) + R1_charging * x(charging_cond)';
    v_model_discharging = K_I_discharging.*(R0_discharging * i(discharging_cond)) + R1_discharging * x(discharging_cond)';

    error = zeros(1,size(v_measured,2)+1);
    error(charging_cond) = v_measured(charging_cond) - v_model_charging;
    error(discharging_cond) = v_measured(discharging_cond) - v_model_discharging;
    error(end)=lambda_error*tau1;
end