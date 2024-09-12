function varargout = GlobalModel_piecewise_function(data1,data2,I_full,SOC_full,t_full,NoSOC,Nocurrent,static,endglobalmodel)
    
    fieldNames1 = fieldnames(data1);
    fieldNames2 = fieldnames(data2);
    for j = 1:length(fieldNames1)
        fieldName1 = fieldNames1{j};
        eval([fieldName1 ' = data1.' fieldName1 ';']);
    end
    for j = 1:length(fieldNames2)
        fieldName2 = fieldNames2{j};
        eval([fieldName2 ' = data2.' fieldName2 ';']);
    end
    % if ~pybamm
    %     run param/params_LCO
    %     [cn_low,cp_low] = init_cs(p,p.volt_min);
    %     [cn_high,cp_high] = init_cs(p,p.volt_max);
    %     Delta_cn = cn_high-cn_low;
    %     Delta_cp = cp_low-cp_high;
    %     p.OneC = min(p.epsilon_s_n*p.L_n*Delta_cn*p.Faraday/3600, p.epsilon_s_p*p.L_p*Delta_cp*p.Faraday/3600);
    % end
    %load('Parameters1RC.mat')
    %load('OCV_SOC.mat')
    
    if static
        nonempty_idx = find(~cellfun(@isempty, R0_full_c));
        R0c=cell2mat(R0_full_c(nonempty_idx));
        R0d=cell2mat(R0_full_d(nonempty_idx));
        R1c=R1c(nonempty_idx);
        R1d=R1d(nonempty_idx);
        C1c=C1c(nonempty_idx);
        C1d=C1d(nonempty_idx);
    else
    
        for k=1:length(i_ref)
            i_ref_c{k} = i_ref{k}(i_ref{k} > 0);
            temp = i_ref{k}(i_ref{k} <= 0);
            i_ref_d{k} = temp(2:end);
            [i_ref_c_unique{k},idx_c{k}]=unique(i_ref_c{k}','rows');     
            idx_c{k}=idx_c{k}'+1; %Bit of a hack to account for an error from setting v_measured(1)=0 in Test_identification
            [i_ref_d_unique{k},idx_d{k}]=unique(i_ref_d{k}','rows'); 
            idx_d{k}=idx_d{k}';
            if ~Nocurrent
                R0_full_c_u{k}=R0_full_c{k}(idx_c{k});
                R0_full_d_u{k}=R0_full_d{k}(idx_d{k});
            end
        end
        i_ref_c_unique=i_ref_c_unique';
        i_ref_d_unique=i_ref_d_unique';
        if ~Nocurrent
            R0_full_c_u=R0_full_c_u';
            R0_full_d_u=R0_full_d_u';
        end
    
        if NoSOC %if we're not using SOC as a parameter
            nonempty_idx = find(~cellfun(@isempty, R0_full_c_u));
            for k=1:length(R0_full_c_u)
                R0_full_c_u{k}=R0_full_c_u{nonempty_idx};
                R0_full_d_u{k}=R0_full_d_u{nonempty_idx};
                i_ref_c_unique{k}=i_ref_c_unique{nonempty_idx};
                i_ref_d_unique{k}=i_ref_d_unique{nonempty_idx};
                C1c(k)=C1c(nonempty_idx);
                C1d(k)=C1d(nonempty_idx);
                R1c(k)=R1c(nonempty_idx);
                R1d(k)=R1d(nonempty_idx);
            end
        end
        parameters=[R1c,R1d,C1c,C1d];%,OCV_SOC(1,:)'];
    
        
        allValues_c = vertcat(i_ref_c_unique{:}); % Concatenate all cell array elements into one vector
        allValues_d = vertcat(i_ref_d_unique{:});
        i_ref_c_unique2 = [0; unique(allValues_c)]; % Extract unique values and sort them. Add 0 as it's not given in data
        i_ref_d_unique2 = unique(allValues_d);
        i_ref_c_matrix = repmat(i_ref_c_unique2', length(i_ref_c_unique), 1); % Create a matrix by replicating i_ref_unique
        i_ref_d_matrix = repmat(i_ref_d_unique2', length(i_ref_d_unique), 1);
        R0_full_c_u_matrix = zeros(size(i_ref_c_matrix));
        R0_full_d_u_matrix = zeros(size(i_ref_d_matrix));
        
        if Nocurrent
            R0_full_c_array=cell2mat(R0_full_c);
            R0_full_d_array=cell2mat(R0_full_d);
            R0_full_c_u_matrix=repmat(R0_full_c_array,1,size(R0_full_c_u_matrix,2));
            R0_full_d_u_matrix=repmat(R0_full_d_array,1,size(R0_full_d_u_matrix,2));
        end
        if ~Nocurrent 
            for i = 1:length(R0_full_c_u)
                % Find the indices in i_ref_unique that match the currents in i_ref_c_unique{i}
                [~, indices_c] = ismember(i_ref_c_unique{i}, i_ref_c_unique2);
                [~, indices_d] = ismember(i_ref_d_unique{i}, i_ref_d_unique2);
                % Place the K_I_c_u values in the corresponding locations
                R0_full_c_u_matrix(i, indices_c) = R0_full_c_u{i};
                R0_full_d_u_matrix(i, indices_d) = R0_full_d_u{i};
                % Interpolate or extrapolate to fill in missing values
                R0_full_c_u_matrix(i, :) = interp1(i_ref_c_unique{i}, R0_full_c_u{i}, i_ref_c_unique2, 'linear', 'extrap');
                R0_full_d_u_matrix(i, :) = interp1(i_ref_d_unique{i}, R0_full_d_u{i}, i_ref_d_unique2, 'linear', 'extrap');
            end
        end
    end

    % for plotting (3D plots)
    if ~static
        ParamStructc=struct('current',i_ref_c_unique2,'SOC',HPPC_startingSOC,'R0',R0_full_c_u_matrix,'tau', tau1c);
        ParamStructd=struct('current',i_ref_d_unique2,'SOC',HPPC_startingSOC,'R0',R0_full_d_u_matrix,'tau', tau1d);
        varargout{2}=ParamStructc;
        varargout{3}=ParamStructd;
    end

    %%
    all_Volt_ECM=[];
    chunk_size = 5000;
    num_points = length(I_full);
    num_chunks = ceil(num_points / chunk_size);
    %num_chunks=5; %ATT
    V_ECM = zeros(1, num_points);  % Pre-allocate for speed
    i_initialcondition=0;
    x0 = 0; % Initial state for current x
    if exist('endglobalmodel', 'var') == 1
        num_chunks=endglobalmodel;
    end
    for k = 1:num_chunks
        
        start_idx = (k-1) * chunk_size + 1;
        end_idx = min(k * chunk_size, num_points);
        i=I_full(start_idx:end_idx);
        % V=V_full(start_idx:end_idx);
        SOC=SOC_full(start_idx:end_idx);
        t=t_full(start_idx:end_idx);
        % Calculation
        i(1)=i_initialcondition;
        i_initialcondition=i(end);
        if ~static
            interpolated_parameters = zeros(length(i), size(parameters, 2));
            for param = 1:size(parameters, 2)
                interpolated_parameters(:,param) = interp1(HPPC_startingSOC, parameters(:, param), SOC, 'spline', 'extrap');
            end
            R_c_interp = zeros(length(i), 1); 
            R_d_interp = zeros(length(i), 1);
            SOC = min(max(SOC, 0), 1); %
            i = min(max(i, min(i_ref_d_matrix,[],'all')), max(i_ref_c_matrix,[],'all')); %
    
            R_c_interp = interp2(i_ref_c_unique2, HPPC_startingSOC, R0_full_c_u_matrix,-i,SOC,"spline");
            R_d_interp = interp2(i_ref_d_unique2, HPPC_startingSOC, R0_full_d_u_matrix,-i,SOC,'spline');
                 
            R1c = interpolated_parameters(:,1)';
            R1d = interpolated_parameters(:,2)';
            C1c = interpolated_parameters(:,3)';
            C1d = interpolated_parameters(:,4)';
        end
        SOC(SOC < 0) = 0; 
        SOC(SOC > 1) = 1; 
        startingOCV_interpol = interp1(OCV_SOC(2,:), OCV_SOC(1,:), SOC, 'linear', 'extrap');
        tau1=R1d.*C1d; 
        if static
            tau1=tau1+zeros(size(t));
        end
        charging_cond=i<0;
        discharging_cond=i>=0;

        % Solve the ODE for x: (x=[iR1])
        x = ode4(@(t_ode, x) ecm_model1RC_global(t,t_ode, x, interp1(t, i, t_ode), tau1), t, x0);
        if static
            v_model_charging = R0c* i(charging_cond) + R1c* x(charging_cond)';
            v_model_discharging = R0d* i(discharging_cond) + R1d* x(discharging_cond)';
        else
            v_model_charging = R_c_interp(charging_cond).* i(charging_cond) + R1c(charging_cond) .* x(charging_cond)';
            v_model_discharging = R_d_interp(discharging_cond).* i(discharging_cond) + R1d(discharging_cond) .* x(discharging_cond)';
        end
        d=1; % Merge charge & discharge
        c=1;
        v_model_total=zeros(1,size(i,2));
        for idx= 1:length(i)
            if charging_cond(idx)==1
                v_model_total(idx)=-v_model_charging(c);
                c=c+1;
            elseif discharging_cond(idx)==1
                v_model_total(idx)=-v_model_discharging(d);
                d=d+1;
            end
        end
        v_model_total=v_model_total+startingOCV_interpol;
        x0=x(end,:); %update initial conditions
        all_Volt_ECM=[all_Volt_ECM, v_model_total];
    end
    V_ECM=all_Volt_ECM;
    varargout{1}=V_ECM;