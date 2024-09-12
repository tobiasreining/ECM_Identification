%% Preamble
clc
clear
close all
fs=20;
lw=2;
run param/params_LCO
addnoise=true
[cn_low,cp_low] = init_cs(p,p.volt_min);
[cn_high,cp_high] = init_cs(p,p.volt_max);
Delta_cn = cn_high-cn_low;
Delta_cp = cp_low-cp_high;
p.OneC = min(p.epsilon_s_n*p.L_n*Delta_cn*p.Faraday/3600, p.epsilon_s_p*p.L_p*Delta_cp*p.Faraday/3600);
OneC=p.OneC; %potentially overwritten in the following line
load('HPPC_voltage_adjusted.mat')
%HPPC_SOC_levels=[99,98,97,96,95,90,85,80,75,70,65,60,55,50,45,40,35,30,25,20,15,10,5,4,3,2,1]/100;
HPPC_SOC_levels=[97.5:-2.5:2.5]/100;

specific_SOC = [99, 98, 97, 3, 2, 1]/100; % Divided by 100 to match the scale of HPPC_SOC_levels
C_rates_standard=[0.05,0.1,0.15,0.2,0.5,1,1.5,2];
C_rates_extreme={[0.05,0.1,0.15,0.2,0.5,1,1.4,1.5,2],[0.05,0.1,0.15,0.2,0.5,1,1.5,1.6,2],[0.05,0.1,0.15,0.2,0.5,1,1.5,1.8,2],[0.05,0.1,0.15,0.2,0.5,1,1.5,1.8,2],[0.05,0.1,0.15,0.2,0.5,1,1.1,1.5,2],[0.05,0.1,0.15,0.2,0.5,1,1.5,2]};
C_rates=cell(length(HPPC_SOC_levels),1);
I_points=cell(length(HPPC_SOC_levels),1);
N_total=cell(length(HPPC_SOC_levels),1);
k_ex=1;
for k=1:length(HPPC_SOC_levels)
    if ismember(HPPC_SOC_levels(k), specific_SOC) 
        C_rates{k}=C_rates_extreme{k_ex};
        k_ex=k_ex+1;
    else
        C_rates{k}=C_rates_standard;
    end
    I_points{k}=C_rates{k}*OneC;
    N_total{k} = length(I_points{k});
end

% %C_rates=[0.5,1,1.5,2];
% HPPC_current=[];
% HPPC_cycle=[ones(1,10)*-OneC,zeros(1,40),ones(1,10)*OneC,zeros(1,40)];
% for C_rate = C_rates
%     HPPC_current=[HPPC_current,HPPC_cycle*C_rate];
% end
% i=-HPPC_current;
% I_points=C_rates*OneC;
%% Using actual data from the DFN model



time=cell(size(voltage_adjusted,2),1);
for k = 1:length(voltage_adjusted)
    if ~isempty(voltage_adjusted{k}) && length(voltage_adjusted{k}) > 1
        voltage_adjusted{k}(:, 1:end-1) = voltage_adjusted{k}(:, 2:end);
        voltage_adjusted{k}(end)=voltage_adjusted{k}(end-1);
        time{k}=1:length(HPPC_current{k});
    end
end

%%
% create voltage_adjusted_avg to test influence of SOC parameter
valid_cells={};
reference_position = ceil(length(voltage_adjusted)/2);
reference_length = length(voltage_adjusted{reference_position});
for k=1:length(voltage_adjusted)
    if length(voltage_adjusted{k})==reference_length %just use the standard cycles, not the extreme SOC cycles (different length)
        valid_cells{end+1} = voltage_adjusted{k};
    end
end

if ~isempty(valid_cells)
    voltage_adjusted_avg = zeros(1, reference_length); % Initialize the average array
    for t = 1:reference_length
        sum_voltage = 0;
        for k = 1:length(valid_cells)
            sum_voltage = sum_voltage + valid_cells{k}(t);
        end
        voltage_adjusted_avg(t) = sum_voltage / length(valid_cells);
    end
end

voltage_adjusted{reference_position}=voltage_adjusted_avg; %just a quick fix to adhere to the cell format. stuffing the avg into the middle, not very elegant

%%
if addnoise
    figure
    sigma=0.1; %5% noise
    hold on
    
    voltage_noisy = cell(size(voltage_adjusted));
    for i = 10%1:numel(voltage_adjusted)
        current_vector = voltage_adjusted{i};
        proportional_sigma = sigma * abs(current_vector);

        noise = proportional_sigma .* randn(size(current_vector));    

        voltage_noisy{i} = current_vector + noise;
        hold on
        plot(1:length(voltage_noisy{i}),voltage_noisy{i},'LineWidth', lw)
        plot(1:length(voltage_adjusted{i}),voltage_adjusted{i},'LineWidth', lw,'LineStyle','-')

        %legend('DFN with added noise','DFN', 'interpreter', 'latex', 'Location', 'southwest','FontSize', fs);
        ylabel('$\Delta$ Voltage (V)','interpreter', 'latex','FontSize', fs)
        xlabel('Time (s)','interpreter', 'latex','FontSize', fs)
        ax = gca; % Get current axes
        ax.FontSize = fs; % Set font size for tick markers
        %xlim([390,800])
        legend({'DFN + noise','DFN'},'Interpreter','latex','Location','southwest')
        box on
        ax.FontSize = fs;
        set(1, 'units', 'centimeters', 'pos', [0 0 35 12])
        exportgraphics(ax,'HPPC noise_V2.pdf')

    end
end

%%
%voltage_adjusted(:,1:end-1)=voltage_adjusted(:,2:end); %just to shift the voltage by one element
%voltage_adjusted(:,end)=0; %and the last element is basically 0
%load('ratios.mat'); %voltage ratios and I points for Hammerstein Model
%N = length(I_points);


%v_model_total=zeros(size(voltage_adjusted));
v_model_total=cell(size(voltage_adjusted,2,1));
%for charging, for now not discharging:
R1c=zeros(size(voltage_adjusted,2),1);
C1c=zeros(size(voltage_adjusted,2),1);
R0d=zeros(size(voltage_adjusted,2),1);
R1d=zeros(size(voltage_adjusted,2),1);
C1d=zeros(size(voltage_adjusted,2),1);
i_ref=cell(size(voltage_adjusted,2),1);
R0_full_c=cell(size(voltage_adjusted,2),1);
R0_full_d=cell(size(voltage_adjusted,2),1);
RMSE=zeros(size(voltage_adjusted,2),1);
charging_cond=cell(size(voltage_adjusted,2),1);
discharging_cond=cell(size(voltage_adjusted,2),1);
discharging_cond_nonzero=cell(size(voltage_adjusted,2),1);



called_cycles=[]; %so that indiviual cycles can also be plotted
if addnoise
    voltage_adjusted=voltage_noisy;
end
for cycle = 1:size(voltage_adjusted,2) %ATT
    N=N_total{cycle};
    i=-HPPC_current{cycle};
    charging_cond{cycle}=i>0;
    discharging_cond{cycle}=i<=0;
    discharging_cond_nonzero{cycle}=i<0;
    i(:,1)=0; %first current element should be 0
    called_cycles=[called_cycles,cycle]; 
    v_measured=voltage_adjusted{cycle};
    t=1:length(v_measured);
    %a quick fix to account for lower currents at extreme SOCs
    % if v_measured(2)<0.005  %a QUICK fix
    %     i=i*1/10;
    % end

    v_measured(1)=0; %This leads to a bit of an issue with the K_I. K_I(1) has to be set to K_I(2)
    tau1_0_nonn = 20; %nonn: non-normalized
    spline_par_0_nonn = zeros([1,N])+1;%v_ratio*0.7; %spline parameters. guesses
    % Set minimum and maximum bounds for tau
    tau_limit=[10,50];
    spline_par_limit=[0.1,5];
    % Normalized tau initial guess
    tau1_0 = (tau1_0_nonn - tau_limit(1)) / (tau_limit(2) - tau_limit(1));
    spline_par_0=(spline_par_0_nonn-spline_par_limit(1))/(spline_par_limit(2)-spline_par_limit(1));
    params_0=[spline_par_0, tau1_0];
    %lower, upper bound
    lb = [0*ones(N, 1); 0]; % lower bounds 
    ub = [1*ones(N, 1); 1]; % upper bounds 
    % Solve the optimization problem. lsqnonlin squares and minimizes the
    % error. "cost_function" just returns the residual 
    options = optimoptions('lsqnonlin', 'Display','iter','Algorithm', 'levenberg-marquardt','ScaleProblem','jacobian', 'MaxFunctionEvaluations', 1e4, 'MaxIterations',1000);
    [params_estimated, resnorm,residual,exitflag,output] = lsqnonlin(@(params) cost_function(params, t, i, v_measured,tau_limit,spline_par_limit,I_points{cycle}),params_0, lb, ub,[], [],[],[],[], options);
    spline_par_estimated_nonn = params_estimated(1:N)*(spline_par_limit(2)-spline_par_limit(1))+spline_par_limit(1);
    params_estimated_nonn = params_estimated(N+1) * (tau_limit(2) - tau_limit(1)) + tau_limit(1);
    
    fprintf('Estimated parameters: tau1 = %f\n', params_estimated_nonn(1));

    tau1 = params_estimated_nonn(1);
    K_I = interp1(I_points{cycle}, spline_par_estimated_nonn, abs(i), 'linear', 'extrap');
    K_I_charging=K_I(charging_cond{cycle});
    K_I_charging(1)=K_I_charging(2); %This is to correct an error from setting v_measured(1)=0 which conflicts with charging_cond.
    K_I_discharging=K_I(discharging_cond{cycle});

    x0 = 0; % Initial state for current x
    % Solve the ODE for x: (x=[iR1])
    x = ode4(@(t_ode, x) ecm_model(t, x, interp1(t, i, t_ode), tau1), t, x0);
    I_matrix_charging = [i(charging_cond{cycle})', x(charging_cond{cycle})];
    I_matrix_discharging = [i(discharging_cond{cycle})', x(discharging_cond{cycle})];
    R_charging = lsqnonneg(I_matrix_charging, v_measured(charging_cond{cycle})');
    R_discharging = lsqnonneg(I_matrix_discharging, v_measured(discharging_cond{cycle})');
    R0_charging = R_charging(1);
    %R0c(cycle)=R0_charging;
    R0_discharging = R_discharging(1);
    %R0d(cycle)=R0_discharging;
    R1_charging=R_charging(2);
    R1c(cycle)=R1_charging;
    R1_discharging=R_discharging(2);
    R1d(cycle)=R1_discharging;
    C1_charging=tau1/R1_charging;
    C1c(cycle)=C1_charging;
    C1_discharging=tau1/R1_discharging;
    C1d(cycle)=C1_discharging;
    i_ref{cycle}=i;
    
    R0_full_charging=R0_charging*K_I_charging;
    R0_full_discharging=R0_discharging*K_I_discharging;
    R0_full_c{cycle}=R0_full_charging;
    R0_full_d{cycle}=R0_full_discharging; 

    C1_charging=tau1/R1_charging;
    C1_discharging=tau1/R1_discharging;
    fprintf('Cycle %d: Estimated parameters: R0_charging=%f, R0_discharging=%f, R1_charging=%f, R1_discharging=%f,C1_charging=%f,C1_discharging=%f\n', cycle, R0_charging, R0_discharging, R1_charging,R1_discharging, C1_charging,C1_discharging);
    
    v_model_charging = R0_full_charging .* i(charging_cond{cycle}) + R1_charging * x(charging_cond{cycle})';
    v_model_discharging = R0_full_discharging .* i(discharging_cond{cycle}) + R1_discharging * x(discharging_cond{cycle})';
    


    d=1;
    c=1;
    v_model_total_cycle=[];
    for idx= 1:length(v_measured)
        if charging_cond{cycle}(idx)==1
            v_model_total_cycle(idx)=v_model_charging(c);
            c=c+1;
        elseif discharging_cond{cycle}(idx)==1
            v_model_total_cycle(idx)=v_model_discharging(d);
            d=d+1;
        end
    end
    v_model_total{cycle}=v_model_total_cycle;
    MAE = mean(abs(v_measured - v_model_total{cycle}));
    MSE = mean((v_measured - v_model_total{cycle}).^2);
    RMSE(cycle) = sqrt(MSE);
    fprintf('MAE: %.4f, MSE: %.4f, RMSE: %.4f', MAE, MSE, RMSE(cycle));
end

legendEntries1 = cell(1, cycle); 
figure()
hold on
for idx = called_cycles
    plot(time{idx}(2:end),voltage_adjusted{idx}(2:end)-v_model_total{idx}(2:end))
    legendEntries1{idx}=sprintf('SOC %.2f , mean: %.6f', HPPC_startingSOC(idx), mean(voltage_adjusted{idx}(2:end) - v_model_total{idx}(2:end)));
end
legendEntries1 = legendEntries1(~cellfun('isempty', legendEntries1));
legend(legendEntries1, 'interpreter', 'latex', 'Location', 'southwest');
hold off


figure()
ax1=subplot(2,1,1);
ylabel('Voltage','FontSize', fs)
xlabel('Time [sec]','FontSize', fs)


set(gca,'FontSize', fs)
xlim([1, 950])
hold on
colors = lines(cycle);
legendEntries2 = cell(1, 2*cycle);  
% for idx = called_cycles
%     plot(time{idx}, voltage_adjusted{idx}, ':', 'LineWidth', 2, 'Color', colors(idx, :));
%     plot(time{idx}, v_model_total{idx}, '-', 'LineWidth', 2, 'Color', colors(idx, :));
%     legendEntries2{2*idx-1} = sprintf('SOC %.2f (DFN)', HPPC_startingSOC(idx));
%     legendEntries2{2*idx} = sprintf('SOC %.2f (ECM)', HPPC_startingSOC(idx));    
% end
legendEntries2 = legendEntries2(~cellfun('isempty', legendEntries2)); %if individual cycles are plotted

legend(legendEntries2, 'interpreter', 'latex', 'Location', 'southwest', 'FontSize', 6);

hold off
ax2=subplot(2,1,2);
plot(time{5}, -HPPC_current{5},'LineWidth', 2);
%legend({'$$I(t)$$'},'interpreter','latex')
ylabel('Standard Current','FontSize', fs)
xlabel('Standard Time [sec]','FontSize', fs)
linkaxes([ax1,ax2],'x')

figure()
plot(HPPC_startingSOC, R1d,'LineWidth', 2);
% legend('R1c','R1d');
ylabel('R1 [$\Omega$]','FontSize', fs,'Interpreter','latex')
title('R1 Chen2020','FontSize',fs,'Interpreter','latex')
xlabel('SOC','FontSize', fs, 'Interpreter','latex')
hold on

%%
%close all
figure()
legendEntries3 = cell(1, cycle); 
%for idx=1:size(R0_full_c,1)
for idx=[23,24,25,26]%
    plot(abs(-HPPC_current{idx}(charging_cond{idx})),R0_full_c{idx},'o-')
    hold on
    legendEntries3{idx}=sprintf('SOC %.2f', HPPC_startingSOC(idx));
end
legendEntries3 = legendEntries3(~cellfun('isempty', legendEntries3));
legend(legendEntries3, 'interpreter', 'latex', 'Location', 'southwest');
ylabel('Resistance','FontSize', fs)
xlabel('abs(Current)','FontSize', fs)
title('R0 over current and SOC - charging')

figure()
legendEntries3 = cell(1, cycle); 
%for idx=1:size(R0_full_c,1)

for idx=[1:5,7,9]%
    dischargeHPPC=-HPPC_current{idx}(discharging_cond{idx});
    dischargeR0=R0_full_d{idx};
    idx_nonzero=find(dischargeHPPC<0);
    
    plot(dischargeHPPC(idx_nonzero), R0_full_d{idx}(idx_nonzero))
    %plot(-HPPC_current{idx}(discharging_cond{idx}),R0_full_d{idx})
    hold on
    legendEntries3{idx}=sprintf('SOC %.2f', HPPC_startingSOC(idx));
end
legendEntries3 = legendEntries3(~cellfun('isempty', legendEntries3));
legend(legendEntries3, 'interpreter', 'latex', 'Location', 'southwest');
ylabel('Resistance','FontSize', fs)
xlabel('abs(Current)','FontSize', fs)
title('R0 over current and SOC - discharging')



%%
figure()
plot(HPPC_startingSOC,RMSE,'o-');
xlabel('SOC')
ylabel('RMSE')
title('RMSE over SOC')

%% padding for interp2 to work?
R1c = pad_array(R1c);
R1d = pad_array(R1d);
C1c = pad_array(C1c);
C1d = pad_array(C1d);
HPPC_startingSOC = pad_array(HPPC_startingSOC);
HPPC_startingSOC(1)=1;
HPPC_startingSOC(end)=0;
HPPC_startingOCV = pad_array(HPPC_startingOCV);
HPPC_startingOCV(1)=4.48;  %ATT hardcoding
HPPC_startingOCV(end)=3.30;
%i_ref_pad=cell(length(i_ref)+2,1);
HPPC_current = [{[HPPC_current{1}]}; HPPC_current; {[HPPC_current{end}]}];
i_ref = [{[i_ref{1}]}; i_ref; {[i_ref{end}]}];  % i_ref_pad{2:end-1}=i_ref;
R0_full_c = [{[R0_full_c{1}]}; R0_full_c; {[R0_full_c{end}]}]; 
R0_full_d = [{[R0_full_d{1}]}; R0_full_d; {[R0_full_d{end}]}]; 
charging_cond = [{[charging_cond{1}]}; charging_cond; {[charging_cond{end}]}]; 
discharging_cond = [{[discharging_cond{1}]}; discharging_cond; {[discharging_cond{end}]}]; 

tau1c=R1c.*C1c;
tau1d=R1d.*C1d;

%% save
save('Parameters_pybamm_Marquis_noise5.mat', 'OneC','R0_full_c','R0_full_d','R1c','R1d','C1c','C1d','i_ref','tau1c','tau1d','HPPC_startingOCV',"HPPC_startingSOC");

% %Lookup Table Business
% I_points=C_rates*OneC;
% v_indices=[10,110,210,310,410,510,610,710];
% v_measured_points=v_measured(v_indices);
% v_approx_points=v_model_total(cycle,v_indices); %ATT
% v_ratio=v_measured_points./v_approx_points; %Forget the v_ratio thing
% I_range=[0:1:ceil(max(I_points))];
% v_ratio_interp = interp1(I_points, v_ratio, I_range, 'linear', 'extrap');
% figure();
% plot(I_range,v_ratio_interp);
% save('ratios.mat', 'I_points', 'v_ratio');


