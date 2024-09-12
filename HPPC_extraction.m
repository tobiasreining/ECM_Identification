%% Extract HPPC data from results
clear
clc
run param/params_CHEN
[cn_low,cp_low] = init_cs(p,p.volt_min);
[cn_high,cp_high] = init_cs(p,p.volt_max);
Delta_cn = cn_high-cn_low;
Delta_cp = cp_low-cp_high;
p.OneC = min(p.epsilon_s_n*p.L_n*Delta_cn*p.Faraday/3600, p.epsilon_s_p*p.L_p*Delta_cp*p.Faraday/3600);

mat_files = dir('results_*.mat');
all_SOC = [];
all_I = [];
all_Volt = [];

for i = 1:length(mat_files)
    load(mat_files(i).name, 'SOC_save', 'I_save', 'Volt_save');
    if i==1
        I_save=I_save;
        SOC_save=SOC_save';
        Volt_save=Volt_save';
    end
    all_SOC = [all_SOC SOC_save];
    all_I = [all_I, I_save'];
    all_Volt = [all_Volt Volt_save];

end
start_HPPC=1001;
SOC=all_SOC(start_HPPC:end);
Volt=all_Volt(start_HPPC:end);
I=all_I(start_HPPC:end);
t=1:length(I);


%C_rates=[0.5,1,1.5,2];
C_rates=[0.05,0.1,0.15,0.2,0.5,1,1.5,2];

HPPC_current_standard=[];
HPPC_cycle=[ones(1,10)*-p.OneC,zeros(1,40),ones(1,10)*p.OneC,zeros(1,40)];
for C_rate = C_rates
    HPPC_current_standard=[HPPC_current_standard,HPPC_cycle*C_rate];
end
length_HPPC_standard=length(HPPC_current_standard);
HPPC_SOC_levels=[99,98,97,96,95,90,85,80,75,70,65,60,55,50,45,40,35,30,25,20,15,10,5,4,3,2,1]/100;
specific_SOC = [99, 98, 97, 3, 2, 1]/100; % Divided by 100 to match the scale of HPPC_SOC_levels
HPPC_current=HPPC_current_standard;

%% comment out if not extreme SOC

C_rates_extreme={[0.05,0.1,0.15,0.2,0.5,1,1.4,1.5,2],[0.05,0.1,0.15,0.2,0.5,1,1.5,1.6,2],[0.05,0.1,0.15,0.2,0.5,1,1.5,1.8,2],[0.05,0.1,0.15,0.2,0.5,1,1.5,1.8,2],[0.05,0.1,0.15,0.2,0.5,1,1.1,1.5,2],[0.05,0.1,0.15,0.2,0.5,1,1.5,2]};
C_ex=[1.4,1.6,1.8,1.8,1.1,0.5]; %Max C for 99,98,97,3,2,1 SOC
%just for CHEN
% C_rates_extreme={[0.05,0.1,0.15,0.2,0.5,1,1.4,1.5,2],[0.05,0.1,0.15,0.2,0.5,1,1.5,1.6,2],[0.05,0.1,0.15,0.2,0.5,1,1.5,1.8,2],[0.05,0.1,0.15,0.2,0.5,1,1.5,2],[0.05,0.1,0.15,0.2,0.5,1,1.5,2],[0.05,0.1,0.15,0.2,0.5,1,1.5,2]};
% C_ex=[1.4,1.6,1.8,1,0.5,0.2]; %Max C for 99,98,97,3,2,1 SOC
HPPC_current_extreme={};
for i = 1:length(C_rates_extreme)
    HPPC_current_extreme{i}=[];
    for C_rate = C_rates_extreme{i}
        if C_rate <= C_ex(i)
            HPPC_current_extreme{i}=[HPPC_current_extreme{i},HPPC_cycle*C_rate];
        elseif i <=  3 %High SOC
            HPPC_current_extreme{i}=[HPPC_current_extreme{i},ones(1,ceil((C_rate/C_ex(i))*10*2))*-p.OneC/2*C_ex(i),zeros(1,40),ones(1,10)*p.OneC*C_rate,zeros(1,40)];
        elseif i > 3 %Low SOC
            HPPC_current_extreme{i}=[HPPC_current_extreme{i},ones(1,10)*-p.OneC*C_rate,zeros(1,40),ones(1,ceil((C_rate/C_ex(i))*10*2))*p.OneC/2*C_ex(i),zeros(1,40)];    
        end
    end
end
%%
HPPC_current=cell(length(HPPC_SOC_levels),1);
k_ex=1;
for k=1:length(HPPC_SOC_levels)
    if ismember(HPPC_SOC_levels(k), specific_SOC) 
        HPPC_current{k}=HPPC_current_extreme{k_ex};
        k_ex=k_ex+1;
    else
        HPPC_current{k}=HPPC_current_standard;
    end
end
length_HPPC=zeros(1,length(HPPC_SOC_levels))+length_HPPC_standard;

length_HPPC=zeros(1,length(HPPC_SOC_levels));
for k = 1:length(HPPC_SOC_levels)
    if k <= 3 
        length_HPPC(k)=length(HPPC_current_extreme{k});
    elseif k >= length(HPPC_SOC_levels) -2
        length_HPPC(k)=length(HPPC_current_extreme{k-21});
    else
        length_HPPC(k)=length_HPPC_standard;
    end
end

startIndex_unfiltered = find(I == HPPC_current_standard(1));
startIndex=[startIndex_unfiltered(1)];
for k = 2:length(startIndex_unfiltered)
    if ~(startIndex_unfiltered(k)==startIndex_unfiltered(k-1)+1)
        startIndex=[startIndex,startIndex_unfiltered(k)];
    end
end
startIndex=startIndex-1;
%HPPC_voltage=zeros(length_HPPC,length(startIndex));
HPPC_voltage=cell(1,length(startIndex));
for k=1:length(startIndex)
    HPPC_voltage{k}=Volt(startIndex(k)+1:startIndex(k)+length_HPPC(k));
end
HPPC_voltage=HPPC_voltage';
HPPC_startingOCV=Volt(startIndex)';
HPPC_startingSOC=SOC(startIndex);
%%

%load('HPPC_data.mat');

figure()
clf

fs = 16;
color_array = {'r', 'g', 'b', 'm', 'c', 'k'}; % colors for lines

for i = 1:length(HPPC_voltage)
    plot(1:length(HPPC_voltage{i}), HPPC_voltage{i}, 'Color', color_array{mod(i-1,length(color_array))+1}, 'LineWidth',2); 
    hold on;
end

ylabel('Voltage [V]','FontSize', fs)
xlabel('Time [sec]','FontSize', fs)
set(gca,'FontSize', fs)
%legend({'HPPC 1', 'HPPC 2', 'HPPC 3', 'HPPC 4', 'HPPC 5', 'HPPC 6'},'interpreter','latex') % adjust as necessary
xlim([1, 1000])

hold off;

figure(5)
clf

voltage_adjusted=cell(1, length(startIndex));
for i = 1:length(HPPC_voltage)
    voltage_adjusted{i} = HPPC_voltage{i} - HPPC_startingOCV(i); % subtract starting voltage
    plot(1:length(HPPC_voltage{i}), voltage_adjusted{i}, 'Color', color_array{mod(i-1,length(color_array))+1}, 'LineWidth',2); 
    hold on;
end

ylabel('Voltage [V]','FontSize', fs)
xlabel('Time [sec]','FontSize', fs)
set(gca,'FontSize', fs)
xlim([1, 1000])
% average_voltage = mean(voltage_adjusted);  
% plot(1:length(average_voltage), average_voltage, 'k--', 'LineWidth',3);
%legend({'HPPC 1', 'HPPC 2', 'HPPC 3', 'HPPC 4', 'HPPC 5', 'HPPC 6','HPPC average'},'interpreter','latex') % adjust as necessary
hold off;

%save('HPPC_average_data.mat', 'average_voltage');
save('HPPC_voltage_adjusted.mat', 'voltage_adjusted', 'HPPC_startingSOC','HPPC_startingOCV','HPPC_current');
