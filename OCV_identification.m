clc
clear
close all
fs=20;
lw=2;
pybamm=false;


mat_files = dir('results_*.mat'); %Merge data
all_SOC = [];
all_I = [];
all_Volt= [];
for i = 1:length(mat_files)
    load(mat_files(i).name, 'SOC_save','Volt_save', 'I_save');
    if i==1
        SOC_save=SOC_save';
        Volt_save=Volt_save';
    end
    all_SOC = [all_SOC SOC_save];
    all_I = [all_I, I_save'];
    all_Volt = [all_Volt Volt_save];
end
if pybamm 
    load("sol_data.mat", 'I','OneC', 'SoC','T','V');
    t = 1:max(T);
    all_I = interp1(T, I, t, 'linear');
    all_SOC = interp1(T, SoC, t, 'linear');
    all_Volt = interp1(T, V, t, 'linear');

    figure;
    subplot(3,1,1);
    plot(t, all_I);
    title('Interpolated Current (I) vs. Time');
    xlabel('Time');
    ylabel('Current (I)');
    
    subplot(3,1,2);
    plot(t, all_SOC);
    title('Interpolated State of Charge (SoC) vs. Time');
    xlabel('Time');
    ylabel('State of Charge (SoC)');
    
    subplot(3,1,3);
    plot(t, all_Volt);
    title('Interpolated Voltage (V) vs. Time');
    xlabel('Time');
    ylabel('Voltage (V)');
end

OCV_start=1800;
OCV_end=length(all_I);
SOC=all_SOC(OCV_start:OCV_end);
Volt=all_Volt(OCV_start:OCV_end);
I=all_I(OCV_start:OCV_end);
t=0:length(Volt)-1;

figure()

plot(t,Volt,'LineWidth', lw*1,'LineStyle','-')
%xlim([10000,15000])
ylabel('Voltage [V]', 'FontSize', fs, 'interpreter', 'latex');
xlabel('time', 'FontSize', fs, 'interpreter', 'latex');
set(gca,'FontSize', fs)

%% extract individual cycles
clc
deltaI = diff(I);
indices = find(deltaI);
current_at_indices=I(indices);
combined=[indices',current_at_indices'];
v_indices=combined(combined(:,2)==0,1);
v_indices(end+1)=length(Volt);
%odd_indices=[indices(1:2:end), length(I)]; %add the last position of I, too

OCV_SOC=zeros(2,length(v_indices));
OCV_SOC(1,:)=Volt(v_indices);

OCV_SOC(2,:)=SOC(v_indices);

figure()
plot(OCV_SOC(2,:),OCV_SOC(1,:),'LineWidth', lw*4,'LineStyle','-')

ylabel('OCV [V]', 'FontSize', fs, 'interpreter', 'latex');
xlabel('SOC', 'FontSize', fs, 'interpreter', 'latex');
set(gca,'FontSize', fs)
xlim([-0.01,1])
ylim([3.1,4.5])
plot(OCV_SOC(2,:),OCV_SOC(1,:),'o')

%save('OCV_SOC_Marquis.mat','OCV_SOC')
%% something
V_ECM_collection = cell(ceil(length(indices)/2), 1);
V_ECM_collection_norm = cell(ceil(length(indices)/2), 1);
V_DFN_collection = cell(ceil(length(indices)/2), 1);
V_DFN_collection_norm = cell(ceil(length(indices)/2), 1);
I_collection = cell(ceil(length(indices)/2), 1);
start_OCV=V_RCfit(odd_indices);

for k = 1:2:length(indices)
    if k+2 <= length(indices)
        V_ECM_collection{(k+1)/2} = V_ECM(indices(k):indices(k+2));
        V_DFN_collection{(k+1)/2} = V_RCfit(indices(k):indices(k+2));
        I_collection{(k+1)/2} = I_RCfit(indices(k):indices(k+2));
    else
        V_ECM_collection{(k+1)/2} = V_ECM(indices(k):end);
        V_DFN_collection{(k+1)/2} = V_RCfit(indices(k):end);
        I_collection{(k+1)/2} = I_RCfit(indices(k):end);
    end
    V_ECM_collection_norm{(k+1)/2} = V_ECM_collection{(k+1)/2}-start_OCV((k+1)/2);
    V_DFN_collection_norm{(k+1)/2} = V_DFN_collection{(k+1)/2}-start_OCV((k+1)/2);
end

%% plot 
clc
figure()
hold on
for k = 1:length(I_collection)
    plot(1:length(V_ECM_collection{k}),V_ECM_collection{k},'LineWidth',2)
end

%% Normalize

clc
figure()
hold on
for k = 1:length(I_collection)
    plot(1:length(V_ECM_collection_norm{k}),V_ECM_collection_norm{k},'LineWidth',2)
end
