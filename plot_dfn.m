%% Plot Doyle-Fuller-Newman Model Results
%   Created May 23, 2012 by Scott Moura
clc
clear all
close all;
fs =20;
lw=2;
run param/params_LCO
%run params_NMC_Samsung_new_iteration
[cn_low,cp_low] = init_cs(p,p.volt_min);
[cn_high,cp_high] = init_cs(p,p.volt_max);
Delta_cn = cn_high-cn_low;
Delta_cp = cp_low-cp_high;
p.OneC = min(p.epsilon_s_n*p.L_n*Delta_cn*p.Faraday/3600, p.epsilon_s_p*p.L_p*Delta_cp*p.Faraday/3600);
OneC=p.OneC;
%% Merge .mat files
mat_files = dir('results_*.mat');
all_T = [];
all_SOC = [];
all_I = [];
all_Volt = [];
all_c_ss_n = [];
all_c_ss_p = [];
all_c_e_0n = [];
all_c_e_0p = [];

for i = 1:length(mat_files)
    load(mat_files(i).name, 'T_save', 'SOC_save', 'I_save', 'Volt_save', 'c_ss_n_save', 'c_ss_p_save', 'c_e_0n_save', 'c_e_0p_save');
    if i==1
        T_save=T_save';
        SOC_save=SOC_save';
        Volt_save=Volt_save';
        c_e_0n_save=c_e_0n_save';
        c_e_0p_save=c_e_0p_save';

    end
    all_T = [all_T T_save];
    all_SOC = [all_SOC SOC_save];
    if iscolumn(I_save)
        I_save=I_save';
    end
    all_I = [all_I, I_save];
    all_Volt = [all_Volt Volt_save];
    all_c_ss_n = [all_c_ss_n, c_ss_n_save]; % append along columns because it's a 2-row matrix
    all_c_ss_p = [all_c_ss_p, c_ss_p_save]; % append along columns because it's a 2-row matrix
    all_c_e_0n = [all_c_e_0n, c_e_0n_save];
    all_c_e_0p = [all_c_e_0p, c_e_0p_save];
end

t=[0:length(all_I)-1];
all_I_C=all_I/p.OneC;
%% Histogram of current values

histogram(all_I_C, 'BinEdges', -4:0.1:4);

% Add labels and title for clarity
xlabel('Current (C rate)');
ylabel('Frequency');
title('Histogram of Current Values in C Rates');
%% extract rest voltages at 5% SOC increments
volt_switches_discharge = [];% Initialize empty array to store voltage values at switch points
volt_switches_charge = [];
SOC_switches_discharge = [];
SOC_switches_charge = [];
for i = 260:length(all_I)  % Start where SOC 5% cycling happens
    if all_I(i-1) == 0 && all_I(i) < 0  % If the current switches from 0 to negative (charge)
        volt_switches_charge = [volt_switches_charge, all_Volt(i-1)];  % Record the corresponding voltage
        SOC_switches_charge = [SOC_switches_charge, all_SOC(i-1)];
    end
    if all_I(i-1) == 0 && all_I(i) > 0  % If the current switches from 0 to positve (discharge)
        volt_switches_discharge = [volt_switches_discharge, all_Volt(i-1)];  % Record the corresponding voltage
        SOC_switches_discharge = [SOC_switches_discharge, all_SOC(i-1)];
    end
end

volt_switches_charge=[volt_switches_charge, all_Volt(end)]; %add the last value, not captured by loop
SOC_switches_charge=[SOC_switches_charge, all_SOC(end)];
figure()
plot(SOC_switches_charge,volt_switches_charge,'o')
hold on
plot(SOC_switches_discharge,volt_switches_discharge,'o')


%% Plot Current, SOC, Voltage
figure()
clf
lw=3;


ax1=subplot(4,1,1);
plot(t,all_I/p.OneC,'LineWidth',lw)
hold on
ylabel('Current [C-rate]','FontSize', fs,'interpreter', 'latex')
set(gca,'FontSize', fs)
%legend({'$$I(t)$$'},'interpreter','latex')
xlim([t(1), t(end)])
set(gca, 'XTickLabel', {}); 
xlabel(''); 
set(gca, 'Position', [0.1300, 0.77, 0.7750, 0.2157])
ylim([-2.1,2.1])


ax2=subplot(4,1,2);
plot(t,all_SOC,'LineWidth',lw)
ylabel('SOC','FontSize',fs,'interpreter','latex')
set(gca,'FontSize', fs)
%legend({'$$SOC(t)$$'},'interpreter','latex')
xlim([t(1), t(end)])
set(gca, 'XTickLabel', {}); 
xlabel(''); 
set(gca, 'Position', [0.1300, 0.54, 0.7750, 0.2157])
ylim([-0.1,1.1])

ax3=subplot(4,1,3);
plot(t,all_Volt,'LineWidth',lw)
%yyaxis left
%plot(1:length(all_Volt(51100:51500)),all_Volt(51100:51500),'LineWidth',2)
ylabel('Voltage [V]','FontSize', fs,'interpreter','latex')
xlabel('Time [s]','FontSize', fs,'interpreter','latex')
%yyaxis right
%plot(1:length(all_Volt(51100:51500)),all_I(51100:51500)/p.OneC,'LineWidth',2,'LineStyle',':')
%ylabel('Current [C-Rate]','FontSize', fs,'interpreter','latex')
%legend({'$$V(t)$$'},'interpreter','latex')
set(gca,'FontSize', fs)
xlim([1,400])
xlim([t(1), t(end)])
set(gca, 'XTickLabel', {}); 
xlabel(''); 
set(gca, 'Position', [0.1300, 0.31, 0.7750, 0.2157])
yticks([3.2,4,4.8])


ax4=subplot(4,1,4);
plot(t,all_T-273.15,'LineWidth',lw)
ylabel('Temperature [$^\circ$C]','FontSize', fs,'interpreter','latex')
xlabel('Time [s]','FontSize', fs,'interpreter','latex')
%legend({'$$T(t)$$'},'interpreter','latex')
set(gca,'FontSize', fs)
xlim([t(1), t(end)])
set(gca, 'Position', [0.1300, 0.08, 0.7750, 0.2157])

linkaxes([ax1 ax2 ax3 ax4],'x')
%% Plot Current, SOC, Voltage version 2
close all
figure()
clf
include_voltage_plot = true;
include_current_plot = true;
include_temperature_plot = false;
include_SOC_plot = true;
tiles=tiledlayout(4,1); 
start=1;


deltaI = diff(all_I);
indices = find(deltaI);
current_at_indices=all_I(indices);
combined=[indices',current_at_indices'];
v_indices=combined(combined(:,2)==0,1);
v_indices(end+1)=length(all_Volt);


if include_voltage_plot
    ax1=nexttile(1, [2 1]);
    plot(t(start:end)/3600, all_Volt(start:end), 'LineWidth', 2,'Color','b'); 
    hold on
    scatter(t(v_indices)/3600, all_Volt(v_indices), 'LineWidth', 2,'Color','r'); 
    ylabel('Voltage (V)', 'FontSize', fs, 'interpreter', 'latex');
    %xlabel('Time [h]', 'FontSize', fs, 'interpreter', 'latex');
    ax=gca;
    ax.FontSize = fs;
    ylim([3.65,4.15])
end

if include_current_plot
    ax2=nexttile(3);
    plot(t(start:end)/3600, all_I(start:end)/OneC, 'LineWidth', 2);
    ylabel('C-rate', 'FontSize', fs, 'interpreter', 'latex');
    %xlabel('Time [h]', 'FontSize', fs, 'interpreter', 'latex');
    ax=gca;
    ax.FontSize = fs;
end

if include_SOC_plot
    ax3=nexttile(4);
    plot(t(start:end)/3600, all_SOC(start:end), 'LineWidth', 2);
    ylabel('SOC', 'FontSize', fs, 'interpreter', 'latex');
    xlabel('Time (h)', 'FontSize', fs, 'interpreter', 'latex');
    ax=gca;
    ax.FontSize = fs;
end


ax = findobj(gcf,'Type','Axes');
linkaxes(ax, 'x');
xlim([4,6])
xticklabels(ax1,{})
xticklabels(ax2,{})
tiles.TileSpacing = 'compact';
ax1.FontSize = fs;
ax2.FontSize = fs;
ax3.FontSize = fs;
set(1, 'units', 'centimeters', 'pos', [0 0 35 20])
box on
exportgraphics(tiles,'OCV Test.pdf')
%% Plot R0, tau approx scheme
close all
figure
tiles=tiledlayout(3,1); 
start=1;
ax1=nexttile(1, [2 1]);
hold on
start_CC0=14640;
start_CC=14641;
start_R0=14712;
start_tau=14713;
plot(t(start:end)/3600, all_Volt(start:end), 'LineWidth', 2); 
% plot(t(start:start_CC0)/3600, all_Volt(start:start_CC0), 'LineWidth', 2,'Color','k'); 
% plot(t(start_CC0:start_CC)/3600, all_Volt(start_CC0:start_CC), 'LineWidth', 2,'Color','r','LineStyle','--'); 
% plot(t(start_CC:start_R0)/3600, all_Volt(start_CC:start_R0), 'LineWidth', 2,'Color','b','LineStyle',':'); 
% plot(t(start_R0:start_tau)/3600, all_Volt(start_R0:start_tau), 'LineWidth', 2,'Color','r','LineStyle','--'); 
% plot(t(start_tau:end)/3600, all_Volt(start_tau:end), 'LineWidth', 2,'Color','g','LineStyle','-.'); 
xlim([4.05,4.16])
xlabel('Time (min)', 'FontSize', fs, 'interpreter', 'latex');
ylabel('Voltage (V)', 'FontSize', fs, 'interpreter', 'latex');
%xlabel('Time [min]', 'FontSize', fs, 'interpreter', 'latex');
ax=gca;
xlim([4.05,4.16])
ax.FontSize = fs;
box on
xticklabels(ax1,{})
%legend({'Rest Voltage','','Influenced by: OCV, $\tau_1$','Influenced by: $R_0$','Influenced by: $\tau_1$'},"Interpreter","latex",'Location','southeast')

textbox = annotation("textbox",'interpreter','latex','FontSize',fs);
textbox.String='$\approx 4 R_1 \cdot C_1 = 4\tau_1$';
textbox.String='$\Delta v_\infty = (R_0 + R_1) \Delta i $';
textbox.String='$\Delta v_0 = R_0\Delta i$';
ax2=nexttile(3);
plot(t(start:end)/3600, all_I(start:end)/OneC, 'LineWidth', 2); 
ylabel('C-rate', 'FontSize', fs, 'interpreter', 'latex');
xlabel('Time (min)', 'FontSize', fs, 'interpreter', 'latex');
ax=gca;
xlim([4.05,4.16])
box on
ax.FontSize = fs;
set(gca,'XTickLabel',1:1:6)
% annotation('textarrow','String', '$\tau_1$','Interpreter', 'latex','FontSize', fs); 
% textbox = annotation("textbox",'interpreter','latex','FontSize',fs);
% textbox.String='\textbf{1}';

set(1, 'units', 'centimeters', 'pos', [0 0 35 15])
exportgraphics(tiles,'ECM visual3.pdf')
%exportgraphics(tiles,'Initial_condition_estimation_V2.pdf')

%% Plot Current, c_ss(x,t) & c_e(0^\pm,t)
figure()
clf

ax5=subplot(3,1,1);
plot(t,all_I/p.OneC,'LineWidth',2)
ylabel('Current [C-rate, h^{-1}]','FontSize', fs)
set(gca,'FontSize', fs)
legend({'$$I(t)$$'},'interpreter','latex')
xlim([t(1), t(end)])

ax6=subplot(3,1,2);
plot(t,all_c_ss_n(1,:)/p.c_s_n_max,'b-','LineWidth',2); hold on;
plot(t,all_c_ss_n(end,:)/p.c_s_n_max,'b--','LineWidth',2); hold on;
plot(t,all_c_ss_p(1,:)/p.c_s_p_max,'r--','LineWidth',2); hold on;
plot(t,all_c_ss_p(end,:)/p.c_s_p_max,'r-','LineWidth',2); hold on;
ylabel('Solid Conc. [Stoic, -]','FontSize',fs)
set(gca,'FontSize', fs)
legend({'$$c_{ss}^-(0^-,t)$$','$$c_{ss}^-(L^-,t)$$',...
    '$$c_{ss}^+(L^+,t)$$','$$c_{ss}^+(0^+,t)$$'},...
    'interpreter','latex','Location','Best')
xlim([t(1), t(end)])

ax7=subplot(3,1,3);
plot(t,all_c_e_0n,'b-',t,all_c_e_0p,'r-','LineWidth',2)
ylabel('Elec Conc. [mol/m^3]','FontSize', fs)
xlabel('Time [sec]','FontSize', fs)
legend({'$$c_e(0^-,t)$$','$$c_e(0^+,t)$$'},'interpreter','latex')
set(gca,'FontSize', fs)
xlim([t(1), t(end)])

linkaxes([ax5 ax6 ax7],'x')

%% Plot voltage-SOC hysteresis

figure()
clf
% For charge - blue color
plot(SOC_switches_charge, volt_switches_charge, 'x-', 'Color', [0 0 1 0.5], 'LineWidth', 1.5);
hold on
% For discharge - red color
plot(SOC_switches_discharge, volt_switches_discharge, 'x-', 'Color', [1 0 0 0.5], 'LineWidth', 1.5);
hold off

ylabel('OCV at rest (after 1h)','FontSize' , fs)
xlabel('SOC','FontSize', fs)
legend({'Charge','Discharge'},'interpreter','latex')
set(gca,'FontSize', fs)

%% Plot a normal HPPC cycle and an extreme HPPC cycle

close all
HPPC1=[59582:60383];
HPPC2=[126496:127447];
figure()
fs=24;
subplot(3,1,1)
plot(t(3000:length(all_Volt))/3600,all_Volt(3000:end),'LineWidth',lw,'Color','k')
hold on
plot(t(HPPC1)/3600,all_Volt(HPPC1),'LineWidth',lw,'Color','r')
plot(t(HPPC2)/3600,all_Volt(HPPC2),'LineWidth',lw,'Color','b')
set(gca,'FontSize', fs)
ylabel('Voltage (V)','interpreter','latex','FontSize',fs)
xlabel('Time (h)','interpreter','latex','FontSize',fs)
xlim([t(3000)/3600, length(t)/3600])
% yyaxis right
% set(gca,'YColor','#77AC30');
% ylabel('SOC','interpreter','latex','FontSize',fs)
% plot(t(3000:length(all_Volt))/3600,all_SOC(3000:end),'LineWidth',lw,'Color','#77AC30','LineStyle',':')

%set(gca, 'Position', [0.1300, 0.7093, 0.7750, 0.2157]) 
xLimits = xlim;
set(gca, 'XTick', [0:10:35]);
%set(gca, 'XTickLabel', {num2str(0), num2str(xLimits(2))});


ax1=subplot(3,1,2);
plot(t(1:length(HPPC1)),all_Volt(HPPC1)-all_Volt(HPPC1(1)),'LineWidth',lw,'Color','r')
set(gca,'FontSize', fs)
ylabel('Polarisation voltage (V)','interpreter','latex','FontSize',fs)
hold on
plot(t(1:length(HPPC2)),all_Volt(HPPC2)-all_Volt(HPPC2(1)),'LineWidth',lw,'Color','b','LineStyle','-.')
xlim([t(1), length(HPPC2)])
xlabel('Time (s)','interpreter','latex','FontSize',fs)
%set(gca, 'Position', [0.1300, 0.37, 0.7750, 0.2157]) 
set(gca, 'XTick', [0:200:800]);


ax2=subplot(3,1,3);
plot(t(1:length(HPPC1)),all_I(HPPC1)/p.OneC,'LineWidth',lw,'Color','r','LineStyle','-')
set(gca,'FontSize', fs)
ylabel('C-rate','interpreter','latex','FontSize',fs)
xlabel('Time (s)','interpreter','latex','FontSize',fs)
hold on
plot(t(551:length(HPPC2)),all_I(HPPC2(551:end))/p.OneC,'LineWidth',lw,'Color','b','LineStyle','-.')
plot(t(1:length(HPPC2(1:550))),all_I(HPPC2(1:550))/p.OneC,'LineWidth',2,'Color','b','LineStyle','--')
xlim([t(1), length(HPPC2)])
legend({'standard HPPC', 'asymmetric HPPC'},'interpreter','latex','Location','northwest')

%legend({'standard HPPC', 'extreme HPPC'},'interpreter','latex','Location','northwest')
%set(gca, 'Position', [0.1300, 0.1100, 0.7750, 0.2157]) 
set(gca, 'XTick', [0:200:800]);


linkaxes([ax1 ax2],'x')
set(1, 'units', 'centimeters', 'pos', [0 0 35 20])
box on
tiles=figure(1);
exportgraphics(tiles,'HPPC combo4.pdf')

%% plot full HPPC


%% plot CC CV
fs=20;
close all
figure()
tiles=tiledlayout(4,1); 
CCCVstart=22000;
CCCVend=27500;
CCCV=[CCCVstart:CCCVend];
CCstart=22469;
CVstart=25638;
relaxstart=26709;
nexttile(1, [2 1]);

%ax1=subplot(2,1,1);
plot(t(1:CCstart-CCCVstart+1)/3600,all_Volt(CCCVstart:CCstart),'LineWidth',lw,'Color','k')
hold on
plot(t(CCstart+1-CCCVstart:CVstart-CCCVstart)/3600,all_Volt(CCstart:CVstart-1),'LineWidth',lw,'Color','r','LineStyle','-.')
plot(t(CVstart+1-CCCVstart:relaxstart-CCCVstart)/3600,all_Volt(CVstart:relaxstart-1),'LineWidth',lw,'Color','b','LineStyle',':')
plot(t(relaxstart-CCCVstart:CCCVend-CCCVstart)/3600,all_Volt(relaxstart-1:CCCVend-1),'LineWidth',lw,'Color','k')

ylabel('Voltage (V)','interpreter','latex','FontSize',fs)
%xlabel('Time [h]','interpreter','latex','FontSize',fs)
legend({'','Constant current phase','Constant voltage phase',''},'interpreter','latex','Location','northwest')
ax1=gca;
set(gca,'FontSize', fs)
%ax.XAxis.Exponent = 3;
ylim([3.4,4.6])
xlim([t(1), length(CCCV)/3600])



nexttile(3);
plot(t(1:CCstart-CCCVstart+1)/3600,all_I(CCCVstart:CCstart)/OneC,'LineWidth',lw,'Color','k')
hold on
plot(t(CCstart+1-CCCVstart:CVstart-CCCVstart)/3600,all_I(CCstart:CVstart-1)/OneC,'LineWidth',lw,'Color','r','LineStyle','-.')
plot(t(CVstart+1-CCCVstart:relaxstart-CCCVstart)/3600,all_I(CVstart:relaxstart-1)/OneC,'LineWidth',lw,'Color','b','LineStyle',':')
plot(t(relaxstart-CCCVstart:CCCVend-CCCVstart)/3600,all_I(relaxstart-1:CCCVend-1)/OneC,'LineWidth',lw,'Color','k')

ylabel('C-rate', 'FontSize', fs, 'interpreter', 'latex');
%xlabel('Time [h]', 'FontSize', fs, 'interpreter', 'latex');
ax2=gca;
set(gca,'FontSize', fs)
%ax.XAxis.Exponent = 3;
xlim([t(1), length(CCCV)/3600])



% nexttile(5);
% plot(t(1:CCstart-CCCVstart+1)/3600,all_T(CCCVstart:CCstart)-273.15,'LineWidth',lw,'Color','k')
% hold on
% plot(t(CCstart+1-CCCVstart:CVstart-CCCVstart)/3600,all_T(CCstart:CVstart-1)-273.15,'LineWidth',lw,'Color','r','LineStyle','-.')
% plot(t(CVstart+1-CCCVstart:relaxstart-CCCVstart)/3600,all_T(CVstart:relaxstart-1)-273.15,'LineWidth',lw,'Color','b','LineStyle',':')
% plot(t(relaxstart-CCCVstart:CCCVend-CCCVstart)/3600,all_T(relaxstart-1:CCCVend-1)-273.15,'LineWidth',lw,'Color','k')
% 
% ylabel('Temp. [$^\circ$C]', 'FontSize', fs, 'interpreter', 'latex');
% ylim([25,50])
% %xlabel('Time [h]', 'FontSize', fs, 'interpreter', 'latex');
% ax3=gca;
% set(gca,'FontSize', fs)
% %ax.XAxis.Exponent = 3;
% xlim([t(1), length(CCCV)/3600])



nexttile(4);
plot(t(1:CCstart-CCCVstart+1)/3600,all_SOC(CCCVstart:CCstart),'LineWidth',lw,'Color','k')
hold on
plot(t(CCstart+1-CCCVstart:CVstart-CCCVstart)/3600,all_SOC(CCstart:CVstart-1),'LineWidth',lw,'Color','r','LineStyle','-.')
plot(t(CVstart+1-CCCVstart:relaxstart-CCCVstart)/3600,all_SOC(CVstart:relaxstart-1),'LineWidth',lw,'Color','b','LineStyle',':')
plot(t(relaxstart-CCCVstart:CCCVend-CCCVstart)/3600,all_SOC(relaxstart-1:CCCVend-1),'LineWidth',lw,'Color','k')

ylabel('SOC', 'FontSize', fs, 'interpreter', 'latex');
xlabel('Time (h)', 'FontSize', fs, 'interpreter', 'latex');
ax4=gca;
set(gca,'FontSize', fs)
%ax.XAxis.Exponent = 3;
xlim([t(1), length(CCCV)/3600])


xticklabels(ax1,{})
xticklabels(ax2,{})
tiles.TileSpacing = 'compact';
ax1.FontSize = fs;
ax2.FontSize = fs;
ax4.FontSize = fs;
set(1, 'units', 'centimeters', 'pos', [0 0 35 20])
box on
exportgraphics(tiles,'CCCV_mod2.pdf')


%linkaxes([ax1 ax2],'x')


